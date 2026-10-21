#include <DislocationNetwork.h>
#include <DefectiveCrystal/DefectiveCrystal.h>
#include <PeriodicVelocityConstraints.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <map>
#include <iterator>

void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}


using Network=model::DislocationNetwork<3>;
using Triplets=std::deque<Eigen::Triplet<double>>;
using Clock=std::chrono::steady_clock;

Eigen::VectorXd referenceVelocity(Network& network,const Triplets& triplets,const Eigen::VectorXd& F)
{
    Triplets z;
    const auto freeCount=model::periodicVelocityConstraints<3>(network,z);
    Eigen::SparseMatrix<double> K(F.size(),F.size()),Z(F.size(),3*freeCount);
    K.setFromTriplets(triplets.begin(),triplets.end());
    Z.setFromTriplets(z.begin(),z.end());
    const Eigen::SparseMatrix<double> reduced(Z.transpose()*K*Z);
    const Eigen::VectorXd force(Z.transpose()*F);
    Eigen::VectorXd diagonal=Eigen::VectorXd::Zero(Z.cols());
    for(int k=0;k<reduced.outerSize();++k)
    {
        for(Eigen::SparseMatrix<double>::InnerIterator it(reduced,k);it;++it)
        {
            if(it.row()==it.col())
            {
                diagonal(it.row())+=it.value();
            }
            else
            {
                diagonal(it.row())+=0.5*it.value();
                diagonal(it.col())+=0.5*it.value();
            }
        }
    }
    Eigen::VectorXd velocity=Eigen::VectorXd::Zero(Z.cols());
    for(Eigen::Index k=0;k<diagonal.size();++k)
    {
        if(std::abs(diagonal(k))>FLT_EPSILON)
        {
            velocity(k)=force(k)/diagonal(k);
        }
        else
        {
            require(std::abs(force(k))<=FLT_EPSILON,"Force on a zero-mass degree of freedom");
        }
    }
    return Z*velocity;
}

void verify(Network& network)
{
    network.updateGeometry();
    const auto nodeCount=network.networkNodes().size();
    const auto loopCount=network.loops().size();
    std::map<size_t,std::pair<Eigen::Vector3d,Eigen::Vector3d>> positions;
    for(const auto& entry : network.loopNodes())
    {
        const auto ln=entry.second.lock();
        positions.emplace(ln->sID,std::make_pair(ln->get_P(),ln->networkNode->get_P()));
    }
    std::map<std::pair<size_t,size_t>,std::pair<Eigen::Vector3d,std::set<size_t>>> topology;
    for(const auto& entry : network.networkLinks())
    {
        const auto link=entry.second.lock();
        std::set<size_t> loops;
        for(const auto* ll : link->loopLinks())
        {
            loops.insert(ll->loop->sID);
        }
        topology.emplace(entry.first,std::make_pair(link->burgers(),loops));
        link->createQuadraturePoints(false);
    }
    const auto capture=[&]()
    {
        std::vector<Eigen::Matrix<double,18,1>> values;
        for(const auto& entry : network.networkLinks())
        {
            for(const auto& qp : entry.second.lock()->quadraturePoints())
            {
                Eigen::Matrix<double,18,1> value;
                value.head<9>()=Eigen::Map<const Eigen::Matrix<double,9,1>>(qp.stress.data());
                value.segment<3>(9)=qp.pkForce;
                value.segment<3>(12)=qp.lineTensionForce;
                value.tail<3>()=qp.velocity;
                values.push_back(value);
            }
        }
        return values;
    };
    const auto t0=Clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t k=0;k<network.networkLinks().size();++k)
    {
        auto it=network.networkLinks().begin();
        std::advance(it,k);
        it->second.lock()->updateQuadraturePoints(false);
    }
    const double oldFields=std::chrono::duration<double>(Clock::now()-t0).count();
    const auto oldValues=capture();
    require(!oldValues.empty(),"Fixture has no quadrature points");
    Triplets oldK;
    Eigen::VectorXd oldF=Eigen::VectorXd::Zero(3*nodeCount);
    const auto t1=Clock::now();
    for(const auto& entry : network.networkLinks())
    {
        const auto link=entry.second.lock();
        if(!link->isBoundarySegment())
        {
            link->addToGlobalAssembly(oldK,oldF);
        }
    }
    const double oldAssembly=std::chrono::duration<double>(Clock::now()-t1).count();
    const Eigen::VectorXd oldVelocity=referenceVelocity(network,oldK,oldF);

    // Start both traversal methods from identical, freshly created quadrature points.
    for(const auto& entry : network.networkLinks())
    {
        entry.second.lock()->createQuadraturePoints(false);
    }
    const auto t2=Clock::now();
    std::vector<std::shared_ptr<Network::NetworkLinkType>> links;
    links.reserve(network.networkLinks().size());
    for(const auto& entry : network.networkLinks())
    {
        links.push_back(entry.second.lock());
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t k=0;k<links.size();++k)
    {
        links[k]->updateQuadraturePoints(false);
    }
    const double newFields=std::chrono::duration<double>(Clock::now()-t2).count();
    const auto newValues=capture();
    require(newValues.size()==oldValues.size(),"Snapshot changed the quadrature count");
    double fieldError=0;
    for(size_t i=0;i<oldValues.size();++i)
    {
        require(newValues[i].allFinite(),"Nonfinite quadrature result");
        fieldError=std::max(fieldError,(newValues[i]-oldValues[i]).cwiseAbs().maxCoeff());
    }
    require(fieldError==0,"Snapshot changed a quadrature field, force or velocity");
    const auto t3=Clock::now();
    std::map<size_t,size_t> indices;
    for(const auto& entry : network.networkNodes())
    {
        indices.emplace(entry.second.lock()->sID,indices.size());
    }
    Triplets newK;
    Eigen::VectorXd newF=Eigen::VectorXd::Zero(3*nodeCount);
    for(const auto& link : links)
    {
        if(!link->isBoundarySegment())
        {
            link->addToGlobalAssembly(newK,newF,indices.at(link->source->sID),indices.at(link->sink->sID));
        }
    }
    const double newAssembly=std::chrono::duration<double>(Clock::now()-t3).count();
    require(oldK.size()==newK.size(),"Indexed assembly changed matrix sparsity");
    for(size_t k=0;k<oldK.size();++k)
    {
        require(oldK[k].row()==newK[k].row() && oldK[k].col()==newK[k].col() && oldK[k].value()==newK[k].value(),
                "Indexed assembly changed a matrix entry");
    }
    require((oldF.array()==newF.array()).all(),"Indexed assembly changed the force vector");
    const Eigen::VectorXd velocity=network.glideSolver->getNodeVelocities();
    const double velocityError=(velocity-oldVelocity).cwiseAbs().maxCoeff();
    require(velocityError==0,"Production solver differs from legacy assembly");
    require(network.networkNodes().size()==nodeCount && network.loops().size()==loopCount &&
            network.networkLinks().size()==topology.size(),"Field evaluation changed topology");
    for(const auto& entry : network.networkLinks())
    {
        const auto link=entry.second.lock();
        std::set<size_t> loops;
        for(const auto* ll : link->loopLinks())
        {
            loops.insert(ll->loop->sID);
        }
        require(loops==topology.at(entry.first).second &&
                (link->burgers().array()==topology.at(entry.first).first.array()).all(),"Junction sharing changed");
    }
    for(const auto& entry : network.loopNodes())
    {
        const auto ln=entry.second.lock();
        require((ln->get_P().array()==positions.at(ln->sID).first.array()).all() &&
                (ln->networkNode->get_P().array()==positions.at(ln->sID).second.array()).all(),
                "Field evaluation changed the periodic loop/network positions");
    }
    std::cout<<"\nCACHE BENCHMARK nodes="<<nodeCount<<" links="<<links.size()<<" qpoints="<<oldValues.size()
             <<" fields="<<oldFields<<" -> "<<newFields<<" assembly="<<oldAssembly<<" -> "<<newAssembly
             <<" fieldError="<<fieldError<<" velocityError="<<velocityError<<'\n';
}

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2 || argc == 3, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        if(argc==3 && std::filesystem::is_directory(argv[2]))
        {
#ifdef _OPENMP
            omp_set_num_threads(12);
#endif
            model::DefectiveCrystalBase<3> base(argv[2]);
            model::DefectiveCrystalConfigIO<3> saved(base.simulationParameters.traitsIO.configFolder);
            saved.read(base.simulationParameters.runID);
            model::DefectiveCrystal<3> crystal(base);
            crystal.initializeConfiguration(saved);
            const auto network=crystal.getUniqueTypedMicrostructure<Network>();
            require(bool(network) && !network->computeDDinteractions,"Benchmark requires DD with interactions disabled");
            verify(*network);
            return 0;
        }
        const auto folder = std::filesystem::temp_directory_path() /
            ("glide-index-caches-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");

        std::filesystem::copy_file(root / "Library/Meshes/bicrystal_48.msh", input / "mesh.msh");
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        std::filesystem::create_directories(folder / "config");
        const double length=3671.9477120460206;
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
            std::ofstream p(input / "polycrystal.txt");
            p << std::setprecision(17) << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\n"
              << "C2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF=" << length << " 0 0 0 " << length << " 0 0 0 " << length
              << ";\nX0=0 0 0;\nperiodicFaceIDs=-1;\n";
        }
        model::DefectiveCrystalBase<3> base(folder.string());
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        model::DefectiveCrystalConfigIO<3> saved(folder.string());
        std::ifstream stream(argc==3 ? std::filesystem::path(argv[2]) : root/"tests/data/bicrystalSharedConstraints1915.txt");
        require(bool(stream),"Missing bicrystal junction regression fixture");
        saved.readTxtStream(stream);
        network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
        network.computeDDinteractions=false;
#ifdef _OPENMP
        omp_set_num_threads(argc==3 ? 12 : 2);
#endif
        verify(network);
        if(argc==2)
        {
            std::shared_ptr<Network::NetworkLinkType> shared;
            for(const auto& entry : network.networkLinks())
            {
                const auto link=entry.second.lock();
                if(link->loopLinks().size()>1 && link->chordLength()>1 && !link->isBoundarySegment())
                {
                    shared=link;
                    break;
                }
            }
            require(bool(shared),"Missing shared junction for topology-change test");
            const auto memberships=shared->loopLinks().size();
            auto inserted=network.networkNodes().create(shared->get_r(0.5));
            network.expandNetworkLink(shared,inserted);
            shared.reset();
            require(inserted->loopNodes().size()==memberships,"Expansion lost junction sharing");
            verify(network);
            network.removeNetworkNode(inserted->sID);
            inserted.reset();
            verify(network);
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: exact cached-index equivalence before and after shared junction topology changes\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
