#include <FiniteElement.h>
#include <WeakProblem.h>
#include <GlobalRigidMotionPenalty.h>
#include <NonPeriodicExternalBoundary.h>
#include <ElasticDeformation.h>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <iostream>
using namespace model;
void require(bool v,const char* text)
{
    if(!v)
    {
        throw std::runtime_error(text);
    }
}
template<class Trial>
struct Load
{
    std::shared_ptr<Trial> trial;
    Eigen::VectorXd value;
    const auto& field() const
    {
        return trial->field();
    }
    Eigen::VectorXd globalVector() const
    {
        return value;
    }
};
void exercise(const std::filesystem::path& root,int rank)
{
    auto mesh=std::make_shared<SimplicialMesh<3>>((root/"Library/Meshes/bicrystal_48.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d(2,3,4),std::set<int>{-1});
    std::set<int> seeds;
    for(int d=0;d<rank;++d)
    {
        for(const auto& region : mesh->regions())
        {
            for(const auto& face : region.second->faces())
            {
                const auto& pair=face.second->periodicFacePair;
                if(pair.second && pair.first.cross(mesh->periodicBasis().col(d)).norm()<1.e-8)
                {
                    seeds.insert(face.first);
                    break;
                }
            }
        }
    }
    mesh->identifyParallelFaces(seeds);
    using FE=FiniteElement<LagrangeElement<3,2>>;
    auto fe=std::make_shared<FE>(mesh);
    auto u=fe->trial<3>();
    using Trial=typename decltype(u)::element_type;
    auto dV=fe->domain<EntireDomain,4,GaussLegendre>();
    Eigen::Matrix<double,6,6> C=Eigen::Matrix<double,6,6>::Identity();
    C.topLeftCorner<3,3>().array()+=1.0;
    auto strain=def(u);
    auto stiffness=(test(strain),C*strain)*dV;
    auto penalty=std::make_shared<GlobalRigidMotionPenalty<Trial>>(u,2.0);
    auto form=stiffness+penalty;
    const Eigen::MatrixXd F=penalty->lowRankPenaltyRows();
    const Eigen::MatrixXd R=penalty->rigidModes();
    const int modes=rank==0 ? 6 : (rank==1 ? 4 : 3);
    require(F.rows()==modes,"Wrong admissible rigid-mode count");
    Eigen::SparseMatrix<double> sparse(u->gSize(),u->gSize());
    auto triplets=stiffness->globalTriplets();
    sparse.setFromTriplets(triplets.begin(),triplets.end());
    const Eigen::MatrixXd K=sparse;
    require((K*R).norm()<1.e-10*K.norm()*R.norm(),"Rigid mode has elastic strain energy");
    require((F*R).fullPivLu().rank()==modes,"Penalty leaves a rigid mode free");
    // Independent explicit T, g reference, including a nonzero affine jump.
    std::map<size_t,size_t> columns;
    for(size_t n=0;n<fe->nodeSize();++n)
    {
        columns.emplace(fe->periodicRepresentative(n),columns.size());
    }
    Eigen::MatrixXd T=Eigen::MatrixXd::Zero(u->gSize(),3*columns.size());
    Eigen::VectorXd g=Eigen::VectorXd::Zero(u->gSize());
    Eigen::Matrix3d H=Eigen::Matrix3d::Zero();
    H(2,2)=0.02;
    u->setPeriodicGradient(H);
    for(size_t n=0;n<fe->nodeSize();++n)
    {
        T.block<3,3>(3*n,3*columns.at(fe->periodicRepresentative(n))).setIdentity();
        g.segment<3>(3*n)=H*fe->periodicShift(n);
    }
    Eigen::VectorXd q(T.cols());
    for(Eigen::Index i=0;i<q.size();++i)
    {
        q(i)=std::sin(0.31*i);
    }
    Eigen::VectorXd target=T*q+g;
    target-=R*(F*R).fullPivLu().solve(F*target);
    require((F*target).norm()<1.e-10,"Deformation has nonzero global rigid moment");
    require((K*target).norm()>1.e-4,"Test deformation is trivial");
    auto rhs=std::make_shared<Load<Trial>>(Load<Trial>{u,K*target});
    auto problem=weakProblem(form,rhs);
    problem->solveWithMasterSlaveConstraints();
    require((u->dofVector()-target).norm()<1.e-7*target.norm(),"Global penalty altered deformation");
    const Eigen::MatrixXd A=K+F.transpose()*F;
    const Eigen::MatrixXd reduced=T.transpose()*A*T;
    const Eigen::VectorXd b=T.transpose()*(rhs->value-A*g);
    const Eigen::VectorXd expected=T*reduced.ldlt().solve(b)+g;
    require((u->dofVector()-expected).norm()<1.e-7*expected.norm(),"Bordered solve differs from explicit global penalty");
    require((problem->reducedForceVector()-b).norm()<1.e-9*std::max(1.0,b.norm()),"Affine penalty lifting is wrong");
    require(problem->solverMatrix().rows()==T.cols()+modes,"Penalty was not kept low rank");
    rhs->value*=2.0;
    require(problem->solveWithMasterSlaveConstraints().reusedFactorization,"Changing RHS refactorized stiffness");
    auto strong=std::make_shared<GlobalRigidMotionPenalty<Trial>>(u,200.0);
    auto other=weakProblem(stiffness+strong,rhs);
    const Eigen::VectorXd solution=u->dofVector();
    other->solveWithMasterSlaveConstraints();
    require((u->dofVector()-solution).norm()<1.e-7*solution.norm(),"Penalty strength changed deformation");
    const auto external=fe->boundary<ExternalBoundary,3,GaussLegendre>();
    const auto physical=fe->boundary<NonPeriodicExternalBoundary,3,GaussLegendre>();
    for(const auto& face : physical)
    {
        require(!mesh->periodicSimplexPairs().count(&face.first->simplex.child(face.second)),"Periodic face received free-surface traction");
    }
    require((rank==0 && physical.size()==external.size()) || (rank>0 && physical.size()<external.size()),"Physical traction selection is wrong");
    // Fixing one point removes translations but retains compatible rotations about it.
    for(size_t i=0;i<3;++i)
    {
        u->dirichletConditions()[i]=solution(i);
    }
    require(penalty->rigidModes().cols()==(rank==0 ? 3 : (rank==1 ? 1 : 0)),"Point constraints removed the wrong rigid modes");
    // Homogeneous Dirichlet restrictions must also remove modes from the gauge.
    for(size_t i=0;i<u->gSize();++i)
    {
        u->dirichletConditions()[i]=solution(i);
    }
    require(penalty->lowRankPenaltyRows().rows()==0,"Fixed field still penalizes rigid modes");
    problem->solveWithMasterSlaveConstraints();
    require((u->dofVector()-solution).norm()<1.e-10,"Fixed solution changed");
}


struct UniformSource : ElasticDeformation<3>
{
    Eigen::Matrix3d plastic=Eigen::Matrix3d::Zero();
    using ElasticDeformation<3>::ElasticDeformation;
    Eigen::Matrix3d averagePlasticDistortion() const override
    {
        return plastic;
    }
};
void physicsIntegration(const std::filesystem::path& root,int rank,bool rotated,int mode,bool withSource=false,bool withPin=false,bool irregular=false)
{
    using V=Eigen::Matrix<double,6,1>;
    using M=Eigen::Matrix<double,6,6>;
    const std::string meshName=irregular ? "polycrystalCube_10_periodicXY.msh" : "bicrystal_48.msh";
    const Eigen::Matrix3d Q=rotated ? Eigen::AngleAxisd(0.63,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix() : Eigen::Matrix3d::Identity();
    StaticID<Polyhedron<3,3-1>>::force_count(0);
    std::set<int> seeds;
    {
        SimplicialMesh<3> mesh((root/"Library/Meshes"/meshName).string(),Q,Eigen::Vector3d::Zero(),{-1});
        for(const auto& region : mesh.regions())
        {
            for(const auto& face : region.second->faces())
            {
                const auto& pair=face.second->periodicFacePair;
                for(int i=0;i<rank;++i)
                {
                    const int axis=i==0 ? 1 : (i==1 ? 0 : 2);
                    if(pair.second && pair.first.normalized().cross(Q.col(axis)).norm()<1.e-8)
                    {
                        seeds.insert(face.first);
                    }
                }
            }
        }
    }
    StaticID<Polyhedron<3,3-1>>::force_count(0);
    const auto folder=std::filesystem::temp_directory_path()/("elastic-mixed-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    const auto input=folder/"inputFiles";
    std::filesystem::create_directories(input);
    std::filesystem::create_directories(folder/"config");
    std::filesystem::create_directories(folder/"F");
    std::filesystem::copy_file(root/"Library/Materials/W.txt",input/"W.txt");
    std::filesystem::copy_file(root/"Library/Meshes"/meshName,input/"mesh.msh");
    V strain0,strainRate,stress0,stressRate,alpha;
    strain0<<0.002,0.001,-0.0005,0.0007,0.0003,-0.0004;
    strainRate=0.1*strain0;
    stress0<<0.001,-0.0004,0.0007,0.0003,-0.0002,0.0004;
    stressRate=0.2*stress0;
    if(mode==0)
    {
        alpha.setZero();
    }
    else if(mode==1)
    {
        alpha.setConstant(1.e20);
    }
    else
    {
        alpha<<0,1.e20,2,0.5,1.e20,0;
    }
    if(irregular)
    {
        strain0.setZero();
        strainRate.setZero();
        strainRate(2)=0.001;
        stress0.setZero();
        stressRate.setZero();
        alpha.setZero();
        alpha(2)=1.e20;
    }
    const auto writeLoading=[&](bool fem,const V& ratios,const V& stresses,const V& rates)
    {
        std::ofstream out(input/"ElasticDeformation.txt");
        out<<std::setprecision(17)<<"useElasticDeformationFEM="<<fem<<";\ninertiaReliefPenaltyFactor=0.001;\nExternalStrain0="<<strain0.transpose()<<";\nExternalStrainRate="<<strainRate.transpose()<<";\nExternalStress0="<<stresses.transpose()<<";\nExternalStressRate="<<rates.transpose()<<";\nstiffnessRatio="<<ratios.transpose()<<";\n";
    };
    {
        std::ofstream dc(input/"DefectiveCrystal.txt");
        dc<<"physics=ElasticDeformation;\nuseFEM=1;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=";
        for(int i=0;i<rank;++i)
        {
            dc<<"0 ";
        }
        dc<<";\nEwaldLengthFactor=0;\n";
        const auto matrixText=[](const Eigen::Matrix3d& matrix)
        {
            std::ostringstream result;
            result<<std::setprecision(17);
            for(int row=0;row<3;++row)
            {
                for(int col=0;col<3;++col)
                {
                    result<<matrix(row,col)<<" ";
                }
            }
            return result.str();
        };
        std::ofstream poly(input/"polycrystal.txt");
        poly<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1="<<matrixText(Q)<<";\nC2G2="<<matrixText(Q)<<";\nF="<<matrixText(3671.9477120460206*Q)<<";\nX0=0 0 0;\nperiodicFaceIDs=";
        for(int seed : seeds)
        {
            poly<<seed<<" ";
        }
        poly<<";\n";
        if(irregular)
        {
            for(int grain=3;grain<=10;++grain)
            {
                poly<<"C2G"<<grain<<"="<<matrixText(Q)<<";\n";
            }
        }
        writeLoading(true,alpha,stress0,stressRate);
    }
    {
        DefectiveCrystalBase<3> base(folder.string());
        require(base.fe!=nullptr,"Explicit FEM request was disabled in a periodic domain");
        require(base.periodicLatticeBasis.cols()==rank,"Unexpected active periodic rank");
        // Each physics configuration has its own immutable context. Do not
        // change the meaning of an existing base by rewriting its input file.
        std::unique_ptr<DefectiveCrystalBase<3>> sourceBase;
        std::unique_ptr<MicrostructureContainer<3>> sourceContainer;
        MicrostructureContainer<3> container(base);
        ElasticDeformation<3> physics(container);
        physics.initializeSolver();
        auto u=physics.elasticDeformationFEM->u;
        if(withPin)
        {
            const size_t node=base.fe->periodicNodePairs().front().sink;
            for(int c=0;c<3;++c)
            {
                u->dirichletConditions()[3*node+c]=0.0;
            }
        }
        physics.zDot=Eigen::VectorXd::Zero(u->gSize());
        const M C=physics.elasticDeformationFEM->C;
        V eigenstrain=V::Zero();
        V sourceStress=V::Zero();
        if(withSource)
        {
            sourceStress<<0.0002,0.0003,-0.0001,0.0002,0,-0.0001;
            writeLoading(false,V::Zero(),sourceStress,V::Zero());
            StaticID<Polyhedron<3,3-1>>::force_count(0);
            sourceBase=std::make_unique<DefectiveCrystalBase<3>>(folder.string());
            sourceContainer=std::make_unique<MicrostructureContainer<3>>(*sourceBase);
            auto source=std::make_shared<UniformSource>(*sourceContainer);
            eigenstrain<<0.0004,0.0001,0,0.0002,0,0;
            source->plastic=base.voigtTraits.v2m(eigenstrain,true);
            source->uniformLoadController->gs=eigenstrain;
            container.microstructures().push_back(source);
            writeLoading(true,alpha,stress0,stressRate);
        }
        size_t count=0;
        for(double time : {7.0,12.0,12.0})
        {
            base.simulationParameters.totalTime=time;
            if(sourceBase)
            {
                sourceBase->simulationParameters.totalTime=time;
            }
            physics.solve();
            for(const auto& bc : u->dirichletConditions())
            {
                require(std::abs(u->dofVector()(bc.first)-bc.second)<1.e-9,"Mixed elastic controller moved a pinned node");
            }
            M system=C;
            V target=stress0+time*stressRate+C*eigenstrain;
            for(int i=0;i<6;++i)
            {
                if(alpha(i)>=1.e15)
                {
                    system.row(i).setZero();
                    system(i,i)=1.0;
                    target(i)=strain0(i)+time*strainRate(i);
                }
                else
                {
                    const double device=alpha(i)*C(i,i);
                    system(i,i)+=device;
                    target(i)+=device*(strain0(i)+time*strainRate(i));
                }
            }
            const V expectedTotal=system.fullPivLu().solve(target);
            const V offset=eigenstrain+C.ldlt().solve(sourceStress);
            const V expected=expectedTotal-offset;
            require((physics.femLoadController->averageStrain()-expected).norm()<1.e-9,"Mixed FEM average strain differs from uniform controller");
            require((base.voigtTraits.m2v(physics.averageStress(),false)+sourceStress-C*(expectedTotal-eigenstrain)).norm()<1.e-9,"Mixed FEM average stress differs from uniform controller");
            const Eigen::Matrix3d expectedStrain=base.voigtTraits.v2m(expected,true);
            const double length=(base.mesh.xMax()-base.mesh.xMin()).norm();
            Eigen::MatrixXd coordinates(base.fe->nodeSize(),4),values(base.fe->nodeSize(),3);
            for(const auto& node : base.fe->nodes())
            {
                coordinates.block<1,3>(node.gID,0)=node.P0.transpose()/length;
                coordinates(node.gID,3)=1.0;
                values.row(node.gID)=u->dofVector().segment<3>(3*node.gID).transpose();
            }
            const Eigen::MatrixXd fit=coordinates.colPivHouseholderQr().solve(values);
            const Eigen::Matrix3d expectedGradient=fit.topRows(3).transpose()/length;
            require((0.5*(expectedGradient+expectedGradient.transpose())-expectedStrain).norm()<1.e-9,"Mixed loading produced incorrect homogeneous strain");
            if(!withPin)
            {
                require((expectedGradient-expectedGradient.transpose()).norm()<1.e-9,"Global stabilization introduced mean rotation");
                Eigen::Vector3d centroid=Eigen::Vector3d::Zero();
                double volume=0.0;
                for(const auto& entry : base.fe->elements())
                {
                    centroid+=entry.second.simplex.vol0*entry.second.simplex.center();
                    volume+=entry.second.simplex.vol0;
                }
                centroid/=volume;
                const double meanDisplacement=(expectedGradient*centroid+fit.row(3).transpose()).norm();
                const double displacementScale=std::max(1.0,values.norm()/std::sqrt(double(values.rows())));
                // The volume-normalized penalty has weight 1e-3 against dimensional
                // stiffness; test drift relative to the physical displacement scale.
                require(meanDisplacement<1.e-6*displacementScale,"Mean affine displacement is not zero");
                if(irregular)
                {
                    require(std::abs(expectedGradient(2,0))+std::abs(expectedGradient(2,1))<1.e-9,"Horizontal plane has a spurious u3 gradient");
                }
            }
            require((coordinates*fit-values).norm()<2.e-7*std::max(1.0,values.norm()),"Mixed loading introduced nonuniform deformation");
            Eigen::Matrix3d projection=Eigen::Matrix3d::Zero();
            if(rank>0)
            {
                const auto& B=base.periodicLatticeBasis;
                projection=B*(B.transpose()*B).inverse()*B.transpose();
            }
            require((u->periodicGradient()-expectedGradient*projection).norm()<1.e-8,"Solved affine gradient used inactive periodic directions");
            for(const auto& pair : base.fe->periodicNodePairs())
            {
                const Eigen::Vector3d jump=u->dofVector().segment<3>(3*pair.sink)-u->dofVector().segment<3>(3*pair.source);
                require((jump-u->periodicGradient()*pair.shift).norm()<1.e-8,"Solved periodic jump is inconsistent");
            }
            if(count>0)
            {
                require(physics.weakProblemSolver->factorizationCount()==count,"Changing mixed load refactorized stiffness");
            }
            count=physics.weakProblemSolver->factorizationCount();
        }
        if(mode==2)
        {
            // Independently check the finite-ratio case using the actual uniform controller.
            V finiteAlpha=alpha;
            finiteAlpha(1)=3.0;
            finiteAlpha(4)=4.0;
            writeLoading(true,finiteAlpha,stress0,stressRate);
            ElasticDeformation<3> finite(container);
            finite.zDot=Eigen::VectorXd::Zero(u->gSize());
            finite.solve();
            UniformController<6> reference(0.0,C,finiteAlpha,strain0,strainRate,stress0,stressRate);
            reference.gs=eigenstrain;
            const V offset=eigenstrain+C.ldlt().solve(sourceStress);
            require((finite.femLoadController->averageStrain()+offset-reference.grad(base.simulationParameters.totalTime)).norm()<1.e-9,"Finite FEM device stiffness disagrees with UniformController");
        }
    }
    std::filesystem::remove_all(folder);
}
int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository path");
        for(int rank=0;rank<=3;++rank)
        {
            exercise(argv[1],rank);
        }
        for(int rank=0;rank<=3;++rank)
        {
            for(int mode=0;mode<3;++mode)
            {
                physicsIntegration(argv[1],rank,false,mode);
            }
        }
        physicsIntegration(argv[1],1,true,2);
        physicsIntegration(argv[1],2,true,2,true);
        physicsIntegration(argv[1],1,true,2,false,true);
        for(int rank=0;rank<=2;++rank)
        {
            physicsIntegration(argv[1],rank,false,1,false,false,true);
        }
        std::cout<<"Global elastic penalty tests passed\n";
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
