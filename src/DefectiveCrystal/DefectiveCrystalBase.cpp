/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * Additional contributors:
 *  Nicholas Huebner Julian <njulian@lanl.gov>
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalBase_CPP_
#define model_DefectiveCrystalBase_CPP_

#include <ElasticDeformation/ElasticDeformationBase.h>
#include <ClusterDynamics/ClusterDynamicsBase.h>
#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <CTM.h>
#include <StrUtilities.h>
#include <sstream>
#include <algorithm>
namespace model
{

    namespace
    {
        std::vector<std::string> parsePhysics(const std::string& input)
        {
            std::vector<std::string> result;
            std::istringstream stream(input);
            std::string name;
            while(stream >> name)
            {
                name=StrUtilities::lowercase(name);
                if(name=="none")
                {
                    continue;
                }
                if(name!="dislocationdynamics" && name!="clusterdynamics" && name!="elasticdeformation" && name!="inclusionmicrostructure")
                {
                    throw std::runtime_error("Unknown physics '"+name+"'.");
                }
                if(std::find(result.begin(),result.end(),name)==result.end())
                {
                    result.push_back(name);
                }
            }
            return result;
        }
    }

    template <int dim>
    std::vector<Eigen::Matrix<double,dim,1>> getPeriodicShifts(const Eigen::Matrix<double,dim,Eigen::Dynamic>& shiftVectors,const std::vector<int>& periodicImageSize)
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        std::vector<VectorDim> temp;
        temp.push_back(VectorDim::Zero());
        if(periodicImageSize.size())
        {
            std::cout<<"Box periodicity basis ("<<shiftVectors.size()<<") in columns:"<<std::endl<<shiftVectors<<std::endl;
            if(shiftVectors.cols()!=int(periodicImageSize.size()))
            {
                std::cout<<"shiftVectors.size()="<<shiftVectors.size()<<std::endl;
                std::cout<<"periodicImageSize.size()="<<periodicImageSize.size()<<std::endl;
                throw std::runtime_error("shiftVectors.size() must equal periodicImageSize.size()");
            }
            
            for(int k=0;k<shiftVectors.cols();++k)
            {
                const auto shiftVector(shiftVectors.col(k));
                const int imageSize(std::abs(periodicImageSize[k]));
                std::vector<VectorDim> newTemp;
                for(const auto& v:temp)
                {// grab existing shift vectors
                    for(int i=-imageSize;i<=imageSize;++i)
                    {// add current shift times corresponding image size
                        newTemp.push_back(v+i*shiftVector);
                    }
                }
                temp.swap(newTemp);
            }
        }
        return temp;
    }

    template <int _dim>
    DefectiveCrystalBase<_dim>::DefectiveCrystalBase( const std::string& folderName) :
    /* init */ simulationParameters( folderName)
    /* init */,voigtTraits((typename SymmetricVoigtTraits<3>::VoigtSizeMatrixType()<<0,0,1,1,2,2,0,1,1,2,0,2).finished())
    /* init */,meshOwner(std::make_shared<const SimplicialMesh<dim>>(simulationParameters.traitsIO.meshFile,
                    TextFileParser(simulationParameters.traitsIO.polyFile).readMatrix<double>("F",_dim,_dim,true),
                    TextFileParser(simulationParameters.traitsIO.polyFile).readMatrix<double>("X0",1,_dim,true).transpose(),
                    simulationParameters.periodicFaceIDs))
    /* init */,mesh(*meshOwner)
    /* init */,periodicImageSize(mesh.periodicBasis().cols()>0? TextFileParser(simulationParameters.traitsIO.dcFile).readArray<int>("periodicImageSize",true) : std::vector<int>())
    /* init */,periodicLatticeBasis(mesh.periodicBasis())
    /* init */,periodicLatticeReciprocalBasis(periodicLatticeBasis*(periodicLatticeBasis.transpose()*periodicLatticeBasis).inverse())
    /* init */,periodicShifts(getPeriodicShifts(periodicLatticeBasis,periodicImageSize))
    /* init */,polyOwner(std::make_shared<const Polycrystal<dim>>(simulationParameters.traitsIO.polyFile,mesh))
    /* init */,poly(*polyOwner)
    /* init */,fe(simulationParameters.useFEM ? new FiniteElement<ElementType>(meshOwner) : nullptr)
    /* init */,activePhysics(parsePhysics(simulationParameters.physics))
    /* init */,fmmContext(std::make_shared<FMMContext<dim>>(mesh.xMin(),mesh.xMax(),simulationParameters.traitsIO.dcFile))
    /* init */,dislocationDynamicsBase(hasPhysics("DislocationDynamics") ? std::make_shared<DislocationDynamicsBase<dim>>(meshOwner,polyOwner,periodicLatticeBasis,simulationParameters.traitsIO.dcFile) : nullptr)
    /* init */,clusterDynamicsBase(hasPhysics("ClusterDynamics") ? std::make_shared<ClusterDynamicsBase<dim>>(*this) : nullptr)
    /* init */,elasticDeformationBase(hasPhysics("ElasticDeformation") ? std::make_shared<ElasticDeformationBase<dim>>(*this) : nullptr)
    {
        if(!mesh.simplices().size())
        {
            throw std::runtime_error("Mesh is empty");
        }
        
        for(const auto& region : mesh.regions())
        {
            for(const auto& face : region.second->faces())
            {
                if(face.second->isExternal() && face.second->periodicFacePair.second)
                {
                    const Eigen::Matrix<double,Eigen::Dynamic,1> pd(periodicLatticeReciprocalBasis.transpose()*face.second->periodicFacePair.first);
                    const Eigen::Matrix<double,Eigen::Dynamic,1> pn(pd.array().round());
                    if((pd-pn).norm()>FLT_EPSILON*pd.norm())
                    {
                        std::cout<<"region="<<region.second->regionID<<std::endl;
                        std::cout<<"face="<<face.second->sID<<std::endl;
                        std::cout<<"pd="<<pd.transpose()<<std::endl;
                        std::cout<<"pn="<<pn.transpose()<<std::endl;
                        throw std::runtime_error("Face shift must be an intiger combination of the periodicity vectors");
                    }
                }
            }
        }
        
        std::cout<<"periodicLatticeBasis="<<periodicLatticeBasis<<std::endl;
        std::cout<<"periodicLatticeReciprocalBasis="<<periodicLatticeReciprocalBasis<<std::endl;
        std::cout<<"periodicShifts("<<periodicShifts.size()<<")="<<std::endl;
        for(const auto& shift : periodicShifts)
        {
            std::cout<<shift.transpose()<<std::endl;
        }
    }

    template <int _dim>
    Eigen::VectorXd DefectiveCrystalBase<_dim>::periodicCoordinates(const VectorDim& x) const
    {
        const VectorDim dx(x-mesh.xCenter());
        return periodicLatticeReciprocalBasis.transpose()*dx;
    }

    template<int dim>
    bool DefectiveCrystalBase<dim>::hasPhysics(const std::string& name) const
    {
        return std::find(activePhysics.begin(),activePhysics.end(),StrUtilities::lowercase(name))!=activePhysics.end();
    }

    template<int dim>
    auto DefectiveCrystalBase<dim>::dislocationDynamics() const -> DislocationDynamicsBase<dim>&
    {
        if(!dislocationDynamicsBase)
        {
            throw std::runtime_error("DislocationDynamics physics is not enabled.");
        }
        return *dislocationDynamicsBase;
    }

    template<int dim>
    auto DefectiveCrystalBase<dim>::clusterDynamics() const -> ClusterDynamicsBase<dim>&
    {
        if(!clusterDynamicsBase)
        {
            throw std::runtime_error("ClusterDynamics physics is not enabled.");
        }
        return *clusterDynamicsBase;
    }

    template<int dim>
    auto DefectiveCrystalBase<dim>::elasticDeformation() const -> ElasticDeformationBase<dim>&
    {
        if(!elasticDeformationBase)
        {
            throw std::runtime_error("ElasticDeformation physics is not enabled.");
        }
        return *elasticDeformationBase;
    }

    template struct DefectiveCrystalBase<3>;

} // namespace model
#endif
