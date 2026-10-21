#include <SimplicialMeshActor.h>
#include <QApplication>
#include <filesystem>
#include <iostream>
#include <set>
void require(bool condition,const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
int main(int argc,char** argv)
{
    QApplication app(argc,argv);
    try
    {
        require(argc==2,"Expected repository root");
        const std::filesystem::path root(argv[1]);
        for(const auto& file : {"bicrystal_48.msh","unitCube24.msh"})
        {
            for(bool periodic : {false,true})
            {
                model::SimplicialMesh<3> mesh((root/"Library/Meshes"/file).string(),
                    Eigen::Matrix3d::Identity()*100.0,Eigen::Vector3d::Zero(),
                    periodic?std::set<int>{-1}:std::set<int>{});
                auto window=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
                window->SetReadyForRendering(false);
                auto renderer=vtkSmartPointer<vtkRenderer>::New();
                window->AddRenderer(renderer);
                model::SimplicialMeshActor actor(window,renderer,mesh);
                std::set<const model::Simplex<3,2>*> internal,external,sameGrain,crossGrain;
                for(const auto& region : mesh.regions())
                {
                    for(const auto& entry : region.second->faces())
                    {
                        const auto& face=*entry.second;
                        auto& target=!face.isExternal()?internal:
                            !face.periodicFacePair.second?external:
                            face.regionIDs==face.periodicFacePair.second->regionIDs?sameGrain:crossGrain;
                        target.insert(face.simplices().begin(),face.simplices().end());
                    }
                }
                const auto expected=2*internal.size()+external.size()+crossGrain.size();
                require(actor.gbTrianglePolyData->GetNumberOfPolys()==static_cast<vtkIdType>(expected),
                        "Grain surface actor includes same-grain seams or omits physical faces");
                require(actor.gbColors->GetNumberOfTuples()==static_cast<vtkIdType>(expected),"Face colors mismatch");
                if(periodic)
                {
                    require(!sameGrain.empty(),"Missing same-grain periodic test faces");
                    if(std::string(file)=="bicrystal_48.msh")
                    {
                        require(!internal.empty() && !crossGrain.empty(),"Missing bicrystal GB types");
                    }
                    else
                    {
                        require(expected==0,"Fully periodic single grain still has shaded boundaries");
                    }
                }
                else
                {
                    require(!external.empty(),"Missing free surface test faces");
                }
            }
        }
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
}
