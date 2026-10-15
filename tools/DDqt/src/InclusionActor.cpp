/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionActor_cpp_
#define model_InclusionActor_cpp_

#include <Eigen/Dense>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkPolyhedron.h>
#include <vtkIdList.h>
//#include <vtkXMLUnstructuredGridWriter.h>
#include <InclusionActor.h>
#include <EllipsoidalInclusion.h>
#include <EllipsoidalInclusionPolyData.h>
#include <vtkAppendPolyData.h>


// VTK documentation
// http://vtk.1045678.n5.nabble.com/VTK-slow-to-display-300-vtkSphereSource-in-real-time-td5740730.html
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/Scaleglyph3D
// http://hiraku0n.blogspot.com/2012/10/writing-scalar-and-vector-field-in-vtk.html
// http://vtk.1045678.n5.nabble.com/Glyphing-vtkImageData-scalars-3D-as-arrows-td3199837.html

#include <DDqtControls.h>

namespace model
{

    void InclusionActor::updateConfiguration()
    {
        std::cout<<"Updating inclusions..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();

        // Replace the previous frame, including when no ellipsoids remain.
        vtkNew<vtkPolyData> emptyEllipsoids;
        ellipsoidalMapper->SetInputData(emptyEllipsoids);
        vtkNew<vtkAppendPolyData> ellipsoids;

        if(inclusions)
        {
            
                    
            vtkSmartPointer<vtkPoints> points(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkFloatArray> diameters(vtkSmartPointer<vtkFloatArray>::New());
            vtkSmartPointer<vtkFloatArray> colors(vtkSmartPointer<vtkFloatArray>::New());
            diameters->SetName("diameters");
            colors->SetName("colors");
            
            vtkNew<vtkUnstructuredGrid> ugrid;
            vtkNew<vtkPoints> polyhedronPoints;
//            vtkIdType polyhedronPointsIDs[inclusions->polyhedronInclusionNodes().size()];
            std::vector<vtkIdType> polyhedronPointsIDs(inclusions->polyhedronInclusionNodes().size());
            long long k=0;
            for(const auto& node : inclusions->polyhedronInclusionNodes())
            {
                polyhedronPoints->InsertNextPoint(node.second.P.data());
                polyhedronPointsIDs[k]=k;
                k++;
            }
            ugrid->SetPoints(polyhedronPoints);

                        
            for(const auto& inclusion : inclusions->eshelbyInclusions())
            {
                
                auto* sphericalDerived = dynamic_cast<SphericalInclusion<dim>*>(inclusion.second.get());
                if (sphericalDerived)
                {
                    points->InsertNextPoint(sphericalDerived->C(0), sphericalDerived->C(1), sphericalDerived->C(2));
                    diameters->InsertNextValue(2.0*sphericalDerived->a);  // origin of arrow
                    colors->InsertNextValue(sphericalDerived->phaseID);
                }
                
                auto* ellipsoidalDerived = dynamic_cast<EllipsoidalInclusion<dim>*>(inclusion.second.get());
                if(ellipsoidalDerived)
                {
                    ellipsoids->AddInputData(ellipsoidalInclusionPolyData(
                        ellipsoidalDerived->C,ellipsoidalDerived->axes));
                }

                auto* polyhedronDerived = dynamic_cast<PolyhedronInclusion<dim>*>(inclusion.second.get());
                if (polyhedronDerived)
                {
//                    configIO.polyhedronInclusions().emplace_back(*polyhedronDerived);
                    
                    vtkNew<vtkIdList> faces;
//                    for(const auto& pair2 : pair1.second)
//                    {
//                        faces->InsertNextId(pair2.second.size());
//                        for(const auto& nID : pair2.second)
//                        {
//                            faces->InsertNextId(nID);
//                        }
//                    }
//                    ugrid->InsertNextCell(VTK_POLYHEDRON, inclusions->polyhedronInclusionNodes().size(), polyhedronPointsIDs, pair1.second.size(), faces->GetPointer(0));
                    
                    for(const auto& facePair : polyhedronDerived->faces)
                    {
                        faces->InsertNextId(facePair.second.size());

                        for(const auto& vertexPair : facePair.second)
                        {
                            faces->InsertNextId(vertexPair.first);
                        }
                    }
                    ugrid->InsertNextCell(VTK_POLYHEDRON, inclusions->polyhedronInclusionNodes().size(), polyhedronPointsIDs.data(), polyhedronDerived->faces.size(), faces->GetPointer(0));

                }
            }
            
//            for(const auto& inclusion : configFields.configIO.sphericalInclusions())
//            {
//                points->InsertNextPoint(inclusion.C(0), inclusion.C(1), inclusion.C(2));
//                diameters->InsertNextValue(2.0*inclusion.a);  // origin of arrow
//                colors->InsertNextValue(inclusion.phaseID);
//            }
            
            grid->SetPoints(points);
            grid->GetPointData()->AddArray(diameters);
            grid->GetPointData()->SetActiveScalars("diameters"); // to set radius first
            grid->GetPointData()->AddArray(colors);
            
//            vtkNew<vtkUnstructuredGrid> ugrid;
//            vtkNew<vtkPoints> polyhedronPoints;
//            vtkIdType polyhedronPointsIDs[inclusions->polyhedronInclusionNodes().size()];
//            long long k=0;
//            for(const auto& node : inclusions->polyhedronInclusionNodes())
//            {
//                polyhedronPoints->InsertNextPoint(node.second.P.data());
//                polyhedronPointsIDs[k]=k;
//                k++;
//            }
//            ugrid->SetPoints(polyhedronPoints);
            
//            std::map<size_t,std::map<size_t,std::vector<size_t>>> faces;
//            for(const auto& edge : inclusions->polyhedronInclusionEdges())
//            {
//                const size_t& iID(edge.inclusionID);
//                const size_t& fID(edge.faceID);
//                const size_t& sourceID(edge.sourceID);
//                faces[iID][fID].push_back(sourceID);
//            }
//            
//            for(const auto& pair1 : faces)
//            {
//                vtkNew<vtkIdList> faces;
//                for(const auto& pair2 : pair1.second)
//                {
//                    faces->InsertNextId(pair2.second.size());
//                    for(const auto& nID : pair2.second)
//                    {
//                        faces->InsertNextId(nID);
//                    }
//                }
//                ugrid->InsertNextCell(VTK_POLYHEDRON, inclusions->polyhedronInclusionNodes().size(), polyhedronPointsIDs, pair1.second.size(), faces->GetPointer(0));
//            }
            
//            polyhedronMapper->SetInputData(polyhedronPoints->GetOutputPort());
            polyhedronMapper->SetInputData(ugrid);
            if(ellipsoids->GetNumberOfInputConnections(0)>0)
            {
                ellipsoids->Update();
                ellipsoidalMapper->SetInputData(ellipsoids->GetOutput());
            }
            
        }
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }

    InclusionActor::InclusionActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer,const DefectiveCrystal<3>& defectiveCrystal_in) :
    /* init */ renderWindow(renWin)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,showInclusions(new QCheckBox(this))
    /* init */,sliderInclusionOpacity(new QSlider(this))
    /* init */,grid(vtkSmartPointer<vtkUnstructuredGrid>::New())
    /* init */,sphereSource(vtkSmartPointer<vtkSphereSource>::New())
    /* init */,glyph3D(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,actor(vtkSmartPointer<vtkActor>::New())
    /* init */,lookUpColors(vtkSmartPointer<vtkLookupTable>::New())
    /* init */,ellipsoidalMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,ellipsoidalActor(vtkSmartPointer<vtkActor>::New())
    /* init */,polyhedronMapper(vtkSmartPointer<vtkDataSetMapper>::New())
    /* init */,polyhedronActor(vtkSmartPointer<vtkActor>::New())
    /* init */,defectiveCrystal(defectiveCrystal_in)
    /* init */,inclusions(defectiveCrystal.template getUniqueTypedMicrostructure<InclusionMicrostructure<3>>())
    {
        showInclusions->setChecked(true);
        showInclusions->setText("show inclusions");
        auto *inclusionBox =
            controls::box("Inclusions", this, "Show spherical, ellipsoidal and polyhedral inclusions.", showInclusions);
        mainLayout->addWidget(inclusionBox, 0, 0);
        mainLayout->setRowStretch(1, 1);
        this->setLayout(mainLayout);
        connect(showInclusions,SIGNAL(stateChanged(int)), this, SLOT(modify()));
        
        sliderInclusionOpacity->setMinimum(0);
        sliderInclusionOpacity->setMaximum(10);
        sliderInclusionOpacity->setValue(5);
        sliderInclusionOpacity->setOrientation(Qt::Horizontal);
        connect(sliderInclusionOpacity,SIGNAL(valueChanged(int)), this, SLOT(modify()));
        controls::field(controls::grid(inclusionBox), 0, "Opacity", sliderInclusionOpacity,
                        "Set inclusion opacity from transparent (0) to opaque (10).");

        sphereSource->SetPhiResolution(20);
        sphereSource->SetThetaResolution(20);
        glyph3D->SetInputData(grid);
        glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
        mapper->SetInputConnection(glyph3D->GetOutputPort());
        mapper->ScalarVisibilityOff();
        actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
        actor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white
        actor->SetMapper(mapper);
        
        ellipsoidalMapper->ScalarVisibilityOff();
        vtkNew<vtkPolyData> emptyEllipsoids;
        ellipsoidalMapper->SetInputData(emptyEllipsoids);
        ellipsoidalActor->SetMapper(ellipsoidalMapper);
        ellipsoidalActor->GetProperty()->SetColor(0.5,0.5,0.5);
        ellipsoidalActor->GetProperty()->SetOpacity(0.4);

        polyhedronActor->SetMapper(polyhedronMapper);
        polyhedronActor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
        polyhedronActor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white
                
        renderer->AddActor(actor);
        renderer->AddActor(ellipsoidalActor);
        renderer->AddActor(polyhedronActor);
    }

    void InclusionActor::modify()
    {
        actor->SetVisibility(showInclusions->isChecked());
        actor->GetProperty()->SetOpacity(sliderInclusionOpacity->value()/10.0);

        ellipsoidalActor->SetVisibility(showInclusions->isChecked());
        ellipsoidalActor->GetProperty()->SetOpacity(sliderInclusionOpacity->value()/10.0);

        polyhedronActor->SetVisibility(showInclusions->isChecked());
        polyhedronActor->GetProperty()->SetOpacity(sliderInclusionOpacity->value()/10.0);

        
        renderWindow->Render();
    }

}
#endif
