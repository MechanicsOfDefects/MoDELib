/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLoopActor_cpp_
#define model_NetworkLoopActor_cpp_

#include <DislocationNetwork.h>
#include <iostream>
#include <map>
#include <deque>
#include <string>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLabeledDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkLine.h>
#include <vtkCellData.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <QLabel>

#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <MeshPlane.h>
#include <NetworkLoopActor.h>
#include <DislocationDynamics/DislocationLoopPatches.h>
#include <LoopAreaTriangulation.h>

#include <DDqtControls.h>

namespace model
{


/**********************************************************************/
NetworkLoopActor::NetworkLoopActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer,
                                   const DefectiveCrystal<3>& defectiveCrystal_in) :
/* init */ renderWindow(renWin)
/* init */,mainLayout(new QGridLayout(this))
/* init */,showLoops(new QCheckBox(this))
/* init */,colorScheme(new QComboBox(this))
/* init */,manualColorButton(new QPushButton("Choose color…", this))
/* init */,slippedAreaBox(new controls::GroupBox(tr("&Slipped Area")))
/* init */,sliderSlippedArea(new QSlider(this))
/* init */,meshAreaBox(new controls::GroupBox(tr("&Slipped Mesh")))
/* init */,nodePolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,nodeGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
/* init */,nodeMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,nodeActor(vtkSmartPointer<vtkActor>::New())
/* init */,labelPolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,labelMapper(vtkSmartPointer<vtkLabeledDataMapper>::New())
/* init */,labelActor(vtkSmartPointer<vtkActor2D>::New())
/* init */,loopPolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,loopMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,loopActor(vtkSmartPointer<vtkActor>::New())
/* init */,areaPolyData(vtkSmartPointer<vtkPolyData>::New())
/* init */,areaMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,areaActor(vtkSmartPointer<vtkActor>::New())
/* init */,meshPolydata(vtkSmartPointer<vtkPolyData>::New())
/* init */,meshMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
/* init */,meshActor(vtkSmartPointer<vtkActor>::New())
/* init */,globalMeshSizeEdit(new QLineEdit("100000000"))
/* init */,localMeshSizeEdit(new QLineEdit("1.0"))
/* init */,defectiveCrystal(defectiveCrystal_in)
/* init */,dislocationNetwork(defectiveCrystal.template getUniqueTypedMicrostructure<DislocationNetwork<3>>())
{
    
    nodeGlyphs->SetInputData(nodePolyData);
    nodeGlyphs->SetSourceConnection(vtkSmartPointer<vtkSphereSource>::New()->GetOutputPort());
    nodeGlyphs->ScalingOn();
    nodeGlyphs->SetScaleModeToScaleByVector();
    nodeGlyphs->SetScaleFactor(10);
    nodeGlyphs->SetColorModeToColorByScalar();
    nodeGlyphs->Update();
    nodeMapper->SetInputConnection(nodeGlyphs->GetOutputPort());
    nodeActor->SetMapper(nodeMapper);
    nodeActor->SetVisibility(false);

    // Labels
    labelMapper->SetInputData(labelPolyData);
    labelMapper->SetLabelModeToLabelScalars();
    labelMapper->SetLabelFormat("%1.0f");
    labelMapper->GetLabelTextProperty()->SetFontSize(20);
    labelActor->SetMapper(labelMapper);
    labelActor->GetProperty()->SetColor(0.0, 0.0, 0.0); //(R,G,B)
    labelActor->SetVisibility(false);

    
    showLoops->setChecked(false);
    showLoops->setText("Loops");
    
    slippedAreaBox->setCheckable(true);
    slippedAreaBox->setChecked(false);
    
    meshAreaBox->setCheckable(true);
    meshAreaBox->setChecked(false);

    
    //            showSlippedArea->setChecked(false);
    //            showSlippedArea->setText("SlippedArea");
    //            sliderSlippedArea->setEnabled(false);
    sliderSlippedArea->setMinimum(0);
    sliderSlippedArea->setMaximum(10);
    sliderSlippedArea->setValue(5);
    sliderSlippedArea->setOrientation(Qt::Horizontal);
    
    //            colorSlippedArea->setWindowFlags(Qt::Widget );

    QGridLayout *slippedAreaLayout = new QGridLayout();
    controls::field(slippedAreaLayout, 0, "Opacity", sliderSlippedArea,
                    "Opacity of slipped-area polygons: 0 is transparent and 10 is opaque.");

    //            slippedAreaLayout->addWidget(colorSlippedArea);
    slippedAreaBox->setLayout(slippedAreaLayout);
    
    /* a few options that we must set for it to work nicely */
    //            colorSlippedArea->setOptions(
    //                            /* do not use native dialog */
    //                            QColorDialog::DontUseNativeDialog
    //                            /* you don't need to set it, but if you don't set this
    //                                the "OK" and "Cancel" buttons will show up, I don't
    //                                think you'd want that. */
    ////                            | QColorDialog::NoButtons
    //                );
    
    //            QGroupBox *groupBox = new QGroupBox(tr("&Push Buttons"));
    //                groupBox->setCheckable(true);
    //                groupBox->setChecked(true);
    
    QGridLayout *meshAreaLayout = new QGridLayout();
    controls::field(
        meshAreaLayout, 0, "Global size", globalMeshSizeEdit,
        "Maximum triangle area away from the loop boundary, in b squared. Press Enter to rebuild the slipped mesh.");
    controls::field(meshAreaLayout, 1, "Local size", localMeshSizeEdit,
                    "Spacing between refinement points along dislocation segments, in Burgers-vector units. Press "
                    "Enter to rebuild the slipped mesh.");

    meshAreaBox->setLayout(meshAreaLayout);

    auto *loopsBox =
        controls::box("Loops", this, "Show or hide loop outlines, loop nodes and loop-node IDs.", showLoops);
    auto *appearanceBox = controls::box(
        "Appearance", this, "Color settings shared by outlines, node markers, slipped areas and slipped meshes.");
    mainLayout->addWidget(loopsBox, 0, 0, 1, 2);
    colorScheme->setObjectName("loopColorScheme");
    colorScheme->addItems({"Manual", "Slip system"});
    manualColorButton->setObjectName("loopManualColor");
    controls::field(controls::grid(appearanceBox), 0, "Coloring", colorScheme,
                    "Manual uses one chosen color. Slip system follows the SlipSystems palette and its live updates.");
    controls::grid(appearanceBox)->addWidget(manualColorButton, 1, 0, 1, 2);
    manualColorButton->setToolTip("Choose the manual color for all loop geometry. Default is light green.");
    mainLayout->addWidget(appearanceBox, 1, 0, 1, 2);
    mainLayout->addWidget(slippedAreaBox, 2, 0, 1, 2);
    mainLayout->addWidget(meshAreaBox, 3, 0, 1, 2);
    mainLayout->setRowStretch(4, 1);
    slippedAreaBox->setToolTip("Show or hide the area swept by each dislocation loop.");
    meshAreaBox->setToolTip("Show or hide triangulated slipped areas; enabling this group builds the mesh.");
    connect(colorScheme, &QComboBox::currentIndexChanged, this,
            [this](int)
            {
                applyColors();
                modify();
            });
    connect(manualColorButton, &QPushButton::clicked, this,
            [this]()
            {
                const QColor chosen = QColorDialog::getColor(manualColor, this, "Loop color");
                if(chosen.isValid())
                {
                    setManualColor(chosen);
                }
            });

    //            mainLayout->addWidget(sliderSlippedArea,1,1,1,1);
    //            mainLayout->addWidget(colorSlippedArea,1,2,1,1);
    
    
    this->setLayout(mainLayout);
    
    connect(showLoops,SIGNAL(stateChanged(int)), this, SLOT(modify()));
    //            connect(slippedAreaBox,SIGNAL(stateChanged(int)), this, SLOT(modify()));
    connect(slippedAreaBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
    const auto refreshMesh = [this]()
    {
        if(meshAreaBox->isChecked())
        {
            updateConfiguration(slipSystemColors);
        }
        modify();
    };
    connect(meshAreaBox, &QGroupBox::toggled, this, refreshMesh);
    connect(globalMeshSizeEdit, &QLineEdit::returnPressed, this, refreshMesh);
    connect(localMeshSizeEdit, &QLineEdit::returnPressed, this, refreshMesh);

    connect(sliderSlippedArea,SIGNAL(valueChanged(int)), this, SLOT(modify()));
//    connect(globalMeshSizeEdit,SIGNAL(returnPressed()), this, SLOT(updateConfiguration()));
//    connect(localMeshSizeEdit,SIGNAL(returnPressed()), this, SLOT(updateConfiguration()));

    
    loopMapper->SetInputData(loopPolyData);
    loopActor->SetMapper(loopMapper);
    loopActor->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
    loopActor->SetVisibility(false);
    
    areaMapper->SetInputData(areaPolyData);
    areaActor->SetMapper(areaMapper);
//    areaActor->GetProperty()->SetColor(0.0, 1.0, 1.0); //(R,G,B)
    areaActor->SetVisibility(false);
    
    meshPolydata->Allocate();
    meshMapper->SetInputData(meshPolydata);
    meshActor->SetMapper ( meshMapper );
    meshActor->SetVisibility(false);
    meshActor->GetProperty()->SetOpacity(0.8); //Make the mesh have some transparency.

    renderer->AddActor(nodeActor);
    renderer->AddActor(labelActor);
    renderer->AddActor(loopActor);
    renderer->AddActor(areaActor);
    renderer->AddActor(meshActor);
    applyColors();

    //            std::vector<Eigen::Vector2d> clippedPolygonTemp(SutherlandHodgman::clip(polygon2d,box2d));
    
}



void NetworkLoopActor::updateConfiguration(const SlipSystemTab::SlipSystemColorMapType& sscm)
{
    std::cout<<"Updating loops..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    if(dislocationNetwork)
    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkDoubleArray> nodeLabels(vtkSmartPointer<vtkDoubleArray>::New());
        nodeLabels->SetNumberOfComponents(1);
        nodeLabels->SetName("node IDs");

        // Loop occurrences remain distinct even when they share a network node.
        // Record the actual VTK point index; never retain it across frame updates.
        std::map<size_t,vtkIdType> loopNodeIndices;
        nodeSlipSystems.clear();
        for(const auto& weakNode : dislocationNetwork->loopNodes())
        {
            const auto& loopNode(weakNode.second.lock());
            const auto& P(loopNode->get_P());
            loopNodeIndices.emplace(loopNode->sID,points->InsertNextPoint(P.data()));
            nodeLabels->InsertNextTuple1(loopNode->sID);
            nodeSlipSystems.push_back(loopNode->loop()->slipSystem().get());
        }
        
        nodePolyData->SetPoints(points);
//        nodePolyData->GetPointData()->SetScalars(nodeColors);
        nodePolyData->Modified();

        labelPolyData->SetPoints(points);
        labelPolyData->GetPointData()->SetScalars(nodeLabels);
        labelPolyData->Modified();

        
        vtkSmartPointer<vtkCellArray> cells(vtkSmartPointer<vtkCellArray>::New());
        lineSlipSystems.clear();
        areaSlipSystems.clear();
        for(const auto& loopLink : dislocationNetwork->loopLinks())
        {
            vtkSmartPointer<vtkLine> line(vtkSmartPointer<vtkLine>::New());
            line->GetPointIds()->SetId(0, loopNodeIndices.at(loopLink.second.source->sID)); // the second 0 is the index of the Origin in linesPolyData's points
            line->GetPointIds()->SetId(1, loopNodeIndices.at(loopLink.second.sink->sID));
            cells->InsertNextCell(line);
            lineSlipSystems.push_back(loopLink.second.loop->slipSystem().get());


        }
        loopPolyData->SetPoints(points);
        loopPolyData->SetLines(cells);
        loopPolyData->Modified();
        
        // Slipped area
        vtkNew<vtkPoints> areaPoints;
        vtkNew<vtkCellArray> areaPolygons;


        areaPoints->SetDataTypeToDouble();
        for(const auto& weakloop : dislocationNetwork->loops())
        {
            const auto& loop(weakloop.second.lock());
            for(const auto& pair : loop->patches().localPatches())
            {
                for(const auto& triangle : triangulateLoopArea(pair.second))
                {
                    vtkIdType ids[3];
                    for(size_t i=0;i<3;++i)
                    {
                        const Eigen::Vector3d position=pair.first->glidePlane->globalPosition(triangle[i]);
                        ids[i]=areaPoints->InsertNextPoint(position.data());
                    }
                    areaPolygons->InsertNextCell(3,ids);
                    areaSlipSystems.push_back(loop->slipSystem().get());
                }
            }
        }
        areaPolyData->SetPoints(areaPoints);
        areaPolyData->SetPolys(areaPolygons);
        areaPolyData->Modified();
        areaMapper->SetScalarModeToUseCellData();


        
        
        if(meshAreaBox->isChecked())
        {
            vtkSmartPointer<vtkPoints> meshPts(vtkSmartPointer<vtkPoints>::New());
            vtkSmartPointer<vtkCellArray> meshTriangles(vtkSmartPointer<vtkCellArray>::New());
            meshSlipSystems.clear();

            const double meshSize(std::atof(globalMeshSizeEdit->text() .toStdString().c_str()));
            const double localMeshSize(std::atof(localMeshSizeEdit->text() .toStdString().c_str()));
            size_t ptsIncrement(0);
            for(const auto& weakloop : dislocationNetwork->loops())
            {
                const auto& loop(weakloop.second.lock());
                for(const auto& patchMesh : loop->meshed(meshSize,localMeshSize))
                {
                    for(const auto& point3d : patchMesh.points)
                    {
                        meshPts->InsertNextPoint(point3d(0),point3d(1),point3d(2));
                    }
                    for(const auto& tri : patchMesh.triangles)
                    {
                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        triangle->GetPointIds()->SetId (0,tri(0)+ptsIncrement);
                        triangle->GetPointIds()->SetId (1,tri(1)+ptsIncrement);
                        triangle->GetPointIds()->SetId (2,tri(2)+ptsIncrement);
                        meshTriangles->InsertNextCell ( triangle );
                        meshSlipSystems.push_back(loop->slipSystem().get());

                    }
                    ptsIncrement+=patchMesh.points.size();
                }
            }
            meshPolydata->SetPoints ( meshPts );
            meshPolydata->SetPolys ( meshTriangles );
            meshPolydata->Modified();
            meshMapper->SetScalarModeToUseCellData();
//            renWin->Render();

            
        }
//        renderWindow->Render();
    }
    updateColors(sscm);
    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
}

void NetworkLoopActor::applyColors()
{
    const bool manual = colorScheme->currentIndex() == 0;
    manualColorButton->setVisible(manual);
    manualColorButton->setStyleSheet(QString("background-color: %1; color: %2;")
        .arg(manualColor.name(), manualColor.lightness() > 128 ? "black" : "white"));
    const auto apply = [&](vtkPolyData* data, vtkPolyDataMapper* mapper, vtkActor* actor,
                           const std::vector<const SlipSystem*>& systems, bool pointData)
    {
        actor->GetProperty()->SetColor(manualColor.redF(), manualColor.greenF(), manualColor.blueF());
        mapper->SetScalarVisibility(!manual);
        if(!manual)
        {
            auto colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
            colors->SetNumberOfComponents(3);
            colors->SetName("Slip system colors");
            for(const auto* system : systems)
            {
                Eigen::Vector3i color(255, 255, 255);
                const auto found = slipSystemColors.find(system);
                if(found != slipSystemColors.end())
                {
                    color = found->second->rgb();
                }
                colors->InsertNextTuple3(color.x(), color.y(), color.z());
            }
            if(pointData)
            {
                data->GetPointData()->SetScalars(colors);
                mapper->SetScalarModeToUsePointData();
            }
            else
            {
                data->GetCellData()->SetScalars(colors);
                mapper->SetScalarModeToUseCellData();
            }
            data->Modified();
            mapper->SetColorModeToDirectScalars();
        }
    };
    apply(nodePolyData, nodeMapper, nodeActor, nodeSlipSystems, true);
    apply(loopPolyData, loopMapper, loopActor, lineSlipSystems, false);
    apply(areaPolyData, areaMapper, areaActor, areaSlipSystems, false);
    apply(meshPolydata, meshMapper, meshActor, meshSlipSystems, false);

    nodeGlyphs->Update();
}

void NetworkLoopActor::setManualColor(const QColor& color)
{
    if(color.isValid())
    {
        manualColor = color;
        applyColors();
        modify();
    }
}

void NetworkLoopActor::updateColors(const SlipSystemTab::SlipSystemColorMapType& sscm)
{
    slipSystemColors = sscm;
    applyColors();
}

void NetworkLoopActor::modify()
{    
    nodeGlyphs->SetScaleFactor(10);
//    sliderNodeRadius->setEnabled(showNodes->isChecked());
    nodeActor->SetVisibility(showLoops->isChecked());
    labelActor->SetVisibility(showLoops->isChecked());

    
    loopActor->SetVisibility(showLoops->isChecked());
    areaActor->SetVisibility(slippedAreaBox->isChecked());
    meshActor->SetVisibility(meshAreaBox->isChecked());
    areaActor->GetProperty()->SetOpacity(sliderSlippedArea->value()/10.0);
    
    if(renderWindow->GetReadyForRendering())
    {
        renderWindow->Render();
    }
}

} // namespace model
#endif
