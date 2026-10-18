/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDqtVTKwidget_cpp_
#define model_DDqtVTKwidget_cpp_

#include <QFileDialog>
#include <QCoreApplication>
#include <QDir>
#include <QScrollArea>
#include <QVTKInteractor.h>
#include <vtkInteractorStyle.h>
#include <vtkInteractorStyleSwitch.h>

#include <DDqtVTKwidget.h>
#include <DefectiveCrystalActor.h>

namespace model
{
    
    DDqtVTKwidget::DDqtVTKwidget(QWidget*) :
    /* init */ mainLayout(new QGridLayout(this))
    /* init */,tabWidget(new QTabWidget(this))
    /* init */,openglWidget(new QVTKOpenGLStereoWidget(this))
    /* init */,renderWindow(vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New())
    /* init */,renderer(vtkSmartPointer<vtkRenderer>::New())
    /* init */,ddBase(getWorkingDir())
    /* init */,workingDirLabel(new QLabel(QString::fromStdString(ddBase.simulationParameters.traitsIO.simulationFolder)))
    /* init */,meshActor(new SimplicialMeshActor(renderWindow,renderer,ddBase.mesh))
    /* init */,defectiveCrystalActor(new DefectiveCrystalActor(renderWindow,renderer,openglWidget,ddBase))
    {
        renderer->SetBackground(1,1,1);
        renderWindow->AddRenderer(renderer);
        renderWindow->SetMultiSamples(4);
        
        openglWidget->setRenderWindow(renderWindow);
        auto* cameraStyle = vtkInteractorStyle::SafeDownCast(openglWidget->interactor()->GetInteractorStyle());
        if(auto* styleSwitch = vtkInteractorStyleSwitch::SafeDownCast(cameraStyle))
        {
            cameraStyle = styleSwitch->GetCurrentStyle();
        }
        if(cameraStyle)
        {
            // Use smaller zoom steps for precise framing; leave rotation and pan unchanged.
            cameraStyle->SetMouseWheelMotionFactor(0.25);
        }

        
        tabWidget->addTab(defectiveCrystalActor, tr(std::string("Defects").c_str()));
        auto *meshScroll = new QScrollArea(this);
        meshScroll->setWidgetResizable(true);
        meshScroll->setFrameShape(QFrame::NoFrame);
        meshScroll->setWidget(meshActor);
        tabWidget->addTab(meshScroll, tr("Mesh"));

        mainLayout->addWidget(workingDirLabel,0,0,1,2);
        mainLayout->addWidget(tabWidget,1,0,1,1);
        mainLayout->addWidget(openglWidget,1,1,1,1);
        mainLayout->setColumnStretch(0, 3);
        mainLayout->setColumnStretch(1, 7);
        this->setLayout(mainLayout);
        defectiveCrystalActor->updateConfiguration(0);
    }
    

    DDqtVTKwidget::~DDqtVTKwidget()
    {
        defectiveCrystalActor->stopConfigurations();
        delete tabWidget; // Destroy actors while ddBase and renderer still exist.
        renderer->RemoveAllViewProps();
        delete openglWidget;
    }

    std::string DDqtVTKwidget::getWorkingDir() 
    {
        QDir directory(QCoreApplication::applicationDirPath());
        QString initialDirectory=directory.absolutePath();
        // Find tutorials from the executable location, including build/tools/DDqt
        // and application-bundle layouts, independently of the launch directory.
        do
        {
            const QDir tutorials(directory.filePath("tutorials"));
            if(tutorials.exists())
            {
                initialDirectory=tutorials.absolutePath();
                break;
            }
        }
        while(directory.cdUp());

        return QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                                 initialDirectory,
                                                 QFileDialog::ShowDirsOnly
                                                 | QFileDialog::DontUseNativeDialog).toStdString();

    
    }


} // namespace model

#endif







