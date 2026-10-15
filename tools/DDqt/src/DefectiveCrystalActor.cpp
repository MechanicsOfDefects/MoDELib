/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalActor_cpp_
#define model_DefectiveCrystalActor_cpp_

#include <algorithm>
#include <QString>
#include <QApplication>
#include <QFuture>
#include <QtCore>
#include <QtConcurrent/QtConcurrentRun>
#include <QImage>

#include <DefectiveCrystalActor.h>
#include <SimplicialMesh.h>

#include <DDqtControls.h>
#include <QHBoxLayout>
#include <sstream>

namespace model
{
        
        DefectiveCrystalActor::DefectiveCrystalActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,
                                 QVTKOpenGLStereoWidget* const qvtkGLwidget_in,
                                 DislocationDynamicsBase<3>& ddBase_in) :
        /* init */ DDconfigIO<3>(ddBase_in.simulationParameters.traitsIO.evlFolder)
        /* init */,DDauxIO<3>(ddBase_in.simulationParameters.traitsIO.auxFolder)
        /* init */,renderWindow(renWin)
        /* init */,qvtkGLwidget(qvtkGLwidget_in)
        /* init */,ddBase(ddBase_in)
//        /* init */,configFields(ddBase,*this)
        /* init */,defectiveCrystal(ddBase)
        /* init */,nodes(new NetworkNodeActor(renWin,ren))
        /* init */,segments(new NetworkLinkActor(renWin,ren,defectiveCrystal))
        /* init */,loops(new NetworkLoopActor(renWin,ren,defectiveCrystal))
        /* init */,inclusions(new InclusionActor(renWin,ren,defectiveCrystal))
//        /* init */,glidePlanes(new GlidePlaneActor(renWin,ren,ddBase.poly,ddBase.simulationParameters.traitsIO))
        /* init */,glidePlanes(new GlidePlaneActor(renWin,ren,defectiveCrystal))
        /* init */,quadrature(new QuadratureActor(renWin,ren,ddBase.poly,ddBase.simulationParameters.traitsIO))
        /* init */,chartActor(new ChartActor(ddBase.simulationParameters.traitsIO,renderWindow,ren))
        /* init */,ddField(new DDFieldWidget(renderWindow,ren,defectiveCrystal))
        /* init */,surfaceMesh(new SurfaceMeshActor(renderWindow,ren,defectiveCrystal))
        /* init */,slipSystemTab(new SlipSystemTab(defectiveCrystal))
        /* init */,renderer(ren)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,frameIDedit(new QLineEdit("0"))
        /* init */,plusFrameButton(new QPushButton(">"))
        /* init */,playFrameButton(new QPushButton(">>"))
        /* init */,stopFrameButton(new QPushButton("x"))
        /* init */,minusFrameButton(new QPushButton("<"))
        /* init */,frameIncrementEdit(new QLineEdit(QString::number(ddBase.simulationParameters.outputFrequency)))
        /* init */,rotationAxisEdit(new QLineEdit("0 0 1"))
        /* init */,rotationAngleEdit(new QLineEdit("0"))
        /* init */,saveImage(new QCheckBox(this))
        /* init */,tabWidget(new WrappingTabWidget())
        {

            tabWidget->addTab(nodes, tr(std::string("Nodes").c_str()));
            tabWidget->addTab(segments, tr(std::string("Segments").c_str()));
            tabWidget->addTab(loops, tr(std::string("Loops").c_str()));
            tabWidget->addTab(inclusions, tr(std::string("Inclusions").c_str()));
            tabWidget->addTab(glidePlanes, tr(std::string("GlidePlanes").c_str()));
            tabWidget->addTab(quadrature, tr(std::string("Quadrature").c_str()));
            tabWidget->addTab(chartActor, tr(std::string("Chart").c_str()));
            tabWidget->addTab(ddField, tr(std::string("Fields").c_str()));
            tabWidget->addTab(surfaceMesh, tr("SurfaceMesh"));
            tabWidget->addTab(slipSystemTab, tr(std::string("SlipSystems").c_str()));

            saveImage->setText("save PNG");

            auto *framesBox =
                controls::box("Frames", this, "Select simulation outputs, control playback and capture images.");
            auto *framesLayout = controls::grid(framesBox);
            auto *frameIDLayout = new QGridLayout();
            auto *incrementLayout = new QGridLayout();
            framesLayout->addLayout(frameIDLayout, 0, 0);
            framesLayout->addLayout(incrementLayout, 0, 1);
            framesLayout->setColumnStretch(0, 1);
            framesLayout->setColumnStretch(1, 1);
            controls::field(
                frameIDLayout, 0, "Frame ID", frameIDedit,
                "Simulation output runID to display. Press Enter to load it; output IDs may not be consecutive.");
            controls::field(
                incrementLayout, 0, "Increment", frameIncrementEdit,
                "RunID increment used by previous, next and playback controls. Normally matches outputFrequency.");
            auto *axisLayout = new QGridLayout();
            auto *angleLayout = new QGridLayout();
            framesLayout->addLayout(axisLayout, 1, 0);
            framesLayout->addLayout(angleLayout, 1, 1);
            rotationAxisEdit->setObjectName("frameRotationAxis");
            rotationAngleEdit->setObjectName("frameRotationAngle");
            controls::field(axisLayout, 0, "Axis", rotationAxisEdit,
                            "Global x y z rotation axis through the simulation-box center. Automatic rotation keeps the box centered. Enter a nonzero vector; it is normalized automatically. Applies on Enter or when loading a frame.");
            controls::field(angleLayout, 0, "Angle/frame (°)", rotationAngleEdit,
                            "View rotation in degrees equals this value times the simulation frameID, including skipped IDs. Positive angles follow the right-hand rule. Zero disables rotation. Applies on Enter or when loading a frame.");
            connect(rotationAxisEdit, &QLineEdit::returnPressed, this, &DefectiveCrystalActor::updateRotationSettings);
            connect(rotationAngleEdit, &QLineEdit::returnPressed, this, &DefectiveCrystalActor::updateRotationSettings);
            auto *playback = new QWidget(this);
            auto *playbackLayout = new QHBoxLayout(playback);
            playbackLayout->setContentsMargins(0, 0, 0, 0);
            playbackLayout->addWidget(minusFrameButton);
            playbackLayout->addWidget(plusFrameButton);
            playbackLayout->addWidget(playFrameButton);
            playbackLayout->addWidget(stopFrameButton);
            minusFrameButton->setToolTip("Load the previous output using the frame increment.");
            plusFrameButton->setToolTip("Load the next output using the frame increment.");
            playFrameButton->setToolTip("Advance through available frames until stopped.");
            stopFrameButton->setToolTip("Stop automatic playback.");
            framesLayout->addWidget(playback, 2, 0, 1, 2);
            framesLayout->addWidget(saveImage, 3, 0, 1, 2);
            saveImage->setToolTip("Save a PNG of each loaded frame to the simulation output folder.");
            mainLayout->addWidget(framesBox, 0, 0, 1, 6);
            mainLayout->addWidget(tabWidget,2,0,1,6);
            this->setLayout(mainLayout);

            connect(slipSystemTab, &SlipSystemTab::colorsChanged, this,
                    [this]()
                    {
                        loops->updateColors(slipSystemTab->slipSystemColorMap);
                        if(renderWindow->GetReadyForRendering())
                        {
                            renderWindow->Render();
                        }
                    });

            connect(frameIDedit,SIGNAL(returnPressed()), this, SLOT(updateConfiguration()));
            connect(plusFrameButton,SIGNAL(pressed()), this, SLOT(nextConfiguration()));
            connect(minusFrameButton,SIGNAL(pressed()), this, SLOT(prevConfiguration()));
            connect(playFrameButton,SIGNAL(pressed()), this, SLOT(playConfigurations()));
            connect(stopFrameButton,SIGNAL(pressed()), this, SLOT(stopConfigurations()));

            connect(frameIDedit,SIGNAL(returnPressed()), ddField, SLOT(compute()));
            connect(plusFrameButton,SIGNAL(pressed()), ddField, SLOT(compute()));
            connect(minusFrameButton,SIGNAL(pressed()), ddField, SLOT(compute()));
            connect(playFrameButton,SIGNAL(pressed()), ddField, SLOT(compute()));
            connect(stopFrameButton,SIGNAL(pressed()), ddField, SLOT(compute()));

            
            QApplication::processEvents();
            updateConfiguration(0);
        }

        bool DefectiveCrystalActor::readRotationSettings()
        {
            Eigen::Vector3d axis;
            std::istringstream input(rotationAxisEdit->text().toStdString());
            std::string extra;
            const bool validAxis = static_cast<bool>(input >> axis[0] >> axis[1] >> axis[2]) &&
                                   !(input >> extra) && axis.allFinite() && axis.stableNorm() > 0.0;
            bool validAngle = false;
            const double angle = rotationAngleEdit->text().toDouble(&validAngle);
            validAngle = validAngle && std::isfinite(angle);
            rotationAxisEdit->setStyleSheet(validAxis ? "" : "QLineEdit { border: 1px solid #c44; }");
            rotationAngleEdit->setStyleSheet(validAngle ? "" : "QLineEdit { border: 1px solid #c44; }");
            if(!validAxis || !validAngle)
            {
                return false;
            }
            rotationAxis = axis.stableNormalized();
            rotationDegreesPerFrame = angle;
            return true;
        }

        void DefectiveCrystalActor::updateRotationSettings()
        {
            if(readRotationSettings() && frameLoaded)
            {
                frameRotation.apply(renderer->GetActiveCamera(),
                                    0.5 * (ddBase.mesh.xMin() + ddBase.mesh.xMax()),
                                    rotationAxis, rotationDegreesPerFrame, loadedFrameID);
                renderer->ResetCameraClippingRange();
                if(renderWindow->GetReadyForRendering())
                {
                    renderWindow->Render();
                }
            }
        }

        const DDconfigIO<3>& DefectiveCrystalActor::configIO() const
        {
            return *this;
        }

        DDconfigIO<3>& DefectiveCrystalActor::configIO()
        {
            return *this;
        }

        const DDauxIO<3>& DefectiveCrystalActor::auxIO() const
        {
            return *this;
        }

        DDauxIO<3>& DefectiveCrystalActor::auxIO()
        {
            return *this;
        }

        /*
        The timer's timeout() signal calls the lambda repeatedly
        and each call executes nextConfiguration(). When nextConfiguration()
        returns false, the timer deletes itself (deleteLater)
        */
        void DefectiveCrystalActor::playConfigurations()
        {
            m_stopRequested = false;

            // Get target frame rate (hardcoded, maybe implement a slider in the UI)
            const int fps = 10;
            const int delayMs = 1000 / std::max(1, fps);

            // Timer Setup
            QTimer* timer = new QTimer(this);
            // Lambda Connection
            connect(timer, &QTimer::timeout, [this, timer]() {
                // also, stops the iteration if stop is requested (m_stopRequested)
                if (m_stopRequested || !nextConfiguration()) {
                    timer->stop();
                    timer->deleteLater();
                    return;
                }
                // Let VTK handle rendering at its own pace
                qvtkGLwidget->update(); // Queues a paint event (non-blocking)
                timer->start(delayMs); 
            });
            // Triggers timeout() signal after event processing and wait for delayMs
            timer->start(delayMs); 
        }

        /*
        stop animation play button
        private class definition is added in the header file
        */
        void DefectiveCrystalActor::stopConfigurations()
        {
            m_stopRequested = true;
        }

        bool DefectiveCrystalActor::nextConfiguration()
        {
            const long int currentFrameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            const long int currentIncrement(std::atoi(frameIncrementEdit->text() .toStdString().c_str()));
            const long int nextFrameID(std::max((long int)0,currentFrameID+currentIncrement));
            frameIDedit->setText(QString::fromStdString(std::to_string(nextFrameID)));
            const bool updated(updateConfiguration());
            if(!updated)
            {
                frameIDedit->setText(QString::fromStdString(std::to_string(currentFrameID)));
            }
            return updated;
        }

        bool DefectiveCrystalActor::prevConfiguration()
        {
            const long int currentFrameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            const long int currentIncrement(std::atoi(frameIncrementEdit->text() .toStdString().c_str()));
            const long int nextFrameID(std::max((long int)0,currentFrameID-currentIncrement));
            frameIDedit->setText(QString::fromStdString(std::to_string(nextFrameID)));
            const bool updated(updateConfiguration());
            if(!updated)
            {
                frameIDedit->setText(QString::fromStdString(std::to_string(currentFrameID)));
            }
            return updated;
        }

        bool DefectiveCrystalActor::updateConfiguration()
        {
            const size_t frameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            if(updateConfiguration(frameID))
            {
                renderWindow->Render();
                return true;
            }
            else
            {
                return false;
            }
        }

        bool DefectiveCrystalActor::updateConfiguration(const size_t& frameID)
        {
            try
            {
                DislocationNode<3>::force_count(0);
                DislocationLoopNode<3>::force_count(0);
                DislocationLoop<3>::force_count(0);
                EshelbyInclusionBase<3>::force_count(0);

                configIO().read(frameID);
                auxIO().read(frameID);
                //configFields.updateConfiguration();
                defectiveCrystal.initializeConfiguration(*this);
                nodes->updateConfiguration(*this);
                segments->updateConfiguration(nodes->nodePolyData,slipSystemTab->slipSystemColorMap);
                loops->updateConfiguration(slipSystemTab->slipSystemColorMap);
                inclusions->updateConfiguration();
                surfaceMesh->updateConfiguration();
//                glidePlanes->updateConfiguration(*this);
                glidePlanes->updateConfiguration();
                quadrature->updateConfiguration(*this);
                chartActor->updateConfiguration(frameID);
//                ddField->compute(); // this causes segFault at DDFieldWidget::plotField::renWin->Render();
                
                loadedFrameID = frameID;
                frameLoaded = true;
                // Navigation must use the visible values, even when the user has
                // not pressed Enter. Invalid edits retain the last valid settings.
                readRotationSettings();
                frameRotation.apply(renderer->GetActiveCamera(),
                                    0.5 * (ddBase.mesh.xMin() + ddBase.mesh.xMax()),
                                    rotationAxis, rotationDegreesPerFrame, frameID);
                renderer->ResetCameraClippingRange();
                if(saveImage->isChecked())
                {
                    renderWindow->Render();
                    //qDebug() << "Saving frame:" << frameID 
                    //<< "| UI shows:" << frameIDedit->text();
                    //qvtkGLwidget->repaint(); // Force render update
                    //QCoreApplication::processEvents(); // Ensure paint completes
                    QImage img=qvtkGLwidget->grabFramebuffer();
                    img.save(QString::fromStdString(ddBase.simulationParameters.traitsIO.evlFolder+"/img_"+std::to_string(frameID)+".png"), "PNG", -1);
                }
                return true;
            }
            catch(const std::exception& e)
            {
                std::cout<<e.what()<<std::endl;
                return false;
            }
        }
        
        void DefectiveCrystalActor::modify()
        {
//            nodes.modify();
//            segments.modify();
//            loops.modify();
        }
        
}
#endif
