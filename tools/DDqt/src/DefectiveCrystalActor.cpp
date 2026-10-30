#include <DislocationDynamics/DislocationDynamicsBase.h>
/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalActor_cpp_
#define model_DefectiveCrystalActor_cpp_

#include <DislocationNetwork.h>
#include <algorithm>
#include <QString>
#include <QApplication>

#include <QtCore>

#include <QImage>
#include <QStyle>
#include <QScopedValueRollback>
#include <limits>

#include <DefectiveCrystalActor.h>
#include <SimplicialMesh.h>

#include <DDqtControls.h>
#include <QHBoxLayout>
#include <sstream>

namespace model
{
        
        DefectiveCrystalActor::DefectiveCrystalActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,
                                 QVTKOpenGLStereoWidget* const qvtkGLwidget_in,
                                 DefectiveCrystalBase<3>& ddBase_in) :
        /* init */ renderWindow(renWin)
        /* init */,qvtkGLwidget(qvtkGLwidget_in)
        /* init */,ddBase(ddBase_in)
//        /* init */,configFields(ddBase,*this)
        /* init */,defectiveCrystal(ddBase)
        /* init */,nodes(new NetworkNodeActor(renWin,ren,(ddBase.mesh.xMax()-ddBase.mesh.xMin()).norm()))
        /* init */,segments(new NetworkLinkActor(renWin,ren,defectiveCrystal))
        /* init */,loops(new NetworkLoopActor(renWin,ren,defectiveCrystal))
        /* init */,inclusions(new InclusionActor(renWin,ren,defectiveCrystal))
//        /* init */,glidePlanes(new GlidePlaneActor(renWin,ren,ddBase.poly,ddBase.simulationParameters.traitsIO))
        /* init */,glidePlanes(new GlidePlaneActor(renWin,ren,defectiveCrystal))
        /* init */,quadrature(new QuadratureActor(renWin,ren,ddBase.poly,ddBase.simulationParameters.traitsIO))
        /* init */,chartActor(new ChartActor(ddBase.simulationParameters.traitsIO,renderWindow,ren))
        /* init */,ddField(new DDFieldWidget(renderWindow,ren,defectiveCrystal))
        /* init */,surfaceMesh(new SurfaceMeshActor(renderWindow,ren,defectiveCrystal))
        /* init */,boxes(new BoxTreeActor(renderWindow,ren,ddBase.fmmContext->boxTree,
            [base = ddBase.dislocationDynamicsBase]()
            {
                std::vector<std::shared_ptr<GlidePlane<3>>> planes;
                if(base)
                {
                    for(const auto& [key, weakPlane] : base->glidePlaneFactory.glidePlanes())
                    {
                        if(const auto plane = weakPlane.lock())
                        {
                            planes.push_back(plane);
                        }
                    }
                }
                return planes;
            }))
        /* init */,slipSystemTab(new SlipSystemTab(defectiveCrystal))
        /* init */,configuration(std::make_unique<DefectiveCrystalConfigIO<3>>(ddBase.simulationParameters.traitsIO.configFolder))
        /* init */,auxiliary(std::make_unique<DefectiveCrystalAuxIO<3>>(ddBase.simulationParameters.traitsIO.auxFolder))
        /* init */,playbackTimer(new QTimer(this))
        /* init */,renderer(ren)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,frameIDedit(new QLineEdit("0"))
        /* init */,plusFrameButton(new QToolButton(this))
        /* init */,playFrameButton(new QToolButton(this))
        /* init */,stopFrameButton(new QToolButton(this))
        /* init */,minusFrameButton(new QToolButton(this))
        /* init */,frameIncrementEdit(new QLineEdit(QString::number(ddBase.simulationParameters.outputFrequency)))
        /* init */,rotationAxisEdit(new QLineEdit("0 0 1"))
        /* init */,rotationAngleEdit(new QLineEdit("0"))
        /* init */,cameraAdvanceEdit(new QLineEdit("0"))
        /* init */,saveImage(new QCheckBox(this))
        /* init */,frameStatus(new QLabel(this))
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
            tabWidget->addTab(boxes, tr("Boxes"));
            tabWidget->addTab(slipSystemTab, tr(std::string("SlipSystems").c_str()));

            saveImage->setText("save PNG");

            auto *framesBox =
                controls::box("Navigation", this, "Select simulation outputs, control playback and capture images.");
            auto *framesLayout = controls::grid(framesBox);
            auto *frameIDLayout = new QGridLayout();
            auto *incrementLayout = new QGridLayout();
            framesLayout->addLayout(frameIDLayout, 0, 0);
            framesLayout->addLayout(incrementLayout, 0, 1);
            framesLayout->setColumnStretch(0, 1);
            framesLayout->setColumnStretch(1, 1);
            controls::field(
                frameIDLayout, 0, "Frame ID", frameIDedit,
                "Simulation output runID to display. Press Enter or leave this field to load it and apply camera rotation and advance; output IDs may not be consecutive.");
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
            auto *advanceLayout = new QGridLayout();
            framesLayout->addLayout(advanceLayout, 2, 0, 1, 2);
            cameraAdvanceEdit->setObjectName("frameCameraAdvance");
            controls::field(advanceLayout, 0, "Advance/frame (b)", cameraAdvanceEdit,
                            "Camera travel toward and through the box center, in simulation length units b times frameID, including skipped IDs. Positive advances, negative retreats, zero disables travel. Combines with rotation and uses perspective projection. Applies on Enter or when loading a frame.");
            connect(cameraAdvanceEdit, &QLineEdit::returnPressed, this, &DefectiveCrystalActor::updateRotationSettings);
            auto *playback = new QWidget(this);
            auto *playbackLayout = new QHBoxLayout(playback);
            playbackLayout->setContentsMargins(0, 0, 0, 0);
            minusFrameButton->setIcon(style()->standardIcon(QStyle::SP_MediaSkipBackward));
            playFrameButton->setIcon(style()->standardIcon(QStyle::SP_MediaPlay));
            stopFrameButton->setIcon(style()->standardIcon(QStyle::SP_MediaStop));
            plusFrameButton->setIcon(style()->standardIcon(QStyle::SP_MediaSkipForward));
            minusFrameButton->setAccessibleName("Previous frame");
            playFrameButton->setAccessibleName("Play");
            stopFrameButton->setAccessibleName("Stop");
            plusFrameButton->setAccessibleName("Next frame");
            for(auto *button : {minusFrameButton, plusFrameButton, playFrameButton, stopFrameButton})
            {
                button->setToolButtonStyle(Qt::ToolButtonIconOnly);
                button->setIconSize(QSize(20, 20));
                button->setFixedSize(32, 32);
                if(button == playFrameButton)
                {
                    playbackLayout->addSpacing(12);
                }
                playbackLayout->addWidget(button);
            }
            playbackLayout->addStretch();
            minusFrameButton->setToolTip("Load the previous output using the frame increment.");
            plusFrameButton->setToolTip("Load the next output using the frame increment.");
            playFrameButton->setToolTip("Advance through available frames until stopped.");
            stopFrameButton->setToolTip("Stop automatic playback.");
            framesLayout->addWidget(playback, 3, 0, 1, 2);
            framesLayout->addWidget(saveImage, 4, 0, 1, 2);
            frameStatus->setObjectName("frameStatus");
            frameStatus->setWordWrap(true);
            framesLayout->addWidget(frameStatus, 5, 0, 1, 2);
            saveImage->setToolTip("Save a PNG of each loaded frame to the simulation output folder.");
            mainLayout->addWidget(framesBox, 0, 0, 1, 6);
            periodicImages = new PeriodicImagesWidget(renderWindow, renderer, ddBase.mesh.periodicBasis(), this);
            mainLayout->addWidget(periodicImages, 1, 0, 1, 6);
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

            connect(frameIDedit, &QLineEdit::editingFinished, this, [this]()
            {
                if(frameIDedit->isModified())
                {
                    frameIDedit->setModified(false);
                    updateConfiguration();
                }
            });
            connect(plusFrameButton,SIGNAL(pressed()), this, SLOT(nextConfiguration()));
            connect(minusFrameButton,SIGNAL(pressed()), this, SLOT(prevConfiguration()));
            connect(playFrameButton,SIGNAL(pressed()), this, SLOT(playConfigurations()));
            connect(stopFrameButton,SIGNAL(pressed()), this, SLOT(stopConfigurations()));

            playbackTimer->setObjectName("framePlaybackTimer");
            playbackTimer->setSingleShot(true);
            playbackTimer->setInterval(100);
            stopFrameButton->setEnabled(false);
            connect(playbackTimer, &QTimer::timeout, this, [this]()
            {
                if(playing && nextConfiguration() && playing)
                {
                    playbackTimer->start();
                }
                else
                {
                    stopConfigurations();
                }
            });

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
            bool validAdvance = false;
            const double advance = cameraAdvanceEdit->text().toDouble(&validAdvance);
            validAdvance = validAdvance && std::isfinite(advance);
            cameraAdvanceEdit->setStyleSheet(validAdvance ? "" : "QLineEdit { border: 1px solid #c44; }");
            if(!validAxis || !validAngle || !validAdvance)
            {
                return false;
            }
            rotationAxis = axis.stableNormalized();
            rotationDegreesPerFrame = angle;
            cameraAdvancePerFrame = advance;
            return true;
        }

        void DefectiveCrystalActor::updateRotationSettings()
        {
            if(readRotationSettings() && frameLoaded)
            {
                frameRotation.apply(renderer->GetActiveCamera(),
                                    0.5 * (ddBase.mesh.xMin() + ddBase.mesh.xMax()),
                                    rotationAxis, rotationDegreesPerFrame, loadedFrameID, cameraAdvancePerFrame);
                renderer->ResetCameraClippingRange();
                if(renderWindow->GetReadyForRendering())
                {
                    renderWindow->Render();
                }
            }
        }

        const DefectiveCrystalConfigIO<3>& DefectiveCrystalActor::configIO() const
        {
            return *configuration;
        }

        DefectiveCrystalConfigIO<3>& DefectiveCrystalActor::configIO()
        {
            return *configuration;
        }

        const DefectiveCrystalAuxIO<3>& DefectiveCrystalActor::auxIO() const
        {
            return *auxiliary;
        }

        DefectiveCrystalAuxIO<3>& DefectiveCrystalActor::auxIO()
        {
            return *auxiliary;
        }

        DefectiveCrystalActor::~DefectiveCrystalActor()
        {
            stopConfigurations();
            // Actors and field samplers reference defectiveCrystal and ddBase.
            // Destroy them before those members and their data are released.
            delete periodicImages;
            delete tabWidget;
        }

        void DefectiveCrystalActor::playConfigurations()
        {
            if(playing)
            {
                return;
            }
            playing = true;
            playFrameButton->setEnabled(false);
            stopFrameButton->setEnabled(true);
            playbackTimer->start();
        }

        void DefectiveCrystalActor::stopConfigurations()
        {
            playing = false;
            playbackTimer->stop();
            playFrameButton->setEnabled(true);
            stopFrameButton->setEnabled(false);
        }

        bool DefectiveCrystalActor::nextConfiguration()
        {
            bool ok = false;
            const auto increment = frameIncrementEdit->text().toULongLong(&ok);
            const size_t current = frameLoaded ? loadedFrameID : 0;
            if(!ok || increment == 0 || increment > std::numeric_limits<size_t>::max() - current)
            {
                frameStatus->setText("Enter a positive frame increment within the runID range.");
                stopConfigurations();
                return false;
            }
            return updateConfiguration(current + static_cast<size_t>(increment));
        }

        bool DefectiveCrystalActor::prevConfiguration()
        {
            bool ok = false;
            const auto increment = frameIncrementEdit->text().toULongLong(&ok);
            if(!ok || increment == 0)
            {
                frameStatus->setText("Enter a positive frame increment.");
                stopConfigurations();
                return false;
            }
            const size_t current = frameLoaded ? loadedFrameID : 0;
            return updateConfiguration(increment > current ? 0 : current - static_cast<size_t>(increment));
        }

        bool DefectiveCrystalActor::updateConfiguration()
        {
            bool ok = false;
            const auto requested = frameIDedit->text().toULongLong(&ok);
            if(!ok || requested > std::numeric_limits<size_t>::max())
            {
                frameStatus->setText("Enter a nonnegative integer frame ID.");
                frameIDedit->setText(QString::number(loadedFrameID));
                stopConfigurations();
                return false;
            }
            return updateConfiguration(static_cast<size_t>(requested));
        }

        void DefectiveCrystalActor::refreshFrame(size_t frameID)
        {
            DislocationNode<3>::force_count(0);
            DislocationLoopNode<3>::force_count(0);
            DislocationLoop<3>::force_count(0);
            EshelbyInclusionBase<3>::force_count(0);
            defectiveCrystal.initializeConfiguration(configIO());
            nodes->updateConfiguration(configIO(), loops->dislocationNetwork.get());
            segments->updateConfiguration(nodes->nodePolyData, slipSystemTab->slipSystemColorMap);
            loops->updateConfiguration(slipSystemTab->slipSystemColorMap);
            inclusions->updateConfiguration();
            surfaceMesh->updateConfiguration();
            boxes->updateConfiguration();
            glidePlanes->updateConfiguration();
            quadrature->updateConfiguration(auxIO());
            chartActor->updateConfiguration(frameID);
            ddField->compute();
        }

        bool DefectiveCrystalActor::updateConfiguration(const size_t& frameID)
        {
            if(updatingFrame)
            {
                return false;
            }
            QScopedValueRollback<bool> updating(updatingFrame, true);
            auto nextConfig = std::make_unique<DefectiveCrystalConfigIO<3>>(ddBase.simulationParameters.traitsIO.configFolder);
            auto nextAux = std::make_unique<DefectiveCrystalAuxIO<3>>(ddBase.simulationParameters.traitsIO.auxFolder);
            try
            {
                // Read both files before touching the displayed frame. DefectiveCrystalAuxIO's
                // missing-file path asserts, so check availability explicitly.
                nextConfig->read(frameID);
                if(!nextAux->isBinGood(frameID) && !nextAux->isTxtGood(frameID))
                {
                    throw std::runtime_error("No auxiliary file for frame " + std::to_string(frameID));
                }
                nextAux->read(frameID);
            }
            catch(const std::exception& error)
            {
                frameStatus->setText(QString::fromStdString(error.what()));
                frameIDedit->setText(QString::number(loadedFrameID));
                stopConfigurations();
                return false;
            }

            const bool ready = renderWindow->GetReadyForRendering();
            renderWindow->SetReadyForRendering(false);
            configuration.swap(nextConfig);
            auxiliary.swap(nextAux);
            try
            {
                refreshFrame(frameID);
            }
            catch(const std::exception& error)
            {
                configuration.swap(nextConfig);
                auxiliary.swap(nextAux);
                QString message = QString::fromStdString(error.what());
                try
                {
                    refreshFrame(loadedFrameID);
                }
                catch(const std::exception& restoreError)
                {
                    message += " Restoration failed: " + QString::fromStdString(restoreError.what());
                }
                renderWindow->SetReadyForRendering(ready);
                frameIDedit->setText(QString::number(loadedFrameID));
                frameStatus->setText(message);
                stopConfigurations();
                return false;
            }

            // Commit the frame only after every actor and Fields has been updated.
            loadedFrameID = frameID;
            frameLoaded = true;
            frameIDedit->setText(QString::number(frameID));
            readRotationSettings();
            frameRotation.apply(renderer->GetActiveCamera(),
                                0.5 * (ddBase.mesh.xMin() + ddBase.mesh.xMax()),
                                rotationAxis, rotationDegreesPerFrame, frameID, cameraAdvancePerFrame);
            renderer->ResetCameraClippingRange();
            renderWindow->SetReadyForRendering(ready);
            frameStatus->clear();
            if(ready)
            {
                renderWindow->Render();
            }
            // Capture only the completed frame, after Fields and camera updates.
            // A capture failure does not undo an otherwise valid frame load.
            if(saveImage->isChecked())
            {
                const QString filename = QString::fromStdString(ddBase.simulationParameters.traitsIO.configFolder +
                                                               "/img_" + std::to_string(frameID) + ".png");
                if(!ready || !qvtkGLwidget || !qvtkGLwidget->grabFramebuffer().save(filename, "PNG"))
                {
                    frameStatus->setText("Frame loaded, but PNG capture failed.");
                    stopConfigurations();
                }
            }
            emit frameChanged(static_cast<qulonglong>(frameID));
            return true;
        }

        void DefectiveCrystalActor::modify()
        {
//            nodes.modify();
//            segments.modify();
//            loops.modify();
        }
        
}
#endif
