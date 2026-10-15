/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalActor_h_
#define model_DefectiveCrystalActor_h_

#include <QGroupBox>
#include <QGridLayout>
#include <QCheckBox>
#include <QLabel>
#include <WrappingTabWidget.h>
#include <QTextEdit>
#include <QLineEdit>
#include <QToolButton>
#include <QVTKOpenGLStereoWidget.h>

#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>

#include <DislocationDynamicsBase.h>
//#include <ConfigurationFields.h>
#include <DefectiveCrystal.h>
#include <DDauxIO.h>
#include <DDFieldWidget.h>
#include <SurfaceMeshActor.h>
#include <FrameCameraRotation.h>
#include <PeriodicImagesWidget.h>
#include <SimplicialMesh.h>
#include <NetworkNodeActor.h>
#include <NetworkLinkActor.h>
#include <NetworkLoopActor.h>
#include <InclusionActor.h>
#include <GlidePlaneActor.h>
#include <QuadratureActor.h>
#include <ChartActor.h>
#include <SlipSystemTab.h>
#include <QCheckBox>

// added std library
#include <memory>
#include <QTimer>

namespace model
{
    struct DefectiveCrystalActor : public QWidget
    {
        
        Q_OBJECT
        
        
    public:
        vtkGenericOpenGLRenderWindow* const renderWindow;
        QVTKOpenGLStereoWidget* const qvtkGLwidget;
//        const DDtraitsIO& traitsIO;
//        SlipSystemTab::SlipSystemColorMapType slipSystemColorMap;
        DislocationDynamicsBase<3>& ddBase;
        DefectiveCrystal<3> defectiveCrystal;
        NetworkNodeActor* nodes;
        NetworkLinkActor* segments;
        NetworkLoopActor* loops;
        InclusionActor* inclusions;
        GlidePlaneActor* glidePlanes;
        QuadratureActor* quadrature;
        ChartActor* chartActor;
        DDFieldWidget* ddField;
        SurfaceMeshActor* surfaceMesh;
        SlipSystemTab* slipSystemTab;
        
    private:
        std::unique_ptr<DDconfigIO<3>> configuration;
        std::unique_ptr<DDauxIO<3>> auxiliary;
        QTimer *const playbackTimer;
        bool playing = false;
        bool updatingFrame = false;
        void refreshFrame(size_t frameID);
        vtkRenderer *const renderer;
        PeriodicImagesWidget* periodicImages;
        FrameCameraRotation frameRotation;
        Eigen::Vector3d rotationAxis = Eigen::Vector3d::UnitZ();
        double rotationDegreesPerFrame = 0.0;
        size_t loadedFrameID = 0;
        bool frameLoaded = false;
        bool readRotationSettings();
        void updateRotationSettings();

        private slots:
        bool updateConfiguration();
        bool nextConfiguration();
        bool prevConfiguration();
        void playConfigurations();

    public slots:
        void stopConfigurations();

    signals:
        void frameChanged(qulonglong frameID);

        public:
        
        QGridLayout* mainLayout;
        QLineEdit* frameIDedit;
        QToolButton* plusFrameButton;
        QToolButton* playFrameButton;
        QToolButton* stopFrameButton;
        QToolButton* minusFrameButton;
        QLineEdit* frameIncrementEdit;
        QLineEdit* rotationAxisEdit;
        QLineEdit* rotationAngleEdit;
        QCheckBox* saveImage;
        QLabel* frameStatus;
        WrappingTabWidget* tabWidget;

        DefectiveCrystalActor(vtkGenericOpenGLRenderWindow* const,
                    vtkRenderer* const ren,
                    QVTKOpenGLStereoWidget* const qvtkGLwidget_in,
                    DislocationDynamicsBase<3>& ddBase_in);
        ~DefectiveCrystalActor() override;
        bool updateConfiguration(const size_t& frameID);
        void modify();
        const DDconfigIO<3>& configIO() const;
        DDconfigIO<3>& configIO();
        const DDauxIO<3>& auxIO() const;
        DDauxIO<3>& auxIO();

    };
    
}
#endif
