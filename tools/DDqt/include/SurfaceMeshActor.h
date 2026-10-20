/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_SurfaceMeshActor_H_
#define model_SurfaceMeshActor_H_
#include <QWidget>
#include <QGroupBox>
#include <QCheckBox>
#include <QLineEdit>
#include <QSlider>
#include <QDoubleSpinBox>
#include <QLabel>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <vtkActor.h>
#include <vtkLight.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <DefectiveCrystal/DefectiveCrystal.h>
#include <SurfaceMeshGeometry.h>

namespace model
{
    class SurfaceMeshActor : public QWidget
    {
        vtkGenericOpenGLRenderWindow *const renderWindow;
        vtkRenderer *const renderer;
        const DefectiveCrystal<3> &defectiveCrystal;
        SurfaceMeshGeometry referenceGeometry;
        SurfaceMeshGeometry geometry;
        std::vector<Eigen::Vector3d> displacement;
        bool frameAvailable = false;
        bool dirty = true;
        bool validDisplacement = false;
        bool clippingEnabled = false;
        Eigen::Vector3d appliedClipCenter = Eigen::Vector3d::Zero();
        Eigen::Vector3d appliedClipNormal = Eigen::Vector3d::UnitX();
        double clipTravel = 0.0;

        std::map<size_t, QCheckBox *> externalFaceChecks;
        std::map<size_t, QCheckBox *> internalFaceChecks;
        std::set<size_t> selectedFaceIDs() const;
        void compute();
        void applyDisplacement();
        void render();
        void updateLighting();
        void moveClipPlane();
        void reportError(const std::exception &error);

      public:
        QCheckBox *const showSurface;
        QLineEdit *const meshSizeEdit;
        QDoubleSpinBox *const amplification;
        QCheckBox *const showEdges;
        QLabel *const status;
        QCheckBox *const clipSurface;
        QLineEdit *const clipCenterEdit;
        QLineEdit *const clipNormalEdit;
        QLabel *const clipStatus;
        QSlider *const clipPositionSlider;
        QLabel *const clipOffsetLabel;
        QGroupBox *const externalFaces;
        QGroupBox *const internalFaces;
        std::vector<QCheckBox *> microstructuresCheck;
        const vtkSmartPointer<vtkPolyData> polyData;
        const vtkSmartPointer<vtkPlane> clipPlane;
        const vtkSmartPointer<vtkClipPolyData> clipper;
        const vtkSmartPointer<vtkCleanPolyData> shadingMesh;
        const vtkSmartPointer<vtkPolyDataNormals> normals;
        const vtkSmartPointer<vtkPolyDataMapper> mapper;
        const vtkSmartPointer<vtkActor> actor;
        const vtkSmartPointer<vtkLight> reliefLight;

        SurfaceMeshActor(vtkGenericOpenGLRenderWindow *window, vtkRenderer *renderer, const DefectiveCrystal<3> &crystal);
        ~SurfaceMeshActor() override;
        void rebuild();
        void updateClipping();
        // Called after initializeConfiguration for every successful frame load, including playback.
        // Rendering is left to the owning viewer, avoiding reentrant frame-load rendering.
        void updateConfiguration();
    };
} // namespace model
#endif
