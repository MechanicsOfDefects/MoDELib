/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_BoxTreeActor_H_
#define model_BoxTreeActor_H_

#include <FMM/BoxTree.h>
#include <FMM/PlaneBoxIntersection.h>
#include <functional>
#include <QWidget>
#include <vtkSmartPointer.h>
#include <memory>
#include <cstdint>

class QGroupBox;
class QCheckBox;
class QComboBox;
class QSpinBox;
class QSlider;
class QDoubleSpinBox;
class QLabel;
class vtkGenericOpenGLRenderWindow;
class vtkRenderer;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkActor;

namespace model
{
    template<int dim> struct GlidePlane;
    class BoxTreeActor : public QWidget
    {
        Q_OBJECT
      public:
        using PlaneProvider = std::function<std::vector<std::shared_ptr<GlidePlane<3>>>()>;
      private:
        const PlaneProvider planeProvider;
        QCheckBox* cutsEdit;
        QComboBox* planeEdit;
        std::map<size_t, std::weak_ptr<GlidePlane<3>>> availablePlanes;
        std::shared_ptr<const PlaneBoxIntersection<3>> displayedCut;
        vtkSmartPointer<vtkPolyData> cutPolydata;
        vtkSmartPointer<vtkPolyDataMapper> cutMapper;
        vtkSmartPointer<vtkActor> cutActor;
        void refreshCuts();
        vtkGenericOpenGLRenderWindow* const renderWindow;
        vtkRenderer* const renderer;
        const std::shared_ptr<BoxTree<3>> tree;
        QGroupBox* group;
        QLabel* depthEdit;
        QSpinBox* levelEdit;
        QSlider* opacityEdit;
        QDoubleSpinBox* widthEdit;
        QLabel* summary;
        QLabel* status;
        vtkSmartPointer<vtkPolyData> polydata;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;
        std::uint64_t displayedRevision = 0;
        int displayedLevel = -1;
        void modify();
        void refreshGeometry();

      public:
        BoxTreeActor(vtkGenericOpenGLRenderWindow*, vtkRenderer*, std::shared_ptr<BoxTree<3>>, PlaneProvider = {});
        ~BoxTreeActor() override;
        void updateConfiguration();
        vtkPolyData* data() const;
        vtkPolyData* cutData() const;
    };
}
#endif
