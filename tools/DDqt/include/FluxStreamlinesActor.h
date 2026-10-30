#ifndef model_FluxStreamlinesActor_H_
#define model_FluxStreamlinesActor_H_
#include <QWidget>
#include <vtkSmartPointer.h>
#include <DefectiveCrystal/DefectiveCrystal.h>
class QGroupBox;
class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QSpinBox;
class QLabel;
class QLineEdit;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkActor;
class vtkRenderer;
class vtkGenericOpenGLRenderWindow;
namespace model
{
    class FluxStreamlinesActor : public QWidget
    {
        const DefectiveCrystal<3>& crystal;
        vtkGenericOpenGLRenderWindow* window;
        vtkRenderer* renderer;
        QGroupBox* group;
        QCheckBox* logColors;
        void updateColors();
        QComboBox *species,*direction;
        QDoubleSpinBox *spacing,*offset,*step,*length,*width;
        QLineEdit* minimum;
        QSpinBox *seedLimit,*stepLimit;
        QLabel* status;
        vtkSmartPointer<vtkPolyData> data;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;
    public:
        FluxStreamlinesActor(vtkGenericOpenGLRenderWindow*,vtkRenderer*,const DefectiveCrystal<3>&);
        ~FluxStreamlinesActor();
        void compute();
        void invalidate();
    };
}
#endif
