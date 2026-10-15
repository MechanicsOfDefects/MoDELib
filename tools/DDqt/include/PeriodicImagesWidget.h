/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_PeriodicImagesWidget_H_
#define model_PeriodicImagesWidget_H_

#include <DDqtControls.h>
#include <QSpinBox>
#include <Eigen/Dense>
#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkMatrix4x4.h>
#include <vtkPropCollection.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <array>
#include <map>
#include <set>
#include <vector>

namespace model
{
    // Synchronize immediately before rendering, including actors added or removed
    // by Fields controls. Each copy shares its source mapper and properties.
    class PeriodicImages
    {
        struct Images
        {
            vtkSmartPointer<vtkActor> source;
            std::vector<vtkSmartPointer<vtkActor>> copies;
        };
        vtkRenderer* renderer;
        vtkSmartPointer<vtkCallbackCommand> callback;
        unsigned long observer = 0;
        std::map<vtkActor*, Images> images;
        std::vector<Eigen::Vector3d> translations;

      public:
        explicit PeriodicImages(vtkRenderer* rendererIn) : renderer(rendererIn),
            callback(vtkSmartPointer<vtkCallbackCommand>::New())
        {
            callback->SetClientData(this);
            callback->SetCallback([](vtkObject*, unsigned long, void* data, void*)
            {
                static_cast<PeriodicImages*>(data)->synchronize();
            });
            observer = renderer->AddObserver(vtkCommand::StartEvent, callback);
        }

        ~PeriodicImages()
        {
            renderer->RemoveObserver(observer);
            clear();
        }
        PeriodicImages(const PeriodicImages&) = delete;
        PeriodicImages& operator=(const PeriodicImages&) = delete;

        void clear()
        {
            for(auto& entry : images)
            {
                for(auto& copy : entry.second.copies)
                {
                    renderer->RemoveActor(copy);
                }
            }
            images.clear();
        }

        void setTranslations(std::vector<Eigen::Vector3d> values)
        {
            clear();
            translations = std::move(values);
            synchronize();
        }

        void synchronize()
        {
            if(translations.empty())
            {
                return;
            }
            std::set<vtkActor*> copies;
            for(const auto& entry : images)
            {
                for(const auto& copy : entry.second.copies)
                {
                    copies.insert(copy);
                }
            }
            std::set<vtkActor*> sources;
            auto* props = renderer->GetViewProps();
            vtkCollectionSimpleIterator iterator;
            props->InitTraversal(iterator);
            while(auto* prop = props->GetNextProp(iterator))
            {
                auto* actor = vtkActor::SafeDownCast(prop);
                // Overlays, labels and coordinate axes are not periodic geometry.
                if(actor && actor->GetMapper() && !actor->IsA("vtkCubeAxesActor") && !copies.count(actor))
                {
                    sources.insert(actor);
                }
            }
            for(auto it = images.begin(); it != images.end();)
            {
                if(!sources.count(it->first))
                {
                    for(auto& copy : it->second.copies)
                    {
                        renderer->RemoveActor(copy);
                    }
                    it = images.erase(it);
                }
                else
                {
                    ++it;
                }
            }
            for(auto* source : sources)
            {
                auto& entry = images[source];
                entry.source = source;
                if(entry.copies.empty())
                {
                    for(size_t i = 0; i < translations.size(); ++i)
                    {
                        auto copy = vtkSmartPointer<vtkActor>::New();
                        renderer->AddActor(copy);
                        entry.copies.push_back(copy);
                    }
                }
                for(size_t i = 0; i < translations.size(); ++i)
                {
                    auto& copy = entry.copies[i];
                    copy->ShallowCopy(source);
                    // Compose translation after the complete source transform.
                    auto matrix = vtkSmartPointer<vtkMatrix4x4>::New();
                    matrix->DeepCopy(source->GetMatrix());
                    for(int row = 0; row < 3; ++row)
                    {
                        matrix->SetElement(row, 3, matrix->GetElement(row, 3) + translations[i][row]);
                    }
                    copy->SetOrigin(0, 0, 0);
                    copy->SetPosition(0, 0, 0);
                    copy->SetOrientation(0, 0, 0);
                    copy->SetScale(1, 1, 1);
                    copy->SetUserMatrix(matrix);
                    // Picking must continue to refer to original simulation entities.
                    copy->PickableOff();
                }
            }
        }
    };

    class PeriodicImagesWidget : public QWidget
    {
        const Eigen::Matrix<double, 3, Eigen::Dynamic> basis;
        PeriodicImages images;
        vtkGenericOpenGLRenderWindow* window;
        vtkRenderer* renderer;
        QGroupBox* group;
        std::array<QSpinBox*, 3> counts{{nullptr, nullptr, nullptr}};

        void updateImages()
        {
            std::vector<Eigen::Vector3d> shifts;
            if(group->isChecked())
            {
                std::array<int, 3> n{{0, 0, 0}};
                for(int i = 0; i < basis.cols(); ++i)
                {
                    n[i] = counts[i]->value();
                }
                for(int i = -n[0]; i <= n[0]; ++i)
                {
                    for(int j = -n[1]; j <= n[1]; ++j)
                    {
                        for(int k = -n[2]; k <= n[2]; ++k)
                        {
                            if(i != 0 || j != 0 || k != 0)
                            {
                                Eigen::Vector3d shift = Eigen::Vector3d::Zero();
                                const std::array<int, 3> indices{{i, j, k}};
                                for(int axis = 0; axis < basis.cols(); ++axis)
                                {
                                    shift += indices[axis] * basis.col(axis);
                                }
                                shifts.push_back(shift);
                            }
                        }
                    }
                }
            }
            images.setTranslations(std::move(shifts));
            renderer->ResetCameraClippingRange();
            if(window->GetReadyForRendering())
            {
                window->Render();
            }
        }

      public:
        PeriodicImagesWidget(vtkGenericOpenGLRenderWindow* windowIn, vtkRenderer* rendererIn,
                             const Eigen::Matrix<double, 3, Eigen::Dynamic>& basisIn, QWidget* parent) :
            QWidget(parent), basis(basisIn), images(rendererIn), window(windowIn), renderer(rendererIn),
            group(controls::box("Periodic images", this,
                "Display translated copies of the spatial actors along the active periodicity vectors. Fields and clipping are copied from the original cell, without recomputing physics."))
        {
            auto* layout = new QVBoxLayout(this);
            layout->setContentsMargins(0, 0, 0, 0);
            layout->addWidget(group);
            group->setObjectName("periodicImagesGroup");
            group->setCheckable(true);
            group->setChecked(false);
            group->setEnabled(basis.cols() > 0);
            for(int i = 0; i < basis.cols(); ++i)
            {
                counts[i] = new QSpinBox(group);
                counts[i]->setObjectName(QString("periodicImageCount%1").arg(i + 1));
                counts[i]->setRange(0, 5);
                counts[i]->setKeyboardTracking(false);
                controls::field(controls::grid(group), i,
                    QString("v%1: (%2, %3, %4)").arg(i + 1).arg(basis(0, i), 0, 'g', 5)
                        .arg(basis(1, i), 0, 'g', 5).arg(basis(2, i), 0, 'g', 5),
                    counts[i], "Additional cells on each side along this vector (shown in simulation length units). 0: original only; 1: indices -1, 0, +1. Counts combine across directions. Maximum 5 per side.");
                connect(counts[i], &QSpinBox::valueChanged, this, [this](int)
                {
                    updateImages();
                });
            }
            connect(group, &QGroupBox::toggled, this, [this](bool)
            {
                updateImages();
            });
        }
    };
}
#endif
