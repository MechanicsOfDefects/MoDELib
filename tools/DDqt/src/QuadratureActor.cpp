/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_QuadratureActor_cpp_
#define model_QuadratureActor_cpp_

#include <vtkLine.h>
#include <vtkStructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkStructuredGridAppend.h>

#include <QuadratureActor.h>

#include <DDqtControls.h>

namespace model
{

    QuadratureActor::QuadratureActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,const Polycrystal<3>&,const DefectiveCrystalTraitsIO&) :
    /* init */ renderWindow(renWin)
    /* init */,renderer(ren)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,pkBox(new controls::GroupBox(tr("&PK forces")))
    /* init */,pkScaleEdit(new QLineEdit("10000"))
    /* init */,sfBox(new controls::GroupBox(tr("&stacking fault forces")))
    /* init */,sfScaleEdit(new QLineEdit("10000"))
    /* init */,ltBox(new controls::GroupBox(tr("&line tension forces")))
    /* init */,ltScaleEdit(new QLineEdit("10000"))
    /* init */,velocityBox(new controls::GroupBox(tr("&velocities")))
    /* init */,velocityScaleEdit(new QLineEdit("100"))
    /* init */,pkPolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,pkGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,pkMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,pkActor(vtkSmartPointer<vtkActor>::New())
    /* init */,sfPolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,sfGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,sfMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,sfActor(vtkSmartPointer<vtkActor>::New())
    /* init */,ltPolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,ltGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,ltMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,ltActor(vtkSmartPointer<vtkActor>::New())
    /* init */,velocityPolydata(vtkSmartPointer<vtkPolyData>::New())
    /* init */,velocityGlyphs(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,velocityMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,velocityActor(vtkSmartPointer<vtkActor>::New())
{
        
        pkBox->setCheckable(true);
        pkBox->setChecked(false);
        pkActor->SetVisibility(false);
        QGridLayout *pkLayout = new QGridLayout();
        controls::field(
            pkLayout, 0, "Scale", pkScaleEdit,
            "Display multiplier for Peach-Koehler force arrows at quadrature points. Press Enter to apply.");
        pkBox->setToolTip("Show or hide Peach-Koehler force arrows at quadrature points.");
        pkBox->setLayout(pkLayout);
    
        sfBox->setCheckable(true);
        sfBox->setChecked(false);
        sfActor->SetVisibility(false);
        QGridLayout *sfLayout = new QGridLayout();
        controls::field(
            sfLayout, 0, "Scale", sfScaleEdit,
            "Display multiplier for stacking-fault force arrows at quadrature points. Press Enter to apply.");
        sfBox->setToolTip("Show or hide stacking-fault force arrows at quadrature points.");
        sfBox->setLayout(sfLayout);

    ltBox->setCheckable(true);
    ltBox->setChecked(false);
    ltActor->SetVisibility(false);
    QGridLayout *ltLayout = new QGridLayout();
    controls::field(ltLayout, 0, "Scale", ltScaleEdit,
                    "Display multiplier for line-tension force arrows at quadrature points. Press Enter to apply.");
    ltBox->setToolTip("Show or hide line-tension force arrows at quadrature points.");
    ltBox->setLayout(ltLayout);

    
        velocityBox->setCheckable(true);
        velocityBox->setChecked(false);
        velocityActor->SetVisibility(false);
        QGridLayout *velocityLayout = new QGridLayout();
        controls::field(velocityLayout, 0, "Scale", velocityScaleEdit,
                        "Display multiplier for velocity arrows at quadrature points. Press Enter to apply.");
        velocityBox->setToolTip("Show or hide velocity arrows at quadrature points.");
        velocityBox->setLayout(velocityLayout);

        
        mainLayout->addWidget(pkBox,0,0,1,1);
        mainLayout->addWidget(sfBox,1,0,1,1);
        mainLayout->addWidget(ltBox,2,0,1,1);
        mainLayout->addWidget(velocityBox,3,0,1,1);
        mainLayout->setRowStretch(4, 1);

        this->setLayout(mainLayout);

        connect(pkBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
        connect(pkScaleEdit,SIGNAL(returnPressed()), this, SLOT(modify()));
        connect(sfBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
        connect(sfScaleEdit,SIGNAL(returnPressed()), this, SLOT(modify()));
        connect(ltBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
        connect(ltScaleEdit,SIGNAL(returnPressed()), this, SLOT(modify()));
        connect(velocityBox,SIGNAL(toggled(bool)), this, SLOT(modify()));
        connect(velocityScaleEdit,SIGNAL(returnPressed()), this, SLOT(modify()));

        
        pkGlyphs->SetInputData(pkPolydata);
        pkGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
        pkGlyphs->ScalingOn();
        pkGlyphs->SetScaleModeToScaleByVector();
        pkGlyphs->OrientOn();
        pkGlyphs->ClampingOff();
        pkGlyphs->SetVectorModeToUseVector();
        pkGlyphs->SetIndexModeToOff();
        pkMapper->SetInputConnection(pkGlyphs->GetOutputPort());
        pkMapper->ScalarVisibilityOff();
        pkActor->SetMapper(pkMapper);
        pkActor->GetProperty()->SetColor(0.0, 0.0, 1.0); //(R,G,B)
        renderer->AddActor(pkActor);

        sfGlyphs->SetInputData(sfPolydata);
        sfGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
        sfGlyphs->ScalingOn();
        sfGlyphs->SetScaleModeToScaleByVector();
        sfGlyphs->OrientOn();
        sfGlyphs->ClampingOff();
        sfGlyphs->SetVectorModeToUseVector();
        sfGlyphs->SetIndexModeToOff();
        sfMapper->SetInputConnection(sfGlyphs->GetOutputPort());
        sfMapper->ScalarVisibilityOff();
        sfActor->SetMapper(sfMapper);
        sfActor->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
        renderer->AddActor(sfActor);

    ltGlyphs->SetInputData(ltPolydata);
    ltGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
    ltGlyphs->ScalingOn();
    ltGlyphs->SetScaleModeToScaleByVector();
    ltGlyphs->OrientOn();
    ltGlyphs->ClampingOff();
    ltGlyphs->SetVectorModeToUseVector();
    ltGlyphs->SetIndexModeToOff();
    ltMapper->SetInputConnection(ltGlyphs->GetOutputPort());
    ltMapper->ScalarVisibilityOff();
    ltActor->SetMapper(ltMapper);
    ltActor->GetProperty()->SetColor(1.0, 0.0, 0.0); //(R,G,B)
    renderer->AddActor(ltActor);

    
    
        velocityGlyphs->SetInputData(velocityPolydata);
        velocityGlyphs->SetSourceConnection(vtkSmartPointer<vtkArrowSource>::New()->GetOutputPort());
        velocityGlyphs->ScalingOn();
        velocityGlyphs->SetScaleModeToScaleByVector();
        velocityGlyphs->OrientOn();
        velocityGlyphs->ClampingOff();
        velocityGlyphs->SetVectorModeToUseVector();
        velocityGlyphs->SetIndexModeToOff();
        velocityMapper->SetInputConnection(velocityGlyphs->GetOutputPort());
        velocityMapper->ScalarVisibilityOff();
        velocityActor->SetMapper(velocityMapper);
        velocityActor->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
        renderer->AddActor(velocityActor);

    }

    void QuadratureActor::updateConfiguration(const DefectiveCrystalAuxIO<3>& auxIO)
    {
        
        std::cout<<"Updating QuadraturePoints..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        
        vtkSmartPointer<vtkPoints> nodePoints = vtkSmartPointer<vtkPoints>::New();
        
        vtkSmartPointer<vtkDoubleArray> pkForces(vtkSmartPointer<vtkDoubleArray>::New());
        pkForces->SetNumberOfComponents(3);

        vtkSmartPointer<vtkDoubleArray> sfForces(vtkSmartPointer<vtkDoubleArray>::New());
        sfForces->SetNumberOfComponents(3);

        vtkSmartPointer<vtkDoubleArray> ltForces(vtkSmartPointer<vtkDoubleArray>::New());
        ltForces->SetNumberOfComponents(3);

        
        vtkSmartPointer<vtkDoubleArray> velocities(vtkSmartPointer<vtkDoubleArray>::New());
        velocities->SetNumberOfComponents(3);
        
        for(const auto& qp : auxIO.quadraturePoints())
        {
            nodePoints->InsertNextPoint(qp.r.data());
            pkForces->InsertNextTuple(qp.pkForce.data()); // arrow vector
            sfForces->InsertNextTuple(qp.stackingFaultForce.data()); // arrow vector
            ltForces->InsertNextTuple(qp.lineTensionForce.data()); // arrow vector
            velocities->InsertNextTuple(qp.velocity.data()); // arrow vector
        }
        
        pkPolydata->SetPoints(nodePoints);
        pkPolydata->GetPointData()->SetVectors(pkForces);
        pkPolydata->Modified();
        
        sfPolydata->SetPoints(nodePoints);
        sfPolydata->GetPointData()->SetVectors(sfForces);
        sfPolydata->Modified();

        ltPolydata->SetPoints(nodePoints);
        ltPolydata->GetPointData()->SetVectors(ltForces);
        ltPolydata->Modified();

        
        velocityPolydata->SetPoints(nodePoints);
        velocityPolydata->GetPointData()->SetVectors(velocities);
        velocityPolydata->Modified();

        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }

    void QuadratureActor::modify()
    {
        pkActor->SetVisibility(pkBox->isChecked());
        const double pkScaling(std::atof(pkScaleEdit->text() .toStdString().c_str()));
        pkGlyphs->SetScaleFactor(pkScaling);

        sfActor->SetVisibility(sfBox->isChecked());
        const double sfScaling(std::atof(sfScaleEdit->text() .toStdString().c_str()));
        sfGlyphs->SetScaleFactor(sfScaling);

        ltActor->SetVisibility(ltBox->isChecked());
        const double ltScaling(std::atof(ltScaleEdit->text() .toStdString().c_str()));
        ltGlyphs->SetScaleFactor(ltScaling);

        velocityActor->SetVisibility(velocityBox->isChecked());
        const double velocityScaling(std::atof(velocityScaleEdit->text() .toStdString().c_str()));
        velocityGlyphs->SetScaleFactor(velocityScaling);

        
        renderWindow->Render();
    }


} // namespace model
#endif
