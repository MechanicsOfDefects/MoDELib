/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ChartActor_cpp_
#define model_ChartActor_cpp_

#include <fstream>
#include <string>
#include <vector>

#include <vtkTable.h>
#include <vtkFloatArray.h>
#include <vtkPlot.h>
#include <vtkPoints.h>
#include <vtkAxis.h>

#include <ChartActor.h>
#include <vtkTextProperty.h>

#include <DDqtControls.h>
#include <cmath>

namespace model
{

        ChartActor::ChartActor(const DefectiveCrystalTraitsIO& traitsIO_in,vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_in,vtkRenderer* const renderer_in) :
        /* init */ traitsIO(traitsIO_in)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,showChart(new QCheckBox(this))
        /* init */,chartUpdate(new QPushButton("Reload",this))
        /* init */,xComboBox(new QComboBox(this))
        /* init */,yComboBox(new QComboBox(this))
        /* init */,xAxisLog(new QCheckBox(this))
        /* init */,yAxisLog(new QCheckBox(this))
        /* init */,customScaleBoxX(new controls::GroupBox(tr("&Custom scale X")))
        /* init */,minScaleX(new QLineEdit(tr("0")))
        /* init */,maxScaleX(new QLineEdit(tr("1")))
        /* init */,customScaleBoxY(new controls::GroupBox(tr("&Custom scale Y")))
        /* init */,minScaleY(new QLineEdit(tr("0")))
        /* init */,maxScaleY(new QLineEdit(tr("1")))
        /* init */,renderWindow(renderWindow_in)
        /* init */,renderer(renderer_in)
        /* init */,chart(vtkSmartPointer<vtkChartXY>::New())
        /* init */,chartScene(vtkSmartPointer<vtkContextScene>::New())
        /* init */,chartActor(vtkSmartPointer<vtkContextActor>::New())
        /* init */,table(vtkSmartPointer<vtkTable>::New())
        /* init */,currentTable(vtkSmartPointer<vtkTable>::New())
        /* init */,colors(vtkSmartPointer<vtkNamedColors>::New())
        /* init */,points(chart->AddPlot(vtkChart::LINE))
        /* init */,currentPoints(chart->AddPlot(vtkChart::LINE))
        {

            showChart->setChecked(false);
            showChart->setText("chart");
            chartUpdate->setEnabled(false);
            xComboBox->setEnabled(false);
            yComboBox->setEnabled(false);
            chartActor->SetVisibility(false);
            
            xAxisLog->setChecked(false);
            xAxisLog->setEnabled(false);
            xAxisLog->setText("log-scale");

            yAxisLog->setChecked(false);
            yAxisLog->setEnabled(false);
            yAxisLog->setText("log-scale");
            
            vtkColor3d color3d1 = colors->GetColor3d("black");
            points->SetColorF(color3d1.GetRed(), color3d1.GetGreen(), color3d1.GetBlue());
            points->SetWidth(5.0);

            vtkColor3d color3d2 = colors->GetColor3d("magenta");
            currentPoints->SetColorF(color3d2.GetRed(), color3d2.GetGreen(), color3d2.GetBlue());
            currentPoints->SetWidth(5.0);

            std::ifstream flabFile(traitsIO.flabFile);
            if(flabFile)
            {
                std::cout<<"Reading "<<traitsIO.flabFile<<" ... "<<std::flush;
                const auto t0= std::chrono::system_clock::now();

                std::string line;
                while (std::getline(flabFile, line))
                {
                    if(!line.empty())
                    {
//                        std::cout<<line<<std::endl;
                        xComboBox->addItem(QString::fromStdString(line));
                        yComboBox->addItem(QString::fromStdString(line));
                        
//                        vtkNew<vtkFloatArray> farr;
//                        farr->SetName(line.c_str());
//                        table->AddColumn(farr);
                    }
                }
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            else
            {
                std::cout<<"Unable to read "<<traitsIO.flabFile<<std::endl;
            }
            
            QGridLayout* autoscaleLayoutX = new QGridLayout();
            customScaleBoxX->setToolTip("Override automatic X-axis bounds with the entered limits.");
            customScaleBoxX->setCheckable(true);
//            customScaleBoxX->setEnabled(false);
            customScaleBoxX->setChecked(false);
            controls::field(autoscaleLayoutX, 0, "Minimum", minScaleX,
                            "X-axis bound used when custom scaling is enabled. Press Enter to apply.");

            controls::field(autoscaleLayoutX, 1, "Maximum", maxScaleX,
                            "X-axis bound used when custom scaling is enabled. Press Enter to apply.");

            customScaleBoxX->setLayout(autoscaleLayoutX);

            QGridLayout* autoscaleLayoutY = new QGridLayout();
            customScaleBoxY->setToolTip("Override automatic Y-axis bounds with the entered limits.");
            customScaleBoxY->setCheckable(true);
//            customScaleBoxY->setEnabled(false);
            customScaleBoxY->setChecked(false);
            controls::field(autoscaleLayoutY, 0, "Minimum", minScaleY,
                            "Y-axis bound used when custom scaling is enabled. Press Enter to apply.");

            controls::field(autoscaleLayoutY, 1, "Maximum", maxScaleY,
                            "Y-axis bound used when custom scaling is enabled. Press Enter to apply.");

            customScaleBoxY->setLayout(autoscaleLayoutY);

            auto *chartBox = controls::box("Chart", this, "Show or hide the simulation-history chart.", showChart);
            auto *chartLayout = controls::grid(chartBox);
            chartLayout->addWidget(chartUpdate, 0, 0, 1, 2);
            chartUpdate->setToolTip("Reload the simulation history file and refresh the chart.");
            auto *xBox =
                controls::box("X axis", this, "Select the plotted quantity, axis scaling and optional bounds.");
            controls::field(controls::grid(xBox), 0, "Quantity", xComboBox,
                            "Simulation-history column to plot on this axis.");
            controls::grid(xBox)->addWidget(xAxisLog, 1, 0, 1, 2);
            xAxisLog->setToolTip("Use a logarithmic axis. Displayed values must be positive.");
            controls::grid(xBox)->addWidget(customScaleBoxX, 2, 0, 1, 2);
            chartLayout->addWidget(xBox, 1, 0, 1, 2);
            auto *yBox =
                controls::box("Y axis", this, "Select the plotted quantity, axis scaling and optional bounds.");
            controls::field(controls::grid(yBox), 0, "Quantity", yComboBox,
                            "Simulation-history column to plot on this axis.");
            controls::grid(yBox)->addWidget(yAxisLog, 1, 0, 1, 2);
            yAxisLog->setToolTip("Use a logarithmic axis. Displayed values must be positive.");
            controls::grid(yBox)->addWidget(customScaleBoxY, 2, 0, 1, 2);
            chartLayout->addWidget(yBox, 2, 0, 1, 2);
            mainLayout->addWidget(chartBox, 0, 0);
            mainLayout->setRowStretch(1, 1);
            this->setLayout(mainLayout);

            connect(showChart,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(xAxisLog, SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(yAxisLog, SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(xComboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
            connect(yComboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
            connect(chartUpdate,SIGNAL(released()), this, SLOT(reload()));
            connect(customScaleBoxX, &QGroupBox::toggled, this, &ChartActor::modify);
            connect(customScaleBoxY, &QGroupBox::toggled, this, &ChartActor::modify);
            for(auto *edit : {minScaleX, maxScaleX, minScaleY, maxScaleY})
            {
                connect(edit, &QLineEdit::returnPressed, this, &ChartActor::modify);
            }

            chartScene->AddItem(chart);
            chartActor->SetScene(chartScene);
            renderer->AddActor(chartActor);
            chartScene->SetRenderer(renderer);
            
        }
        
    void ChartActor::reload()
    {
        table->Initialize();

        const int cols(xComboBox->count());

        for(int c=0;c<cols;++c)
        {
            vtkNew<vtkFloatArray> farr;
            farr->SetName(xComboBox->itemText(c).toStdString().c_str());
            table->AddColumn(farr);
        }
        
        
        std::ifstream fFile(traitsIO.fFile);
        if(fFile)
        {
            std::cout<<"Reading "<<traitsIO.fFile<<" ... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            double temp(0.0);
            std::string line;
            while (std::getline(fFile, line))
            {
                if(!line.empty())
                {
                    const auto row=table->InsertNextBlankRow(cols);
                    std::stringstream ss(line);
                    for(int col=0;col<cols;++col)
                    {
                        ss >> temp;
                        table->SetValue(row, col, temp);
                    }
                }
            }
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        else
        {
            std::cout<<"Unable to read "<<traitsIO.fFile<<std::endl;
        }

        table->Modified();
        updatePlots();
    }

    void ChartActor::updateConfiguration(const size_t& frameID)
    {
        
        if(frameID==0)
        {
            reload();
        }
        
        std::cout<<"Updating chart..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();

        
        currentTable->Initialize();
        
        for(long long k=0;k<table->GetNumberOfColumns();++k)
        {
            vtkNew<vtkFloatArray> farr;
            farr->SetName(table->GetColumnName(k));
            currentTable->AddColumn(farr);
        }
        
        for(long long k=0;k<table->GetNumberOfRows();++k)
        {
            if(table->GetValue(k,0)<=frameID)
            {
                const auto row=currentTable->InsertNextBlankRow(table->GetNumberOfColumns());
                for(int col=0;col<table->GetNumberOfColumns();++col)
                {
                    currentTable->SetValue(row, col, table->GetValue(row,col));
                }
            }
        }
        currentTable->Modified();
        updatePlots();
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;

    }

void ChartActor::updatePlots()
{
    if(xComboBox->count() && yComboBox->count())
    {
        if(   table->GetNumberOfColumns()>xComboBox->currentIndex()
           && table->GetNumberOfColumns()>yComboBox->currentIndex())
        {
            points->SetInputData(table, xComboBox->currentIndex(), yComboBox->currentIndex());
        }
        
        if(   currentTable->GetNumberOfColumns()>xComboBox->currentIndex()
           && currentTable->GetNumberOfColumns()>yComboBox->currentIndex())
        {
            currentPoints->SetInputData(currentTable, xComboBox->currentIndex(), yComboBox->currentIndex());
        }

        const auto applyRange = [](vtkAxis *axis, QGroupBox *group, QLineEdit *minimum, QLineEdit *maximum)
        {
            bool minimumOK = false;
            bool maximumOK = false;
            const double lower = minimum->text().toDouble(&minimumOK);
            const double upper = maximum->text().toDouble(&maximumOK);
            const bool valid = minimumOK && maximumOK && std::isfinite(lower) && std::isfinite(upper) &&
                               lower < upper && (!axis->GetLogScale() || lower > 0.0);
            const bool invalid = group->isChecked() && !valid;
            for(auto *edit : {minimum, maximum})
            {
                edit->setStyleSheet(invalid ? "QLineEdit { border: 1px solid #c44; }" : "");
                edit->setToolTip("Enter finite limits with minimum < maximum; logarithmic limits must be positive. "
                                 "Press Enter to apply. Invalid limits use automatic scaling.");
            }
            axis->SetBehavior(group->isChecked() && valid ? vtkAxis::FIXED : vtkAxis::AUTO);
            if(group->isChecked() && valid)
            {
                axis->SetUnscaledRange(lower, upper);
            }
        };
        applyRange(chart->GetAxis(1), customScaleBoxX, minScaleX, maxScaleX);
        applyRange(chart->GetAxis(0), customScaleBoxY, minScaleY, maxScaleY);
        chart->RecalculateBounds();
        chart->Modified();

    }
}

 
        
        void ChartActor::modify()
        {
//            std::cout<<"modifying chart..."<<std::endl;

            xComboBox->setEnabled(showChart->isChecked());
            yComboBox->setEnabled(showChart->isChecked());
            customScaleBoxX->setEnabled(showChart->isChecked());
            customScaleBoxY->setEnabled(showChart->isChecked());

            xAxisLog->setEnabled(showChart->isChecked());
            yAxisLog->setEnabled(showChart->isChecked());
            chartUpdate->setEnabled(showChart->isChecked());


            chart->GetAxis(0)->SetTitle(yComboBox->currentText().toStdString());
            chart->GetAxis(1)->SetTitle(xComboBox->currentText().toStdString());

            chart->GetAxis(0)->SetLogScale(yAxisLog->isChecked());
            chart->GetAxis(1)->SetLogScale(xAxisLog->isChecked());

            chart->GetAxis(0)->SetGridVisible(false);
            chart->GetAxis(1)->SetGridVisible(false);

            chart->GetAxis(0)->SetNumberOfTicks(10);
            chart->GetAxis(1)->SetNumberOfTicks(10);

            chart->GetAxis(0)->GetLabelProperties()->SetFontSize(16);
            chart->GetAxis(1)->GetLabelProperties()->SetFontSize(16);

            chart->GetAxis(0)->GetTitleProperties()->SetFontSize(16);
            chart->GetAxis(1)->GetTitleProperties()->SetFontSize(16);
            
            
            updatePlots();

            
//            std::cout<<"xComboBox->currentIndex()="<<xComboBox->currentIndex()<<std::endl;
//            std::cout<<"yComboBox->currentIndex()="<<yComboBox->currentIndex()<<std::endl;
//            std::cout<<"table->GetNumberOfRows()="<<table->GetNumberOfRows()<<std::endl;
//            std::cout<<"table->GetNumberOfColumns()="<<table->GetNumberOfColumns()<<std::endl;
//            std::cout<<"currentTable->GetNumberOfRows()="<<currentTable->GetNumberOfRows()<<std::endl;
//            std::cout<<"currentTable->GetNumberOfColumns()="<<currentTable->GetNumberOfColumns()<<std::endl;

            
//            chart->ClearPlots();

//            if(xComboBox->count() && yComboBox->count())
//            {
//                if(   table->GetNumberOfColumns()>xComboBox->currentIndex()
//                   && table->GetNumberOfColumns()>yComboBox->currentIndex()
////                   && table->GetNumberOfRows()>2
//                   )
//                {
////                    vtkPlot* points = chart->AddPlot(vtkChart::LINE);
//                    points->SetInputData(table, xComboBox->currentIndex(), yComboBox->currentIndex());
////                    vtkColor3d color3d1 = colors->GetColor3d("black");
////                    points->SetColor(color3d1.GetRed(), color3d1.GetGreen(), color3d1.GetBlue());
////                    points->SetWidth(2.0);
//                }
//
//                if(   currentTable->GetNumberOfColumns()>xComboBox->currentIndex()
//                   && currentTable->GetNumberOfColumns()>yComboBox->currentIndex()
////                   && currentTable->GetNumberOfRows()>2
//                   )
//                {
////                    vtkPlot* currentPoints = chart->AddPlot(vtkChart::LINE);
//                    currentPoints->SetInputData(currentTable, xComboBox->currentIndex(), yComboBox->currentIndex());
////                    vtkColor3d color3d2 = colors->GetColor3d("magenta");
////                    currentPoints->SetColor(color3d2.GetRed(), color3d2.GetGreen(), color3d2.GetBlue());
////                    currentPoints->SetWidth(3.0);
//                }
//
//                chart->RecalculateBounds();
//                chart->Modified();
//
//            }
            
            chartActor->SetVisibility(showChart->isChecked());
            renderWindow->Render(); // this makes code crash if called before the app window is created
            
//            std::cout<<"Min/Max="<<chart->GetAxis(0)->GetMinimum()<<", "<<chart->GetAxis(0)->GetMaximum()<<std::endl;

        }



} // namespace model
#endif







