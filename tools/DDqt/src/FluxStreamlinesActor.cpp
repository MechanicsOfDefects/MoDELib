#include <FluxStreamlinesActor.h>
#include <FluxStreamlineTrace.h>
#include <DDqtControls.h>
#include <ClusterDynamics.h>
#include <ClusterDynamics/ClusterDynamicsBase.h>
#include <DislocationNetwork.h>
#include <QComboBox>
#include <QCheckBox>
#include <limits>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QDoubleValidator>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include <cmath>
#include <exception>
namespace model
{
    FluxStreamlinesActor::FluxStreamlinesActor(vtkGenericOpenGLRenderWindow* w,vtkRenderer* r,
                                               const DefectiveCrystal<3>& c):
        crystal(c),window(w),renderer(r),data(vtkSmartPointer<vtkPolyData>::New()),
        mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),actor(vtkSmartPointer<vtkActor>::New())
    {
        auto* layout=new QVBoxLayout(this);
        layout->setContentsMargins(0,0,0,0);
        group=controls::box("Flux streamlines",this,"Three-dimensional total mobile-species flux seeded near dislocations. Physical surfaces and grain boundaries terminate traces.");
        group->setObjectName("fluxStreamlinesGroup");
        group->setCheckable(true);
        group->setChecked(false);
        layout->addWidget(group);
        auto* grid=controls::grid(group);
        species=new QComboBox(this);
        species->setObjectName("fluxSpecies");
        if(c.ddBase.clusterDynamicsBase)
        {
            for(const auto value : c.ddBase.clusterDynamics().msVector)
            {
                const int n=static_cast<int>(value);
                species->addItem("J_"+(std::abs(n)==1?QString():QString::number(std::abs(n)))+(n<0?"v":"i"));
            }
        }
        direction=new QComboBox(this);
        direction->addItems({"Forward","Backward","Both"});
        direction->setCurrentIndex(2);
        controls::field(grid,0,"Species",species,"Trace the selected species' flux, including the CD correction and discrete-dislocation concentration.");
        controls::field(grid,1,"Direction",direction,"Forward follows flux; backward shows where absorbed species originate. Both traces each seed in both directions.");
        logColors=new QCheckBox(this);
        logColors->setObjectName("fluxLogColors");
        logColors->setChecked(true);
        controls::field(grid,2,"Logarithmic colors",logColors,"Spread flux magnitudes over a logarithmic color scale. Equal ratios have equal color spacing. Changes recolor existing traces without recomputing them.");
        connect(logColors,&QCheckBox::toggled,this,[this](){updateColors();window->Render();});
        const auto number=[&](int row,const char* label,double value,const char* help)
        {
            auto* edit=new QDoubleSpinBox(this);
            edit->setDecimals(8);
            edit->setRange(1.e-8,1.e12);
            edit->setValue(value);
            controls::field(grid,row,label,edit,help);
            return edit;
        };
        spacing=number(3,"Seed spacing [b]",100,"Spacing along each nonzero-Burgers segment. Four offset seeds surround each sample.");
        offset=number(4,"Seed offset [b]",10,"Radial distance from the dislocation chord. Seeds outside the mesh are skipped. Use an offset larger than the core width.");
        step=number(5,"Step length [b]",10,"Maximum arclength step. Reduce to check trajectory convergence; tetrahedron faces are crossed through mesh connectivity.");
        length=number(6,"Trace length [b]",1000,"Maximum length per seed and tracing direction, including periodic wraps.");
        minimum=new QLineEdit("0",this);
        minimum->setValidator(new QDoubleValidator(0.0,1.e300,16,minimum));
        controls::field(grid,8,"Minimum flux",minimum,"Stop below this flux magnitude in simulation units. Scientific notation such as 1e-25 is accepted; zero stops only at zero flux.");
        width=number(9,"Line width",2,"Rendered streamline width in pixels.");
        width->setMaximum(20);
        seedLimit=new QSpinBox(this);
        seedLimit->setObjectName("fluxSeedLimit");
        seedLimit->setRange(1,10000);
        seedLimit->setValue(64);
        controls::field(grid,10,"Maximum seeds",seedLimit,"Bound the cost of tracing; samples are distributed across the network when this limit is reached.");
        stepLimit=new QSpinBox(this);
        stepLimit->setRange(1,100000);
        stepLimit->setValue(2000);
        controls::field(grid,11,"Maximum steps",stepLimit,"Maximum integration steps per seed and direction; bounds closed or slowly progressing traces.");
        auto* computeButton=new QPushButton("Compute streamlines",this);
        computeButton->setObjectName("computeFluxStreamlines");
        grid->addWidget(computeButton,12,0,1,2);
        status=new QLabel(this);
        status->setWordWrap(true);
        grid->addWidget(status,13,0,1,2);
        mapper->SetInputData(data);
        auto lut=vtkSmartPointer<vtkLookupTable>::New();
        lut->SetHueRange(0.667,0.0);
        lut->Build();
        mapper->SetLookupTable(lut);
        mapper->SetScalarModeToUsePointFieldData();
        mapper->SelectColorArray("Relative flux magnitude");
        mapper->SetColorModeToMapScalars();
        mapper->ScalarVisibilityOn();
        mapper->SetScalarRange(0.0,1.0);
        actor->SetMapper(mapper);
        actor->SetVisibility(false);
        renderer->AddActor(actor);
        connect(computeButton,&QPushButton::clicked,this,[this](){compute();window->Render();});
        const auto invalidatePlot=[this](){invalidate();window->Render();};
        connect(group,&QGroupBox::toggled,this,invalidatePlot);
        connect(species,&QComboBox::currentIndexChanged,this,invalidatePlot);
        connect(direction,&QComboBox::currentIndexChanged,this,invalidatePlot);
        for(auto* control : {spacing,offset,step,length})
        {
            connect(control,&QDoubleSpinBox::valueChanged,this,invalidatePlot);
        }
        connect(minimum,&QLineEdit::textChanged,this,invalidatePlot);
        connect(seedLimit,&QSpinBox::valueChanged,this,invalidatePlot);
        connect(stepLimit,&QSpinBox::valueChanged,this,invalidatePlot);
        status->setText("Choose parameters, then click Compute streamlines.");
        connect(width,&QDoubleSpinBox::valueChanged,this,[this](double v){actor->GetProperty()->SetLineWidth(v);window->Render();});
        if(!c.ddBase.clusterDynamicsBase)
        {
            group->setEnabled(false);
            group->setToolTip("Enable ClusterDynamics to display species flux streamlines.");
        }
    }
    FluxStreamlinesActor::~FluxStreamlinesActor()
    {
        renderer->RemoveActor(actor);
    }
    void FluxStreamlinesActor::updateColors()
    {
        auto* scalars=data->GetPointData()->GetArray("Flux magnitude");
        if(!scalars || !scalars->GetNumberOfTuples())
        {
            return;
        }
        double range[2];
        scalars->GetRange(range);
        double positiveMin=std::numeric_limits<double>::infinity();
        for(vtkIdType i=0;i<scalars->GetNumberOfTuples();++i)
        {
            const double value=scalars->GetTuple1(i);
            if(value>0)
            {
                positiveMin=std::min(positiveMin,value);
            }
        }
        const bool logarithmic=logColors->isChecked() && std::isfinite(positiveMin);
        const long double low=logarithmic?std::log(static_cast<long double>(positiveMin)):range[0];
        const long double high=logarithmic?std::log(static_cast<long double>(range[1])):range[1];
        auto relative=vtkSmartPointer<vtkDoubleArray>::New();
        relative->SetName("Relative flux magnitude");
        for(vtkIdType i=0;i<scalars->GetNumberOfTuples();++i)
        {
            const double physical=scalars->GetTuple1(i);
            const long double value=logarithmic?(physical>0?std::log(static_cast<long double>(physical)):low):physical;
            relative->InsertNextValue(high>low?static_cast<double>((value-low)/(high-low)):0.5);
        }
        data->GetPointData()->AddArray(relative);
        data->Modified();
    }

    void FluxStreamlinesActor::invalidate()
    {
        data->Initialize();
        actor->SetVisibility(false);
        status->setText("Parameters or frame changed. Click Compute streamlines.");
    }
    void FluxStreamlinesActor::compute()
    {
        data->Initialize();
        actor->SetVisibility(false);
        if(!group->isChecked())
        {
            return;
        }
        try
        {
            const auto cd=crystal.getUniqueTypedMicrostructure<ClusterDynamics<3>>();
            const auto dn=crystal.getUniqueTypedMicrostructure<DislocationNetwork<3>>();
            if(!cd || !dn || species->currentIndex()<0)
            {
                status->setText("Requires ClusterDynamics and DislocationDynamics.");
                return;
            }
            const int k=species->currentIndex();
            const auto& mesh=crystal.ddBase.mesh;
            // Local snapshots are rebuilt for each explicit Compute, never across frames.
            using Source=DislocationSegment<3>::MobileClusterFluxSource;
            std::map<size_t,std::vector<Source>> ddSources;
            if(dn->climbSolver)
            {
                for(const auto& entry:dn->networkLinks())
                {
                    const auto segment=entry.second.lock();
                    if(segment && !segment->hasZeroBurgers() && segment->grains().size()==1)
                    {
                        auto source=segment->mobileClusterFluxSource(dn->climbSolver->CD->cdp,k);
                        if(source.qA!=0 || source.qB!=0)
                        {
                            ddSources[(*segment->grains().begin())->grainID].push_back(std::move(source));
                        }
                    }
                }
            }
            const auto& sources=ddSources;

            const FluxEvaluator evaluate=[&](const Eigen::Vector3d& x,const Simplex<3,3>& cell)
            {
                using Element=MicrostructureBase<3>::ElementType;
                const Element* element=crystal.ddBase.fe?&crystal.ddBase.fe->elements().at(cell.xID):nullptr;
                FluxSample sample;
                for(const auto& physics : crystal.microstructures())
                {
                    if(physics.get()!=dn.get())
                    {
                        sample.flux+=physics->mobileClusterFlux(x,nullptr,element,&cell).col(k);
                    }
                    else
                    {
                        Eigen::Vector3d ddFlux=Eigen::Vector3d::Zero();
                        const auto found=sources.find(cell.region->regionID);
                        if(found!=sources.end())
                        {
                            for(const auto& source:found->second)
                            {
                                ddFlux+=source.evaluate(x,crystal.ddBase.periodicShifts);
                            }
                        }
                        sample.flux+=ddFlux;
                    }
                }
                return sample;
            };
            struct Segment { Eigen::Vector3d a,b; size_t count; };
            std::vector<Segment> segments;
            size_t total=0;
            for(const auto& entry : dn->networkLinks())
            {
                const auto owned=entry.second.lock();
                if(!owned)
                {
                    continue;
                }
                const auto& link=*owned;
                const Eigen::Vector3d a=link.source->get_P(),b=link.sink->get_P();
                if(link.burgers().squaredNorm()>0.0 && (b-a).norm()>0.0)
                {
                    const size_t count=static_cast<size_t>(std::min(1.e6,std::max(1.0,std::ceil((b-a).norm()/spacing->value()))));
                    segments.push_back({a,b,count});
                    total+=4*count;
                }
            }
            const size_t used=std::min(total,static_cast<size_t>(seedLimit->value()));
            auto points=vtkSmartPointer<vtkPoints>::New();
            points->SetDataTypeToDouble();
            auto lines=vtkSmartPointer<vtkCellArray>::New();
            auto scalars=vtkSmartPointer<vtkDoubleArray>::New();
            scalars->SetName("Flux magnitude");
            FluxTraceOptions options;
            options.step=step->value();options.length=length->value();
            bool validMinimum=false;
            options.minimumFlux=minimum->text().toDouble(&validMinimum);
            if(!validMinimum || !std::isfinite(options.minimumFlux) || options.minimumFlux<0.0)
            {
                throw std::runtime_error("Minimum flux must be a finite nonnegative number.");
            }
            options.maxSteps=stepLimit->value();
            struct TraceJob
            {
                Eigen::Vector3d seed;
                int direction;
            };
            std::vector<TraceJob> jobs;
            jobs.reserve(2*used);
            size_t segmentID=0,base=0;
            for(size_t seed=0;seed<used;++seed)
            {
                const size_t index=static_cast<size_t>((static_cast<long double>(seed)+0.5L)*total/used);
                while(index>=base+4*segments[segmentID].count)
                {
                    base+=4*segments[segmentID].count;
                    ++segmentID;
                }
                const auto& segment=segments[segmentID];
                const size_t local=index-base;
                const Eigen::Vector3d tangent=(segment.b-segment.a).normalized();
                const Eigen::Vector3d u=tangent.unitOrthogonal(),v=tangent.cross(u);
                const double angle=(seed%4)*1.5707963267948966;
                const Eigen::Vector3d x=segment.a+(double(local/4)+0.5)/segment.count*(segment.b-segment.a)+
                    offset->value()*(std::cos(angle)*u+std::sin(angle)*v);
                for(const int sign : {1,-1})
                {
                    if((direction->currentIndex()==0 && sign<0) || (direction->currentIndex()==1 && sign>0))
                    {
                        continue;
                    }
                    jobs.push_back({x,sign});
                }
            }
            // Freeze controls and geometry on the UI thread. Workers only read
            // the current physics state and write to their own result slots.
            std::vector<FluxTrace> traces(jobs.size());
            std::vector<std::exception_ptr> failures(jobs.size());
            const auto count=static_cast<std::ptrdiff_t>(jobs.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) if(count>1)
#endif
            for(std::ptrdiff_t i=0;i<count;++i)
            {
                try
                {
                    traces[i]=traceFlux(mesh,jobs[i].seed,evaluate,options,jobs[i].direction);
                }
                catch(...)
                {
                    failures[i]=std::current_exception();
                }
            }
            // Propagate failures outside OpenMP before publishing any geometry.
            for(const auto& failure:failures)
            {
                if(failure)
                {
                    std::rethrow_exception(failure);
                }
            }
            // VTK remains on the main thread, in deterministic seed order.
            for(const auto& trace:traces)
            {
                for(const auto& part:trace)
                {
                    if(part.size()<2)
                    {
                        continue;
                    }
                    lines->InsertNextCell(part.size());
                    for(const auto& point:part)
                    {
                        const auto id=points->InsertNextPoint(point.position.data());
                        lines->InsertCellPoint(id);
                        scalars->InsertNextValue(point.magnitude);
                    }
                }
            }
            data->SetPoints(points);data->SetLines(lines);data->GetPointData()->SetScalars(scalars);
            data->Modified();
            double range[2]={0,0};
            if(scalars->GetNumberOfTuples()>0)
            {
                scalars->GetRange(range);

            }
            updateColors();
            actor->GetProperty()->SetLineWidth(width->value());
            actor->SetVisibility(lines->GetNumberOfCells()>0);
            status->setText(QString("%1 / %2 seeds; %3 trace pieces. %4 range: %5 to %6. Click Compute again after changing settings or frames.")
                .arg(used).arg(total).arg(lines->GetNumberOfCells()).arg("Flux magnitude")
                .arg(range[0],0,'g',4).arg(range[1],0,'g',4));
        }
        catch(const std::exception& error)
        {
            data->Initialize();
            status->setText(QString("Streamlines: %1").arg(error.what()));
        }
    }
}
