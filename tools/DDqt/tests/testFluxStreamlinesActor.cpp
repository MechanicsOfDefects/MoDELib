#include "ClusterPeriodicFixture.h"
#include <FluxStreamlinesActor.h>
#include <MicrostructureGenerator.h>
#include <QApplication>
#include <QGroupBox>
#include <QComboBox>
#include <QSpinBox>
#include <QPushButton>
#include <QCheckBox>
#include <vtkUnsignedCharArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActorCollection.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkScalarsToColors.h>
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc,char** argv)
{
    QApplication app(argc,argv);
    try
    {
        require(argc==2,"Expected repository root");
        const fs::path root(argv[1]);
        const auto folder=fixture(root,0,0,false);
        auto text=read(folder/"inputFiles/DefectiveCrystal.txt");
        parameter(text,"physics","ClusterDynamics DislocationDynamics");
        std::ofstream(folder/"inputFiles/DefectiveCrystal.txt")<<text;
        auto dd=read(root/"Library/DislocationDynamics/DD.txt");
        parameter(dd,"glideSolverType","none");
        std::ofstream(folder/"inputFiles/DD.txt")<<dd;
        {
            DefectiveCrystalBase<3> base(folder.string());
            MicrostructureGenerator generator(base);
            const Eigen::Vector3d normal=Eigen::Vector3d::UnitZ();
            const auto grain=base.poly.grain(1);
            const auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.getFromKey(
                GlidePlaneKey<3>(Eigen::Vector3d::Zero(),grain->reciprocalLatticeDirection(normal)));
            std::vector<Eigen::Vector3d> polygon;
            for(int i=0;i<8;++i)
            {
                const double angle=2*M_PI*i/8;
                polygon.emplace_back(6*std::cos(angle),6*std::sin(angle),0);
            }
            require(generator.insertJunctionLoop(polygon,plane,std::sqrt(8.0/3.0)*normal,normal,Eigen::Vector3d::Zero(),1,DislocationLoopIO<3>::SESSILELOOP),"Loop insertion failed");
            DefectiveCrystal<3> crystal(base);
            crystal.initializeConfiguration(generator.config());
            auto cd=crystal.getUniqueTypedMicrostructure<ClusterDynamics<3>>();
            auto& field=cd->clusterDynamicsFEM->mobileClusters;
            auto values=field->dofVector();
            for(const auto& node : base.fe->nodes())
            {
                values(node.gID*ClusterDynamics<3>::mSize)=1.e-6+1.e-9*node.P0.x();
            }
            field->setFullDofVector(values);
            const auto& element=base.fe->elements().begin()->second;
            const auto& meshNode=*base.fe->nodes().begin();
            const Eigen::Vector3d center=element.simplex.center();
            const Eigen::Vector3d expectedFlux=-cd->cdp.D.at(1)[0]*Eigen::Vector3d(1.e-9,0,0)/cd->cdp.omega;
            const MicrostructureBase<3>& physics=*cd;
            for(const Eigen::Vector3d& actual:std::vector<Eigen::Vector3d>{
                physics.mobileClusterFlux(center,nullptr,&element,&element.simplex).col(0),
                physics.mobileClusterFlux(center,nullptr,nullptr,&element.simplex).col(0),
                physics.mobileClusterFlux(meshNode.P0,&meshNode,nullptr,nullptr).col(0)})
            {
                require((actual-expectedFlux).norm()<1.e-8*expectedFlux.norm(),"Physics flux point information changed the FEM result");
            }
            auto window=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
            window->SetReadyForRendering(false);
            auto renderer=vtkSmartPointer<vtkRenderer>::New();
            window->AddRenderer(renderer);
            FluxStreamlinesActor widget(window,renderer,crystal);
            auto* group=widget.findChild<QGroupBox*>("fluxStreamlinesGroup");
            auto* species=widget.findChild<QComboBox*>("fluxSpecies");
            auto* limit=widget.findChild<QSpinBox*>("fluxSeedLimit");
            require(group && species && limit,"Missing flux controls");
            require(!group->isChecked(),"Flux traces should start disabled");
            require(species->count()==ClusterDynamics<3>::mSize && species->itemText(0)=="J_v","Incorrect species labels");
            limit->setValue(4);
            group->setChecked(true);
            auto* actors=renderer->GetActors();
            actors->InitTraversal();
            auto* actor=actors->GetNextActor();
            require(actor,"Missing streamline actor");
            auto* mapper=vtkPolyDataMapper::SafeDownCast(actor->GetMapper());
            require(mapper,"Missing streamline mapper");
            auto* data=mapper->GetInput();
            require(data,"Missing streamline polydata");
            require(!actor->GetVisibility() && data->GetNumberOfLines()==0,"Enabling flux unexpectedly computed streamlines");
            limit->setValue(5);
            require(data->GetNumberOfLines()==0,"Parameter change unexpectedly computed streamlines");
            auto* compute=widget.findChild<QPushButton*>("computeFluxStreamlines");
            require(compute,"Missing Compute streamlines button");
            compute->click();
            require(actor->GetVisibility() && data->GetNumberOfLines()>0,"Compute button produced no streamlines");
            const double expected=(cd->cdp.D.at(1)[0]*Eigen::Vector3d(1.e-9,0,0)/cd->cdp.omega).norm();
            auto* scalars=data->GetPointData()->GetScalars();
            for(vtkIdType i=0;i<scalars->GetNumberOfTuples();++i)
            {
                require(std::abs(scalars->GetTuple1(i)-expected)<1.e-8*expected,"FEM flux magnitude incorrect");
            }
            // Recomputing a frame with zero concentration removes stale lines.
            values.setZero();
            field->setFullDofVector(values);
            widget.compute();
            require(!actor->GetVisibility() && data->GetNumberOfLines()==0,"Stale streamlines survived zero flux frame");
            // The common API must expose DD's analytical flux and aggregate it.
            const auto dn=crystal.getUniqueTypedMicrostructure<DislocationNetwork<3>>();
            for(const auto& entry:dn->networkNodes())
            {
                entry.second.lock()->climbVelocityScalar.setConstant(0.2);
            }
            const Eigen::Vector3d probe(0,0,5);
            DislocationSegment<3>::ClusterFluxType segmentSum=DislocationSegment<3>::ClusterFluxType::Zero();
            for(const auto& entry:dn->networkLinks())
            {
                segmentSum+=entry.second.lock()->mobileClusterFlux(probe,cd->cdp);
            }
            require(segmentSum.norm()>0,"DD flux fixture has no source");
            const MicrostructureBase<3>& ddPhysics=*dn;
            require((ddPhysics.mobileClusterFlux(probe,nullptr,nullptr,nullptr)-segmentSum).norm()<1.e-12*segmentSum.norm(),"DD flux dispatch differs from segment sum");
            require((crystal.mobileClusterFlux(probe,nullptr,nullptr,nullptr)-segmentSum).norm()<1.e-12*segmentSum.norm(),"Container flux does not aggregate physics");
#ifdef _OPENMP
            omp_set_dynamic(0);
            omp_set_num_threads(1);
#endif
            widget.compute();
            require(actor->GetVisibility() && data->GetNumberOfLines()>0,"Analytical DD flux produced no streamlines");
            auto serial=vtkSmartPointer<vtkPolyData>::New();
            serial->DeepCopy(data);
#ifdef _OPENMP
            omp_set_num_threads(4);
#endif
            widget.compute();
            require(data->GetNumberOfPoints()==serial->GetNumberOfPoints() && data->GetNumberOfLines()==serial->GetNumberOfLines(),"Parallel tracing changed geometry sizes");
            for(vtkIdType i=0;i<data->GetNumberOfPoints();++i)
            {
                double a[3],b[3];
                data->GetPoint(i,a);
                serial->GetPoint(i,b);
                for(int axis=0;axis<3;++axis)
                {
                    require(a[axis]==b[axis],"Parallel tracing changed point positions");
                }
                require(data->GetPointData()->GetScalars()->GetTuple1(i)==serial->GetPointData()->GetScalars()->GetTuple1(i),"Parallel tracing changed flux scalars");
            }
            auto* actualCells=data->GetLines();
            auto* expectedCells=serial->GetLines();
            actualCells->InitTraversal();
            expectedCells->InitTraversal();
            vtkIdType na,nb;
            const vtkIdType* a;
            const vtkIdType* b;
            while(actualCells->GetNextCell(na,a))
            {
                require(expectedCells->GetNextCell(nb,b) && na==nb,"Parallel tracing changed line connectivity");
                for(vtkIdType j=0;j<na;++j)
                {
                    require(a[j]==b[j],"Parallel tracing changed line point indices");
                }
            }
            require(!widget.findChild<QComboBox*>("fluxColor"),"Obsolete concentration coloring remains");
            for(const auto& entry:dn->networkNodes())
            {
                entry.second.lock()->climbVelocityScalar.setConstant(1.e-30);
            }
            widget.compute();
            require(data->GetNumberOfPoints()>0,"Tiny nonzero flux lost streamlines");
            auto* relative=data->GetPointData()->GetArray("Relative flux magnitude");
            require(relative,"Missing normalized coloring values");
            double normalizedRange[2];
            relative->GetRange(normalizedRange);
            require(normalizedRange[0]==0 && normalizedRange[1]==1,"Tiny flux failed to span color range");
            require(std::string(mapper->GetArrayName())=="Relative flux magnitude","Mapper is not using normalized flux");
            double low[3],high[3];
            mapper->GetLookupTable()->GetColor(0,low);
            mapper->GetLookupTable()->GetColor(1,high);
            require(low[2]>low[0] && high[0]>high[2],"Flux endpoints do not map blue to red");
            auto* logarithmic=widget.findChild<QCheckBox*>("fluxLogColors");
            require(logarithmic && logarithmic->isChecked(),"Missing default logarithmic coloring");
            auto* physical=data->GetPointData()->GetScalars();
            // Strongly skewed values: linear coloring makes nearly all points blue.
            for(vtkIdType i=0;i<physical->GetNumberOfTuples();++i)
            {
                physical->SetTuple1(i,1.e-30*std::pow(10.0,12.0*i/(physical->GetNumberOfTuples()-1)));
            }
            physical->Modified();
            const auto pointCount=data->GetNumberOfPoints();
            logarithmic->setChecked(false);
            logarithmic->setChecked(true);
            require(data->GetNumberOfPoints()==pointCount,"Recoloring retraced geometry");
            auto* rgba=mapper->MapScalars(1.0);
            require(rgba && rgba->GetNumberOfTuples()==pointCount,"Mapper produced no vertex colors");
            double first[4],middle[4],last[4];
            rgba->GetTuple(0,first);
            rgba->GetTuple(pointCount/2,middle);
            rgba->GetTuple(pointCount-1,last);
            require(first[2]>first[0] && last[0]>last[2] && middle[1]>middle[0] && middle[1]>middle[2],"Actual mapped colors do not span blue-green-red for skewed flux");
            for(const auto& entry:dn->networkNodes())
            {
                entry.second.lock()->climbVelocityScalar.setZero();
            }
            widget.compute();
            require(!actor->GetVisibility(),"Compute reused stale cached DD source rates");
            widget.invalidate();
            require(!actor->GetVisibility() && data->GetNumberOfLines()==0,"Frame change retained stale streamlines");
            group->setChecked(false);
            require(!actor->GetVisibility(),"Disabled streamlines remained visible");
        }
        fs::remove_all(folder);
        std::cout<<"PASS: flux controls, analytic FEM flux, dislocation seeds and frame refresh\n";
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
