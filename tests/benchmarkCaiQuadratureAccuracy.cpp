#include <DislocationDynamics/CaiStressParticleKernel.h>
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <DefectiveCrystal/DefectiveCrystal.h>
#include <DislocationDynamics/DislocationNetwork.h>
#include <DislocationDynamics/DislocationFieldBase.h>
#include <DislocationDynamics/StraightDislocationSegment.h>
#include <chrono>
#include <algorithm>
#include <iomanip>
using Clock=std::chrono::steady_clock;
int main(int argc,char** argv)
{
    try
    {
        if(argc<3 || argc>4) { throw std::invalid_argument("Usage: benchmarkCaiQuadratureAccuracy tutorial frame [targets=1024]"); }
        const std::size_t requested=argc>3?std::stoul(argv[3]):1024;
        if(!requested) { throw std::invalid_argument("Positive target count required"); }
        Kokkos::ScopeGuard guard;
        model::DefectiveCrystalBase<3> base(argv[1]);
        model::DefectiveCrystalConfigIO<3> config(base.simulationParameters.traitsIO.configFolder);
        config.read(std::stol(argv[2]));
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        network.initializeConfiguration(config,std::ofstream(),std::ofstream());
        network.updateGeometry();
        const double density=model::TextFileParser(base.simulationParameters.traitsIO.inputFilesFolder+"/DD.txt").readScalar<double>("quadPerLength",true);
        network.fieldEvaluator.rebuild(density,density,model::DislocationTargetPurpose::glide);
        const auto& buffers=network.fieldEvaluator.buffers();
        model::DislocationTargetArrays<3> targets;
        const auto count=std::min(requested,buffers.targets().size());
        if(!count || !buffers.sources().size()) { throw std::runtime_error("Benchmark requires sources and targets"); }
        for(std::size_t t=0;t<count;++t)
        {
            const auto p=t*buffers.targets().size()/count;
            targets.positions.append({buffers.targets().positions.values[0][p],buffers.targets().positions.values[1][p],buffers.targets().positions.values[2][p]});
        }
        model::KokkosDislocationFieldBuffers<3> device;
        device.uploadSources(buffers.sources()); device.uploadTargets(targets);
        const model::CaiStressParticleKernel kernel(base.poly.mu,base.poly.nu,model::DislocationFieldBase<3>::a);
        std::vector<std::shared_ptr<model::DislocationSegment<3>>> segments;
        std::vector<model::StraightDislocationSegment<3>> referenceSegments;
        double maxLength=0;
        for(const auto& entry:network.networkLinks())
        {
            auto segment=entry.second.lock();
            if(!segment->hasZeroBurgers())
            {
                segments.push_back(segment);
                referenceSegments.emplace_back(base.poly,segment->source->get_P(),segment->sink->get_P(),segment->burgers());
                maxLength=std::max(maxLength,segment->chord().norm());
            }
        }
        std::vector<Eigen::Matrix3d> reference(count,Eigen::Matrix3d::Zero());
        // Same unique physical segments, no images/screen/external/FEM stress.
        const auto analytical=[&]
        {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for(std::size_t t=0;t<count;++t)
            {
                const Eigen::Vector3d x(targets.positions.values[0][t],targets.positions.values[1][t],targets.positions.values[2][t]);
                Eigen::Matrix3d sum=Eigen::Matrix3d::Zero();
                for(const auto& segment:referenceSegments) { sum+=segment.stress(x); }
                reference[t]=sum;
            }
        };
        analytical();
        double ref2=0;
        for(const auto& value:reference) { ref2+=value.squaredNorm(); }
        for(double factor:{1.0,2.0,4.0,8.0,16.0})
        {
            if(factor*density*maxLength>1024) { std::cout<<"CAI_ACCURACY skipped_density="<<factor*density<<" max_length="<<maxLength<<" (quadrature order limit)"<<std::endl; continue; }
            network.fieldEvaluator.rebuild(factor*density,density,model::DislocationTargetPurpose::glide);
            const auto& snapshot=network.fieldEvaluator.buffers();
            device.uploadSources(snapshot.sources());
            const auto particles=[&]
            {
                model::evaluateCaiParticles(device.executionSpace(),kernel,device.sourcePositions(),device.sourceTensors(),device.targetPositions(),device.stresses(),device.sourceCount(),device.targetCount());
                device.executionSpace().fence();
            };
            particles();
            const auto measure=[](const auto& operation)
            {
                const auto start=Clock::now(); operation();
                return std::chrono::duration<double,std::milli>(Clock::now()-start).count();
            };
            std::vector<double> segmentTimes,particleTimes;
            for(int repetition=0;repetition<7;++repetition)
            {
                if(repetition%2) { particleTimes.push_back(measure(particles)); segmentTimes.push_back(measure(analytical)); }
                else { segmentTimes.push_back(measure(analytical)); particleTimes.push_back(measure(particles)); }
            }
            std::sort(segmentTimes.begin(),segmentTimes.end()); std::sort(particleTimes.begin(),particleTimes.end());
            std::cout<<"CAI_TIMING frame="<<argv[2]<<" density="<<factor*density<<" threads="<<device.executionSpace().concurrency()
                     <<" segments="<<referenceSegments.size()<<" sources="<<device.sourceCount()<<" targets="<<count
                     <<" analytical_median_ms="<<segmentTimes[3]<<" particle_median_ms="<<particleTimes[3]
                     <<" ratio="<<particleTimes[3]/segmentTimes[3]<<std::endl;
            model::DislocationStressArrays<3> result;
            device.downloadStresses(result);
            double err2=0,maxAbs=0,nearErr2=0,farErr2=0;
            std::vector<double> relative;
            for(std::size_t t=0;t<count;++t)
            {
                const Eigen::Vector3d x(targets.positions.values[0][t],targets.positions.values[1][t],targets.positions.values[2][t]);
                Eigen::Matrix3d actual;
                for(int i=0;i<3;++i) { for(int j=0;j<3;++j) { actual(i,j)=result.stresses.values[3*i+j][t]; } }
                const Eigen::Matrix3d error=actual-reference[t];
                if(!error.allFinite()) { throw std::runtime_error("Nonfinite stress error"); }
                err2+=error.squaredNorm(); maxAbs=std::max(maxAbs,error.norm());
                relative.push_back(error.norm()/std::max(reference[t].norm(),1.e-30));
                Eigen::Matrix3d nearError=Eigen::Matrix3d::Zero();
                for(std::size_t l=0;l<segments.size();++l)
                {
                    const auto& ref=referenceSegments[l];
                    const Eigen::Vector3d chord=ref.P1-ref.P0;
                    const double u=std::clamp((x-ref.P0).dot(chord)/chord.squaredNorm(),0.0,1.0);
                    if((x-ref.P0-u*chord).norm()>10*model::DislocationFieldBase<3>::a) { continue; }
                    auto key=std::make_pair(segments[l]->source->sID,segments[l]->sink->sID);
                    if(key.second<key.first) { std::swap(key.first,key.second); }
                    const auto range=snapshot.segments().at(key).sources;
                    double sum[6]={};
                    for(std::size_t p=range.begin;p<range.end();++p)
                    {
                        double separation[3],q[9];
                        for(int i=0;i<3;++i) { separation[i]=x[i]-snapshot.sources().positions.values[i][p]; }
                        for(int i=0;i<9;++i) { q[i]=snapshot.sources().weightedBurgersTensors.values[i][p]; }
                        kernel.addSymmetric(separation,q,sum);
                    }
                    Eigen::Matrix3d value; int k=0;
                    for(int i=0;i<3;++i) { for(int j=i;j<3;++j) { value(i,j)=value(j,i)=sum[k++]; } }
                    nearError+=value-ref.stress(x);
                }
                nearErr2+=nearError.squaredNorm(); farErr2+=(error-nearError).squaredNorm();
            }
            std::sort(relative.begin(),relative.end());
            std::cout<<std::setprecision(8)<<"CAI_ACCURACY frame="<<argv[2]<<" density="<<factor*density
                     <<" core="<<model::DislocationFieldBase<3>::a<<" sources="<<device.sourceCount()<<" targets="<<count
                     <<" relative_L2="<<std::sqrt(err2/ref2)<<" median_relative="<<relative[count/2]
                     <<" p95_relative="<<relative[std::min(count-1,count*95/100)]<<" max_relative="<<relative.back()
                     <<" rms_absolute="<<std::sqrt(err2/count)<<" max_absolute="<<maxAbs
                     <<" near_error_over_total_ref="<<std::sqrt(nearErr2/ref2)<<" far_error_over_total_ref="<<std::sqrt(farErr2/ref2)<<std::endl;
        }

    }
    catch(const std::exception& error) { std::cerr<<error.what()<<std::endl; return 1; }
}
