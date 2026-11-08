/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_KokkosDislocationFieldBuffers_H_
#define model_KokkosDislocationFieldBuffers_H_
// Load full floating-point limits before SDK math headers request partial ones.
#include <cfloat>
#include <DislocationDynamics/DislocationFieldBuffers.h>
#include <Kokkos_Core.hpp>
#include <string>

namespace model
{
    // ExecutionSpace is intentionally a template parameter: CPU and accelerator
    // builds use identical kernels. LayoutLeft makes each component contiguous
    // over points, matching the CPU structure-of-arrays representation.
    // The application owns Kokkos::ScopeGuard; destroy all buffers before it.
    template<int dim, class ExecutionSpace=Kokkos::DefaultExecutionSpace>
    class KokkosDislocationFieldBuffers
    {
        static_assert(dim>0,"Positive dimension required");
    public:
        using View=Kokkos::View<double**,Kokkos::LayoutLeft,ExecutionSpace>;
        using ConstView=Kokkos::View<const double**,Kokkos::LayoutLeft,ExecutionSpace>;
        // Infer the mirror type across Kokkos 4/5 alias naming changes.
        using HostView=decltype(Kokkos::create_mirror(std::declval<View>()));
    private:
        ExecutionSpace execution_;
        View sourcePositions_, sourceTensors_, targetPositions_, stresses_;
        HostView hostSourcePositions_, hostSourceTensors_, hostTargetPositions_, hostStresses_;
        std::size_t sourceCount_=0, targetCount_=0;

        static ExecutionSpace checkedExecutionSpace()
        {
            if(!Kokkos::is_initialized() || Kokkos::is_finalized())
            {
                throw std::logic_error("Initialize Kokkos before creating field buffers");
            }
            return ExecutionSpace{};
        }
        template<int components>
        static void checkShape(const DislocationComponentArrays<components>& data, std::size_t count)
        {
            for(const auto& component : data.values)
            {
                if(component.size()!=count)
                {
                    throw std::invalid_argument("Inconsistent source/target component lengths");
                }
            }
        }
        void allocate(View& device, HostView& host, const char* name,
                      std::size_t count, std::size_t components)
        {
            if(count>device.extent(0))
            {
                execution_.fence("Before growing field buffer");
                const auto capacity=std::max(count,2*device.extent(0));
                device=View(std::string(name),capacity,components);
                // A separate staging mirror even on host/unified-memory backends.
                host=Kokkos::create_mirror(device);
            }
        }
        template<int components>
        void upload(const DislocationComponentArrays<components>& data,
                    const View& device, const HostView& host, std::size_t count)
        {
            if(count==0)
            {
                return;
            }
            for(int c=0;c<components;++c)
            {
                for(std::size_t q=0;q<count;++q)
                {
                    host(q,c)=data.values[c][q];
                }
            }
            const auto range=std::make_pair(std::size_t(0),count);
            Kokkos::deep_copy(execution_,Kokkos::subview(device,range,Kokkos::ALL()),
                             Kokkos::subview(host,range,Kokkos::ALL()));
        }
    public:
        KokkosDislocationFieldBuffers():execution_(checkedExecutionSpace()) {}
        KokkosDislocationFieldBuffers(const KokkosDislocationFieldBuffers&)=delete;
        KokkosDislocationFieldBuffers& operator=(const KokkosDislocationFieldBuffers&)=delete;
        const ExecutionSpace& executionSpace() const { return execution_; }
        std::size_t sourceCount() const { return sourceCount_; }
        std::size_t targetCount() const { return targetCount_; }
        std::size_t sourceCapacity() const { return sourcePositions_.extent(0); }
        std::size_t targetCapacity() const { return targetPositions_.extent(0); }
        // Views may have spare capacity. Kernels must use the active counts.
        ConstView sourcePositions() const { return sourcePositions_; }
        ConstView sourceTensors() const { return sourceTensors_; }
        ConstView targetPositions() const { return targetPositions_; }
        View stresses() const { return stresses_; }

        void resetResults()
        {
            if(targetCount_!=0)
            {
                const auto range=std::make_pair(std::size_t(0),targetCount_);
                Kokkos::deep_copy(execution_,Kokkos::subview(stresses_,range,Kokkos::ALL()),0.0);
            }
        }
        // Explicit source update allows many target batches to reuse sources.
        // Uploads are blocking so callers can safely reuse staging/CPU storage.
        void uploadSources(const DislocationSourceArrays<dim>& source)
        {
            const auto count=source.size();
            checkShape(source.positions,count);
            checkShape(source.weightedBurgersTensors,count);
            execution_.fence("Before source staging update");
            allocate(sourcePositions_,hostSourcePositions_,"DD source positions",count,dim);
            allocate(sourceTensors_,hostSourceTensors_,"DD source tensors",count,dim*dim);
            upload(source.positions,sourcePositions_,hostSourcePositions_,count);
            upload(source.weightedBurgersTensors,sourceTensors_,hostSourceTensors_,count);
            sourceCount_=count;
            resetResults();
            execution_.fence("Source upload complete");
        }
        void uploadTargets(const DislocationTargetArrays<dim>& target)
        {
            const auto count=target.size();
            checkShape(target.positions,count);
            execution_.fence("Before target staging update");
            allocate(targetPositions_,hostTargetPositions_,"DD target positions",count,dim);
            allocate(stresses_,hostStresses_,"DD stress results",count,dim*dim);
            upload(target.positions,targetPositions_,hostTargetPositions_,count);
            targetCount_=count;
            resetResults();
            execution_.fence("Target upload complete");
        }
        void downloadStresses(DislocationStressArrays<dim>& result)
        {
            result.reset(targetCount_);
            if(targetCount_==0)
            {
                return;
            }
            const auto range=std::make_pair(std::size_t(0),targetCount_);
            Kokkos::deep_copy(execution_,Kokkos::subview(hostStresses_,range,Kokkos::ALL()),
                             Kokkos::subview(stresses_,range,Kokkos::ALL()));
            execution_.fence("Stress download complete");
            for(int c=0;c<dim*dim;++c)
            {
                for(std::size_t q=0;q<targetCount_;++q)
                {
                    result.stresses.values[c][q]=hostStresses_(q,c);
                }
            }
        }
    };
}
#endif
