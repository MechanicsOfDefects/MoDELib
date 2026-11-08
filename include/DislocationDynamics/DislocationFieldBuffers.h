/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_DislocationFieldBuffers_H_
#define model_DislocationFieldBuffers_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace model
{
    enum class DislocationTargetPurpose { glide, climb, diagnostic };

    // Plain numerical records: no ownership of, or pointers into, the network.
    template<int dim>
    struct DislocationQuadratureGeometry
    {
        std::array<double,dim> position{};
        std::array<double,dim> tangent{}; // Unit tangent in network-link orientation.
        double parameter=0.0;
        double lineWeight=0.0; // Quadrature weight times arc-length Jacobian.
    };

    template<int dim>
    struct DislocationSourcePoint
    {
        DislocationQuadratureGeometry<dim> geometry;
        std::array<double,dim*dim> weightedBurgersTensor{}; // Row-major b_i t_j dL.
    };

    template<int dim>
    struct DislocationFieldPoint
    {
        DislocationQuadratureGeometry<dim> geometry;
        DislocationTargetPurpose purpose=DislocationTargetPurpose::diagnostic;
    };

    struct DislocationPointRange
    {
        std::size_t begin=0;
        std::size_t size=0;
        std::size_t end() const { return begin+size; }
    };

    struct DislocationSegmentFieldSelection
    {
        bool glissile=false;
        bool physicalBoundary=false;
        bool grainBoundary=false;
        bool accepts(DislocationTargetPurpose purpose) const
        {
            if(purpose==DislocationTargetPurpose::diagnostic)
            {
                return true;
            }
            if(physicalBoundary || grainBoundary)
            {
                return false;
            }
            return purpose==DislocationTargetPurpose::glide ? glissile : !glissile;
        }
    };

    template<int components>
    struct DislocationComponentArrays
    {
        std::array<std::vector<double>,components> values;
        void reserve(std::size_t n)
        {
            for(auto& component : values)
            {
                if(n>component.capacity())
                {
                    component.reserve(std::max(n,component.capacity()*2));
                }
            }
        }
        void clear()
        {
            for(auto& component : values)
            {
                component.clear();
            }
        }
        void append(const std::array<double,components>& value)
        {
            for(int i=0;i<components;++i)
            {
                values[i].push_back(value[i]);
            }
        }
    };

    template<int dim>
    struct DislocationSourceArrays
    {
        DislocationComponentArrays<dim> positions;
        DislocationComponentArrays<dim*dim> weightedBurgersTensors;
        std::size_t size() const { return positions.values[0].size(); }
    };

    template<int dim>
    struct DislocationTargetArrays
    {
        DislocationComponentArrays<dim> positions;
        std::size_t size() const { return positions.values[0].size(); }
    };

    template<int dim>
    struct DislocationStressArrays
    {
        // Full row-major tensor; no assumptions about a backend's packing.
        DislocationComponentArrays<dim*dim> stresses;
        void reset(std::size_t count)
        {
            for(auto& component : stresses.values)
            {
                component.assign(count,0.0);
            }
        }
    };

    // An evaluator-owned snapshot. CPU metadata stays separate from GPU arrays.
    // Build once from unique network links, not once per loop link. Rebuild after
    // motion/topology changes; clear() retains allocations. No device API here.
    template<int dim>
    class DislocationFieldBuffers
    {
        static_assert(dim>0,"Positive dimension required");
    public:
        using Geometry=DislocationQuadratureGeometry<dim>;
        using Vector=std::array<double,dim>;
        using SegmentKey=std::pair<std::size_t,std::size_t>;
        struct SegmentRanges
        {
            DislocationPointRange sources;
            DislocationPointRange targets;
        };
        struct TargetBinding
        {
            SegmentKey segment{};
            std::size_t localIndex=0;
            DislocationFieldPoint<dim> point;
            bool external=false;
        };

    private:
        DislocationSourceArrays<dim> sources_;
        DislocationTargetArrays<dim> targets_;
        std::map<SegmentKey,SegmentRanges> segments_;
        std::vector<TargetBinding> bindings_;

        static void validate(const Geometry& point)
        {
            if(!std::isfinite(point.parameter) || point.parameter<0.0 || point.parameter>1.0
               || !std::isfinite(point.lineWeight) || point.lineWeight<0.0)
            {
                throw std::invalid_argument("Invalid dislocation quadrature parameter/weight");
            }
            double tangentNorm2=0.0;
            for(int i=0;i<dim;++i)
            {
                if(!std::isfinite(point.position[i]) || !std::isfinite(point.tangent[i]))
                {
                    throw std::invalid_argument("Non-finite dislocation quadrature geometry");
                }
                tangentNorm2+=point.tangent[i]*point.tangent[i];
            }
            if(std::abs(tangentNorm2-1.0)>1.0e-10)
            {
                throw std::invalid_argument("Dislocation quadrature tangent must be a unit vector");
            }
        }
        void appendTarget(const Vector& position, const TargetBinding& binding)
        {
            targets_.positions.append(position);
            bindings_.push_back(binding);
        }

    public:
        const DislocationSourceArrays<dim>& sources() const { return sources_; }
        const DislocationTargetArrays<dim>& targets() const { return targets_; }
        const std::map<SegmentKey,SegmentRanges>& segments() const { return segments_; }
        const std::vector<TargetBinding>& targetBindings() const { return bindings_; }

        void reserve(std::size_t sourceCount, std::size_t targetCount)
        {
            sources_.positions.reserve(sourceCount);
            sources_.weightedBurgersTensors.reserve(sourceCount);
            targets_.positions.reserve(targetCount);
            if(targetCount>bindings_.capacity())
            {
                bindings_.reserve(std::max(targetCount,bindings_.capacity()*2));
            }
        }
        void clear()
        {
            sources_.positions.clear();
            sources_.weightedBurgersTensors.clear();
            targets_.positions.clear();
            segments_.clear();
            bindings_.clear();
        }

        // Independent rules may supply different source and target samples.
        // Burgers and tangent must both use the supplied endpoint orientation.
        // Boundary flags describe the whole segment, not either endpoint.
        void addSegment(SegmentKey endpoints, const Vector& burgers,
                        const DislocationSegmentFieldSelection& selection,
                        const std::vector<Geometry>& sourcePoints,
                        const std::vector<Geometry>& fieldPoints,
                        DislocationTargetPurpose purpose)
        {
            const SegmentKey key=std::minmax(endpoints.first,endpoints.second);
            if(endpoints.first==endpoints.second || segments_.count(key))
            {
                throw std::invalid_argument("Duplicate or invalid network segment in field buffers");
            }
            bool nonzero=false;
            for(const double b : burgers)
            {
                if(!std::isfinite(b))
                {
                    throw std::invalid_argument("Non-finite Burgers vector");
                }
                nonzero=nonzero || b!=0.0;
            }
            // Validate before any mutation. Zero-Burgers links need no samples.
            if(nonzero)
            {
                if(sourcePoints.empty())
                {
                    throw std::invalid_argument("A nonzero-Burgers segment requires source samples");
                }
                for(const auto& point : sourcePoints)
                {
                    validate(point);
                }
                if(selection.accepts(purpose))
                {
                    for(const auto& point : fieldPoints)
                    {
                        validate(point);
                    }
                }
            }
            SegmentRanges ranges{{sources_.size(),0},{targets_.size(),0}};
            if(nonzero)
            {
                reserve(sources_.size()+sourcePoints.size(),targets_.size()+fieldPoints.size());
                for(const auto& point : sourcePoints)
                {
                    DislocationSourcePoint<dim> source{point,{}};
                    for(int i=0;i<dim;++i)
                    {
                        for(int j=0;j<dim;++j)
                        {
                            source.weightedBurgersTensor[i*dim+j]=burgers[i]*point.tangent[j]*point.lineWeight;
                        }
                    }
                    sources_.positions.append(source.geometry.position);
                    sources_.weightedBurgersTensors.append(source.weightedBurgersTensor);
                }
                if(selection.accepts(purpose))
                {
                    for(std::size_t q=0;q<fieldPoints.size();++q)
                    {
                        appendTarget(fieldPoints[q].position,{endpoints,q,{fieldPoints[q],purpose},false});
                    }
                }
            }
            ranges.sources.size=sources_.size()-ranges.sources.begin;
            ranges.targets.size=targets_.size()-ranges.targets.begin;
            segments_.emplace(key,ranges);
        }

        // External FEM/visualization points do not need tangent or line weight.
        std::size_t addExternalTarget(const Vector& position)
        {
            for(const double x : position)
            {
                if(!std::isfinite(x))
                {
                    throw std::invalid_argument("Non-finite external target");
                }
            }
            const auto index=targets_.size();
            reserve(sources_.size(),index+1);
            TargetBinding binding;
            binding.external=true;
            binding.point.geometry.position=position;
            appendTarget(position,binding);
            return index;
        }
    };
}
#endif
