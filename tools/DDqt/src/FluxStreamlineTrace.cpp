#include <FluxStreamlineTrace.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace model
{
    FluxTrace traceFlux(const SimplicialMesh<3>& mesh,const Eigen::Vector3d& seed,
                        const FluxEvaluator& evaluate,const FluxTraceOptions& options,int direction)
    {
        if(!(options.step>0.0) || !(options.length>0.0) || options.minimumFlux<0.0 ||
           !std::isfinite(options.step) || !std::isfinite(options.length) ||
           !std::isfinite(options.minimumFlux) || options.maxSteps<1 || std::abs(direction)!=1)
        {
            throw std::invalid_argument("Invalid flux tracing controls");
        }
        FluxTrace result;
        const auto found=mesh.search(seed);
        if(!found.first)
        {
            return result;
        }
        const Simplex<3,3>* cell=found.second;
        Eigen::Vector3d x=seed;
        double travelled=0.0;
        result.emplace_back();
        for(int iteration=0;iteration<options.maxSteps && travelled<options.length;++iteration)
        {
            const auto sample=evaluate(x,*cell);
            const double speed=sample.flux.stableNorm();
            if(!std::isfinite(speed) || speed<=options.minimumFlux)
            {
                break;
            }
            result.back().push_back({x,speed});
            const double h=std::min(options.step,options.length-travelled);
            Eigen::Vector3d tangent=direction*sample.flux/speed;
            const Eigen::Vector3d midpoint=x+0.5*h*tangent;
            if(cell->pos2bary(midpoint).minCoeff()>1.e-10)
            {
                const auto mid=evaluate(midpoint,*cell);
                const double midSpeed=mid.flux.stableNorm();
                if(!std::isfinite(midSpeed) || midSpeed<=options.minimumFlux)
                {
                    break;
                }
                tangent=direction*mid.flux/midSpeed;
            }
            const Eigen::Vector3d target=x+h*tangent;
            const auto bary=cell->pos2bary(x);
            const auto endBary=cell->pos2bary(target);
            double fraction=1.0;
            int exitFace=-1;
            for(int f=0;f<4;++f)
            {
                if(endBary(f)<0.0 && endBary(f)<bary(f))
                {
                    const double candidate=std::max(0.0,bary(f)/(bary(f)-endBary(f)));
                    if(candidate<fraction)
                    {
                        fraction=candidate;
                        exitFace=f;
                    }
                }
            }
            const Eigen::Vector3d hit=x+fraction*h*tangent;
            travelled+=(hit-x).norm();
            if(exitFace<0)
            {
                x=target;
                if(travelled>=options.length)
                {
                    result.back().push_back({x,speed});
                }
                continue;
            }
            result.back().push_back({hit,speed});
            const auto& face=cell->child(exitFace);
            const Simplex<3,3>* next=nullptr;
            Eigen::Vector3d shift=Eigen::Vector3d::Zero();
            if(face.isRegionBoundarySimplex())
            {
                break;
            }
            if(face.isBoundarySimplex())
            {
                const auto pair=mesh.periodicSimplexPairs().find(&face);
                if(pair==mesh.periodicSimplexPairs().end())
                {
                    break;
                }
                next=pair->second.partner->parents().begin()->second.get();
                if(next->region->regionID!=cell->region->regionID)
                {
                    break;
                }
                shift=pair->second.shift;
                result.emplace_back();
            }
            else
            {
                for(const auto& parent : face.parents())
                {
                    if(parent.second.get()!=cell)
                    {
                        next=parent.second.get();
                        break;
                    }
                }
            }
            if(!next)
            {
                break;
            }
            x=hit+shift;
            // Move just inside the adjacent tetrahedron to avoid recrossing
            // its face through roundoff; never search across a grain boundary.
            x+=1.e-9*(next->center()-x);
            cell=next;
        }
        return result;
    }
}
