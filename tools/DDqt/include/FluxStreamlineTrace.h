#ifndef model_FluxStreamlineTrace_H_
#define model_FluxStreamlineTrace_H_
#include <SimplicialMesh.h>
#include <functional>
#include <vector>
namespace model
{
    struct FluxSample
    {
        Eigen::Vector3d flux=Eigen::Vector3d::Zero();
    };
    struct FluxTracePoint
    {
        Eigen::Vector3d position;
        double magnitude;
    };
    struct FluxTraceOptions
    {
        double step=10.0;
        double length=1000.0;
        double minimumFlux=0.0;
        int maxSteps=2000;
    };
    using FluxEvaluator=std::function<FluxSample(const Eigen::Vector3d&,const Simplex<3,3>&)>;
    using FluxTrace=std::vector<std::vector<FluxTracePoint>>;
    // Integrate in arclength, following simplex connectivity at every face.
    // Each periodic wrap starts a separate drawable polyline.
    FluxTrace traceFlux(const SimplicialMesh<3>& mesh,const Eigen::Vector3d& seed,
                        const FluxEvaluator& evaluate,const FluxTraceOptions& options,int direction);
}
#endif
