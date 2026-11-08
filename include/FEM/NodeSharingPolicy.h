#ifndef model_NodeSharingPolicy_H_
#define model_NodeSharingPolicy_H_
namespace model
{
    // Connectivity only; shape functions and weak formulations are unchanged.
    enum class NodeSharingPolicy { continuous, grainContinuous, discontinuous };
}
#endif
