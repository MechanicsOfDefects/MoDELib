#ifndef model_CrossSlipBranch_H
#define model_CrossSlipBranch_H
#include <algorithm>
#include <deque>
#include <memory>
#include <optional>
#include <map>
#include <vector>
#include <cmath>
#include <stdexcept>
namespace model
{
    template<int dim> class DislocationLoopNode;
    struct ThermalCrossSlipEvent
    {
        int target;
        double length;
        double centerFraction;
    };

    // Rates are integrated over physical line length, never over node count.
    // The three uniforms select occurrence, competing target, and location.
    std::optional<ThermalCrossSlipEvent> sampleThermalCrossSlip(
        const std::vector<double>& lengths,
        const std::vector<std::map<int,double>>& rates,
        const std::map<int,double>& nucleusLengths,
        const double dt, const double occurrence, const double targetDraw, const double positionDraw);


    template<int dim>
    struct CrossSlipBranch
    {
        using LoopNodeType=DislocationLoopNode<dim>;
        std::deque<std::shared_ptr<LoopNodeType>> first;
        int second;
        std::optional<ThermalCrossSlipEvent> thermal;
        CrossSlipBranch(std::deque<std::shared_ptr<LoopNodeType>> nodes,int target);
    };
}
#endif
