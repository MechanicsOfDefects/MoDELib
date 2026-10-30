#ifndef model_SolvedQuantityRegistry_H_
#define model_SolvedQuantityRegistry_H_

#include <Eigen/Core>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace model
{
    struct SolvedQuantityOptions
    {
        size_t maxIterations=100;
        double relativeTolerance=1.e-3;
        double relaxation=1.0;
        void validate() const;
    };

    struct SolvedQuantityContext
    {
        double time=0.0;
        size_t iteration=0;
    };

    struct SolvedQuantityDependency
    {
        std::string quantity;
        std::string provider; // Empty: all providers, excluding the current operation.
        bool required=false;
        std::function<bool()> enabled;
    };

    // One update may publish several quantities, but is executed only once.
    struct SolvedQuantityContribution
    {
        std::string provider;
        std::string name;
        std::vector<std::string> quantities;
        std::vector<SolvedQuantityDependency> dependencies;
        std::function<void(const SolvedQuantityContext&)> update;
        std::function<Eigen::VectorXd()> state;
        std::function<void(const Eigen::VectorXd&)> restore;
        // Absolute tolerances in state units, one per component (species included).
        std::function<Eigen::VectorXd()> absoluteTolerance;
    };

    struct SolvedQuantityBase
    {
        std::string name;
        std::vector<std::shared_ptr<SolvedQuantityContribution>> contributions;
    };

    struct SolvedQuantityGroup
    {
        std::vector<std::shared_ptr<SolvedQuantityContribution>> contributions;
        bool cyclic=false;
        size_t iterations=0;
        double error=0.0;
        void solve(const SolvedQuantityOptions&,double time);
    };

    class SolvedQuantityRegistry
    {
        std::map<std::string,std::shared_ptr<SolvedQuantityContribution>> operations;
        bool finalized=false;
    public:
        std::map<std::string,SolvedQuantityBase> quantities;
        std::vector<SolvedQuantityGroup> groups;
        SolvedQuantityOptions options;
        void clear();
        void add(SolvedQuantityContribution);
        void finalize();
        void solve(double time);
    };
}
#endif
