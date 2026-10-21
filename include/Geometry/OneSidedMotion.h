#ifndef model_OneSidedMotion_h_
#define model_OneSidedMotion_h_
#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <cmath>

namespace model
{
    // Closest motion satisfying equality normals and outward half-spaces.
    // In dimension d the optimum needs at most d independent active faces.
    template<int dim>
    Eigen::Matrix<double,dim,1> oneSidedMotion(
        const Eigen::Matrix<double,dim,1>& requested,
        const std::vector<Eigen::Matrix<double,dim,1>>& equalities,
        const std::vector<Eigen::Matrix<double,dim,1>>& outward)
    {
        using V=Eigen::Matrix<double,dim,1>;
        std::vector<V> basis;
        const auto addNormal=[](std::vector<V>& normals,V n)
        {
            for(const auto& q : normals)
            {
                n-=n.dot(q)*q;
            }
            const double norm=n.norm();
            if(norm>1.e-12)
            {
                normals.push_back(n/norm);
                return true;
            }
            return false;
        };
        for(const auto& n : equalities)
        {
            addNormal(basis,n);
        }
        V best=V::Zero(); // Always feasible for homogeneous constraints.
        double cost=requested.squaredNorm();
        const double tolerance=1.e-12*requested.norm();
        std::function<void(size_t,const std::vector<V>&)> visit;
        visit=[&](size_t start,const std::vector<V>& active)
        {
            V candidate=requested;
            for(const auto& n : active)
            {
                candidate-=candidate.dot(n)*n;
            }
            bool feasible=true;
            for(const auto& n : outward)
            {
                if(candidate.dot(n)<-tolerance)
                {
                    feasible=false;
                    break;
                }
            }
            const double candidateCost=(candidate-requested).squaredNorm();
            if(feasible && candidateCost<=cost)
            {
                best=candidate;
                cost=candidateCost;
            }
            if(active.size()<dim)
            {
                for(size_t i=start;i<outward.size();++i)
                {
                    auto next=active;
                    if(addNormal(next,outward[i]))
                    {
                        visit(i+1,next);
                    }
                }
            }
        };
        visit(0,basis);
        return best;
    }
}
#endif
