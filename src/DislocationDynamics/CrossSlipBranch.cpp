#include <CrossSlipBranch.h>

namespace model
{
    std::optional<ThermalCrossSlipEvent> sampleThermalCrossSlip(
        const std::vector<double>& lengths,
        const std::vector<std::map<int,double>>& rates,
        const std::map<int,double>& nucleusLengths,
        const double dt, const double occurrence, const double targetDraw, const double positionDraw)
    {
        if(lengths.size()!=rates.size() || !std::isfinite(dt) || dt<0.0)
        {
            throw std::runtime_error("Invalid thermal cross-slip sampling data");
        }
        for(const double draw : {occurrence,targetDraw,positionDraw})
        {
            if(!std::isfinite(draw) || draw<0.0 || draw>=1.0)
            {
                throw std::runtime_error("Thermal cross-slip uniforms must lie in [0,1)");
            }
        }
        double length=0.0;
        for(const double value : lengths)
        {
            if(!std::isfinite(value) || value<0.0)
            {
                throw std::runtime_error("Invalid thermal branch length");
            }
            length+=value;
        }
        std::map<int,double> totals;
        for(const auto& segment : rates)
        {
            for(const auto& [target,rate] : segment)
            {
                const double nucleus=nucleusLengths.at(target);
                if(!std::isfinite(rate) || rate<0.0 || !std::isfinite(nucleus) || nucleus<=0.0)
                {
                    throw std::runtime_error("Invalid thermal cross-slip rate or nucleation length");
                }
                if(length>=nucleus)
                {
                    totals[target]+=rate;
                }
            }
        }
        double total=0.0;
        for(const auto& [target,rate] : totals)
        {
            total+=rate;
        }
        if(!std::isfinite(total))
        {
            throw std::runtime_error("Thermal cross-slip total rate overflow");
        }
        if(total<=0.0 || occurrence>=-std::expm1(-total*dt))
        {
            return std::nullopt;
        }
        double selection=targetDraw*total;
        for(const auto& [target,rate] : totals)
        {
            if(rate>0.0 && selection<rate)
            {
                double location=positionDraw*rate;
                double offset=0.0;
                for(size_t k=0;k<rates.size();++k)
                {
                    const auto found=rates[k].find(target);
                    const double local=found==rates[k].end()? 0.0 : found->second;
                    if(local>0.0 && location<local)
                    {
                        const double nucleus=nucleusLengths.at(target);
                        const double center=std::clamp(offset+lengths[k]*location/local,0.5*nucleus,length-0.5*nucleus);
                        return ThermalCrossSlipEvent{target,nucleus,center/length};
                    }
                    location-=local;
                    offset+=lengths[k];
                }
            }
            selection-=rate;
        }
        throw std::runtime_error("Thermal cross-slip selection failed");
    }

 template<int dim>
 CrossSlipBranch<dim>::CrossSlipBranch(std::deque<std::shared_ptr<LoopNodeType>> nodes,int target) :
     first(std::move(nodes)),second(target)
 {
 }

 template struct CrossSlipBranch<3>;
}
