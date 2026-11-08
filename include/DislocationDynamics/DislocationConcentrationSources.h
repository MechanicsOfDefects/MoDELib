#ifndef model_DislocationConcentrationSources_H_
#define model_DislocationConcentrationSources_H_
#include <DislocationDynamics/ScreenedDislocationAssembly.h>
#include <map>
#include <limits>
namespace model
{
    template<class Network,class Parameters>
    std::vector<ClimbKernelChannel> concentrationSources(const Network& network,const Parameters& cdp,bool excludeBoundary)
    {
        std::vector<ClimbKernelChannel> channels;
        std::map<std::pair<std::size_t,std::size_t>,std::size_t> ids;
        for(const auto& item:network.networkLinks())
        {
            const auto link=item.second.lock();
            if(!link)
            {
                throw std::logic_error("Expired concentration source");
            }
            if(link->hasZeroBurgers()||link->grains().size()!=1||link->chordLength()<=std::numeric_limits<float>::epsilon())
            {
                continue;
            }
            if(excludeBoundary&&(link->isBoundarySegment()||link->isGrainBoundarySegment()))
            {
                continue;
            }
            const auto grain=(*link->grains().begin())->grainID;
            Eigen::Array<double,Parameters::mSize,1> bias=Eigen::Array<double,Parameters::mSize,1>::Ones();
            if(link->loopLinks().size()==1)
            {
                const auto loop=(*link->loopLinks().begin())->loop;
                const double orientation=loop->burgers().dot(loop->rightHandedUnitNormal());
                if(orientation>std::numeric_limits<float>::epsilon())
                {
                    bias=cdp.discreteDislocationBias.row(0).transpose();
                }
                else if(orientation<-std::numeric_limits<float>::epsilon())
                {
                    bias=cdp.discreteDislocationBias.row(1).transpose();
                }
            }
            const Eigen::Vector3d bxt=link->burgers().cross(link->unitDirection());
            for(int species=0;species<Parameters::mSize;++species)
            {
                if(!std::isfinite(bias(species))||bias(species)<=0)
                {
                    throw std::invalid_argument("Invalid concentration source bias");
                }
                const auto key=std::make_pair(std::size_t(grain),std::size_t(species));
                auto found=ids.find(key);
                if(found==ids.end())
                {
                    ClimbKernelChannel ch;ch.species=species;ch.grain=grain;ch.diffusion=cdp.D.at(grain).at(species);
                    found=ids.emplace(key,channels.size()).first;channels.push_back(std::move(ch));
                }
                ClimbSourceSegment segment;
                segment.p0=link->source->get_P();
                segment.p1=link->sink->get_P();
                segment.nodes={link->source->gID(),link->sink->gID()};
                segment.weights={bxt.dot(link->source->climbDirection())/bias(species),
                                 bxt.dot(link->sink->climbDirection())/bias(species)};
                channels[found->second].segments.push_back(segment);
            }
        }
        return channels;
    }
}
#endif
