#include <CrossSlipBranch.h>
#include <random>
#include <iostream>
#include <stdexcept>

void require(bool value,const char* message)
{
    if(!value)
    {
        throw std::runtime_error(message);
    }
}

int main()
{
    try
    {
        const std::map<int,double> nuclei{{2,20.0},{7,30.0}};
        std::mt19937 gen(42);
        std::uniform_real_distribution<double> uniform(0.0,1.0);
        size_t events=0,second=0;
        for(size_t trial=0;trial<100000;++trial)
        {
            const double u=uniform(gen),v=uniform(gen),w=uniform(gen);
            const auto reference=model::sampleThermalCrossSlip({100.0},{{{2,2.0},{7,6.0}}},nuclei,0.1,u,v,w);
            for(const size_t count : {2,10,100})
            {
                const std::vector<double> lengths(count,100.0/count);
                const std::vector<std::map<int,double>> rates(count,{{2,2.0/count},{7,6.0/count}});
                const auto result=model::sampleThermalCrossSlip(lengths,rates,nuclei,0.1,u,v,w);
                require(bool(result)==bool(reference),"Segmentation changed event occurrence");
                if(result)
                {
                    require(result->target==reference->target,"Segmentation changed target competition");
                    require(std::abs(result->centerFraction-reference->centerFraction)<1.e-12,"Segmentation changed nucleation position");
                    require(result->length==reference->length,"Segmentation changed nucleus length");
                }
            }
            if(reference)
            {
                ++events;
                second+=reference->target==7;
            }
        }
        require(std::abs(double(events)/100000.0-(-std::expm1(-0.8)))<0.005,"Incorrect total hazard");
        require(std::abs(double(second)/events-0.75)<0.01,"Incorrect target competition");
        require(!model::sampleThermalCrossSlip({10.0},{{{2,2.0}}},nuclei,100.0,0.0,0.0,0.0),"Nucleated on a branch shorter than physical nucleus");
        require(!model::sampleThermalCrossSlip({100.0},{{{2,2.0}}},nuclei,0.0,0.0,0.0,0.0),"Nucleated in zero time");
        require(!model::sampleThermalCrossSlip({100.0},{{{2,0.0}}},nuclei,1.0,0.0,0.0,0.0),"Nucleated with zero rate");
        for(const double v : {0.1,0.6})
        {
            const auto event=model::sampleThermalCrossSlip({100.0},{{{2,2.0},{7,2.0}}},nuclei,100.0,0.0,v,0.5);
            require(event && event->target==(v<0.5? 2 : 7),"Equal rates discarded a target");
        }
        // Nonuniform physical rate: subdivide each half without changing its density.
        const auto a=model::sampleThermalCrossSlip({40.0,60.0},{{{2,1.0}},{{2,9.0}}},nuclei,1.0,0.1,0.1,0.6);
        const auto b=model::sampleThermalCrossSlip({10.0,30.0,20.0,40.0},{{{2,0.25}},{{2,0.75}},{{2,3.0}},{{2,6.0}}},nuclei,1.0,0.1,0.1,0.6);
        require(a && b && std::abs(a->centerFraction-b->centerFraction)<1.e-12,"Nonuniform rates depend on segmentation");
        std::cout<<"PASS: thermal branch sampling independent of segmentation, competing hazards and physical nucleus lengths\n";
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
