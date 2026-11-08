#include <FrankLoopsDensitySpecification.h>
#include <cmath>
#include <numbers>
#include <iostream>
#include <limits>
#include <stdexcept>
void require(bool ok) { if(!ok) { throw std::runtime_error("Frank loop discretization regression"); } }
int main(int argc,char** argv)
{
    try
    {
        require(argc==2);
        const std::string root(argv[1]);
        model::FrankLoopsDensitySpecification fixed(root+"/Library/Microstructures/frankLoopsDensity.txt");
        require(fixed.targetSegmentLength==0 && fixed.sidesForRadius(100)==20);
        model::FrankLoopsDensitySpecification annealing(root+"/tutorials/annealing/inputFiles/microstructureFile1.txt");
        require(annealing.targetSegmentLength==50 && annealing.sidesForRadius(100)==13);
        model::FrankLoopsDensitySpecification spec;
        spec.numberOfSides=20;
        require(spec.sidesForRadius(10)==20 && spec.sidesForRadius(1000)==20);
        spec.targetSegmentLength=50;
        require(spec.sidesForRadius(10)==8);
        require(spec.sidesForRadius(100)==13);
        require(spec.sidesForRadius(1000)==126);
        for(double radius:{1.,10.,40.,100.,1000.})
        {
            const auto n=spec.sidesForRadius(radius);
            require(n>=8 && 2*radius*std::sin(std::numbers::pi/n)<=50*(1+1e-14));
        }
        for(double bad:{-1.,std::numeric_limits<double>::quiet_NaN()})
        {
            spec.targetSegmentLength=bad; bool rejected=false;
            try { spec.validateDiscretization(); } catch(const std::invalid_argument&) { rejected=true; }
            require(rejected);
        }
        spec.targetSegmentLength=50; bool rejected=false;
        try { spec.sidesForRadius(0); } catch(const std::invalid_argument&) { rejected=true; }
        require(rejected);
        std::cout<<"Frank loop discretization tests passed\n";
    }
    catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
