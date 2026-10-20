#include <DefectiveCrystal/DefectiveCrystal.h>
#include <MicrostructureGenerator.h>
#include <iostream>

int main(int argc,char** argv)
{
    try
    {
        if(argc!=2)
        {
            throw std::runtime_error("Expected the spatialCDtest report directory.");
        }
        model::DefectiveCrystalBase<3> base(argv[1]);
        model::MicrostructureGenerator generator(base);
        generator.readMicrostructureFile();
        generator.writeConfigFiles(0);
        model::DefectiveCrystalConfigIO<3> config(base.simulationParameters.traitsIO.configFolder);
        config.read(0);
        model::DefectiveCrystal<3> crystal(base);
        crystal.initializeConfiguration(config);
        crystal.runSteps();
    }
    catch(const std::exception& error)
    {
        std::cerr<<"spatialCDtest failed: "<<error.what()<<std::endl;
        return 1;
    }
    return 0;
}
