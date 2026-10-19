#ifndef model_DDframeTime_H_
#define model_DDframeTime_H_

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace model
{
    struct DDframeTime
    {
        long int runID=0;
        double totalTime=0.0;
        double dt=10.0; // Same initial step as DefectiveCrystalParameters.

        static DDframeTime read(const std::string& historyPath,const size_t frameID)
        {
            std::ifstream history(historyPath);
            bool hasRecords=false;
            std::string line;
            while(std::getline(history,line))
            {
                std::stringstream row(line);
                DDframeTime entry;
                if(row>>entry.runID>>entry.totalTime>>entry.dt)
                {
                    hasRecords=true;
                    if(entry.runID>=0 && static_cast<size_t>(entry.runID)==frameID)
                    {
                        return entry;
                    }
                }
                else if(line.find_first_not_of(" \t\r")==std::string::npos)
                {
                    continue;
                }
                else
                {
                    throw std::runtime_error("Invalid frame-time record in "+historyPath+".");
                }
            }
            // A generated initial configuration predates the first history row.
            // Never infer the time of a missing later frame from its frame ID.
            if(frameID==0 && !hasRecords)
            {
                return DDframeTime();
            }
            throw std::runtime_error("Cannot find frame "+std::to_string(frameID)
                                     +" in history file "+historyPath+".");
        }
    };
}
#endif
