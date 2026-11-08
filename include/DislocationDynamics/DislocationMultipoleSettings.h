#include <DislocationDynamics/CoreType.h>
#ifndef model_DislocationMultipoleSettings_H_
#define model_DislocationMultipoleSettings_H_
#include <FMM/FMMContext.h>
#include <TextFileParser.h>
#include <cmath>
#include <sstream>
#include <stdexcept>
namespace model
{
    // Physics-owned numerical settings; the physical core remains DD.txt:coreSize.
    struct DislocationMultipoleSettings
    {
        CoreType coreType=CoreType::Helmholtz;
        int ewaldSegmentChebyshevOrder=8,ewaldSegmentTableResolution=64;
        double ewaldWidth=0,ewaldWidthFactor=.125,sourceDensity=.2,preconditionerCutoff=64;
        DislocationMultipoleSettings(const std::string& file,MultipoleMode mode,double boxSize)
        {
            TextFileParser parser(file);
            for(const auto* key:{"fmmCoreWidth","fmmScreenWidth","fmmStressCore"}) {
      if(!parser.readStringVector(key).empty()) {
        throw std::invalid_argument(std::string("Retired DD setting ")+key+"; use coreType, coreSize and ewaldWidthFactor in DD.txt");
      }
    }
    if(!parser.readStringVector("coreType").empty())
            {
                const auto core=parser.readString("coreType",true);
                if(core!="Cai" && core!="Helmholtz") {
        throw std::invalid_argument("coreType must be Cai or Helmholtz");
      }
      coreType=core=="Cai"?CoreType::Cai:CoreType::Helmholtz;
            }
            if(coreType==CoreType::Cai && mode==MultipoleMode::ewaldFMM) {
      throw std::invalid_argument("coreType=Cai supports none and FMM; Cai Ewald screening is not implemented");
    }
    const auto read=[&](const char* key,double fallback)
            {
                const auto values=parser.readStringVector(key);
                if(values.empty()) {
        return fallback;
      }
      double value;std::string extra;std::istringstream stream(values.front().first);
                if(values.size()!=1||!(stream>>value)||(stream>>extra)||!std::isfinite(value)||value<=0) {
        throw std::invalid_argument(std::string("DD.txt requires a positive finite ")+key);
      }
      return value;
            };
            const double core=parser.readScalar<double>("coreSize",false);
            if(!std::isfinite(core)||core<=0) {
      throw std::invalid_argument("coreSize must be positive and finite");
    }
    ewaldWidth=4*core; // Unscreened families retain a valid, unused screening parameter.
            sourceDensity=read("fmmSourcePointsPerLength",.2);
                const double order=read("ewaldSegmentChebyshevOrder",8),resolution=read("ewaldSegmentTableResolution",64);
                if(order!=std::floor(order)||order<4||order>16||resolution!=std::floor(resolution)||resolution<4||resolution>256) {
      throw std::invalid_argument("ewaldSegmentChebyshevOrder must be integer 4..16 and ewaldSegmentTableResolution integer 4..256");
    }
    ewaldSegmentChebyshevOrder=int(order);ewaldSegmentTableResolution=int(resolution);
            if(mode==MultipoleMode::ewaldFMM)
            {
                for(const auto* key:{"fmmEwaldWidth","fmmBoxWidthFactor","fmmShortCutoff"}) {
        if(!parser.readStringVector(key).empty()) {
          throw std::invalid_argument(std::string("Retired ewaldFMM setting ")+key+
                            "; set fmmTreeDepth in DefectiveCrystal.txt and ewaldWidthFactor in DD.txt; h=smallest leaf width*ewaldWidthFactor");
        }
      }
      ewaldWidthFactor=read("ewaldWidthFactor",.125);
                if(ewaldWidthFactor>=1) {
        throw std::invalid_argument("ewaldWidthFactor must satisfy 0 < factor < 1");
      }
      if(!std::isfinite(boxSize)||boxSize<=0) {
        throw std::invalid_argument("Ewald box size must be positive and finite");
      }
      ewaldWidth=boxSize*ewaldWidthFactor;
                if(!std::isfinite(ewaldWidth)||ewaldWidth<4*core) {
        throw std::invalid_argument("Derived Ewald width (smallest leaf width*ewaldWidthFactor) must be >=4*coreSize; increase ewaldWidthFactor or reduce fmmTreeDepth");
      }
    }
            else {
      preconditionerCutoff=read("climbPreconditionerCutoff",64);
    }
  }
    };
}
#endif
