#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
#include <cfloat>
#include <Kokkos_Core.hpp>
#endif
#include <FMM/FMMContext.h>
#include <TextFileParser.h>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <iostream>
namespace model
{
    namespace
    {
        double setting(const std::string& file,const std::string& name,double fallback)
        {
            const auto values=TextFileParser(file).readStringVector(name);
            if(values.empty()) {
    return fallback;
  }
  double value;std::string extra;std::istringstream stream(values.front().first);
            if(values.size()!=1||!(stream>>value)||(stream>>extra)||!std::isfinite(value)) {
    throw std::invalid_argument("Invalid FMM setting: "+name);
  }
  return value;
        }
        int order(const std::string& file)
        {const auto value=setting(file,"fmmInterpolationOrder",4);if(value<2||value>6||std::floor(value)!=value) {
    throw std::invalid_argument("fmmInterpolationOrder must be an integer 2..6");
  }
  return static_cast<int>(value);}
        MultipoleMode readMode(const std::string& file)
        {
            if(file.empty() || !std::filesystem::exists(file))
            {
                throw std::runtime_error("Required fmmTreeDepth input file is missing: " + file);
            }
            TextFileParser parser(file);
            if(!parser.readStringVector("masterFMM").empty()) {
    throw std::invalid_argument("Replace masterFMM with multipoleMode=none or multipoleMode=ewaldFMM");
  }
  for(const auto* key:{"fmmCoreWidth","fmmScreenWidth","fmmEwaldWidth","fmmShortCutoff","fmmBoxWidthFactor","ewaldWidthFactor","fmmSourcePointsPerLength","climbPreconditionerCutoff"}) {
    if(!parser.readStringVector(key).empty()) {
      throw std::invalid_argument(std::string("DD-specific setting ")+key+" does not belong in DefectiveCrystal.txt; use DD.txt: coreSize, ewaldWidthFactor, fmmSourcePointsPerLength, climbPreconditionerCutoff");
    }
  }
  const auto values=parser.readStringVector("multipoleMode");
            if(values.empty()) {
    return MultipoleMode::none;
  }
  std::string value,extra;std::istringstream stream(values.front().first);
            if(values.size()!=1||!(stream>>value)||(stream>>extra)) {
    throw std::invalid_argument("multipoleMode must be none, ewaldFMM, or FMM");
  }
  if(value=="none") {
    return MultipoleMode::none;
  }
  if(value=="ewaldFMM") {
    return MultipoleMode::ewaldFMM;
  }
  if(value=="FMM") {
    return MultipoleMode::FMM;
  }
  throw std::invalid_argument("Unknown multipoleMode: "+value+"; expected none, ewaldFMM, or FMM");
        }
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        struct Runtime
        {
            bool owned=false;
            Runtime(){if(Kokkos::is_finalized()) {
      throw std::logic_error("Kokkos was already finalized");
    }
    if(!Kokkos::is_initialized()){Kokkos::initialize();owned=true;}}
            ~Runtime(){if(owned&&Kokkos::is_initialized()&&!Kokkos::is_finalized()) {
      Kokkos::finalize();
    }
  }
        };
        std::shared_ptr<void> runtime()
        {static auto owner=std::make_shared<Runtime>();if(Kokkos::is_finalized()) {
    throw std::logic_error("Kokkos was finalized");
  }
  return owner;}
#endif

        template<int dim>
        unsigned readFmmTreeDepth(const std::string& file)
        {
            if(file.empty() || !std::filesystem::exists(file))
            {
                throw std::runtime_error("Required fmmTreeDepth input file is missing: " + file);
            }
            TextFileParser parser(file);
            const auto values = parser.readStringVector("fmmTreeDepth");
            if(values.empty())
            {
                throw std::runtime_error("Required fmmTreeDepth is missing in " + file);
            }
            long long depth = -1;
            std::string extra;
            std::istringstream stream(values.front().first);
            if(values.size() != 1 || !(stream >> depth) || (stream >> extra) ||
               depth < 0 || depth > BoxTree<dim>::maxDepth)
            {
                throw std::runtime_error("fmmTreeDepth in " + file +
                    " must be a single integer between 0 and " + std::to_string(BoxTree<dim>::maxDepth));
            }
            return static_cast<unsigned>(depth);
        }
    }


    template<int dim>
    FMMContext<dim>::FMMContext(const Vector& lower, const Vector& upper, const std::string& inputFile)
        : multipoleMode(readMode(inputFile)),
          interpolationOrder(usesFMM()?order(inputFile):4),
          fmmTreeDepth(readFmmTreeDepth<dim>(inputFile)),
          boxTree(std::make_shared<BoxTree<dim>>(lower, upper, fmmTreeDepth))
    {
        if(usesFMM())
        {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
            manager->keepRuntimeAlive(runtime());
#else
            throw std::runtime_error("FMM modes require a Kokkos-enabled build; set multipoleMode=none for analytical segments");
#endif
        }
        std::cout<<"multipoleMode="<<(usesEwaldFMM()?"ewaldFMM":(usesFMM()?"FMM":"none"))<<(usesFMM()?" (multipole stress and climb)":" (analytical segment stress and climb)")<<std::endl;
    }
    template struct FMMContext<1>;
    template struct FMMContext<2>;
    template struct FMMContext<3>;
}
