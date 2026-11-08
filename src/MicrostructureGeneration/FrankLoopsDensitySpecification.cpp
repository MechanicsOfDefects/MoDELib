/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FrankLoopsDensitySpecification_cpp_
#define model_FrankLoopsDensitySpecification_cpp_

#include <FrankLoopsDensitySpecification.h>
#include <cmath>
#include <numbers>
#include <limits>
#include <algorithm>
#include <stdexcept>

namespace model
{
    FrankLoopsDensitySpecification::FrankLoopsDensitySpecification():
    /* init */ MicrostructureSpecificationBase("FrankLoops","Density")
    /* init */,targetDensity(0.0)
    /* init */,numberOfSides(0)
    /* init */,targetSegmentLength(0.0)
    /* init */,radiusDistributionMean(0.0)
    /* init */,radiusDistributionStd(0.0)
    /* init */,areVacancyLoops(true)
    /* init */,allowedGrainIDs({-1})
    /* init */,allowedPlaneIDs({-1})
    {
        
    }

    FrankLoopsDensitySpecification::FrankLoopsDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("FrankLoops","Density",fileName)
    /* init */,targetDensity(this->parser->readScalar<double>("targetDensity_SI",true))
    /* init */,numberOfSides(targetDensity>0.0 && !this->parser->readStringVector("numberOfSides").empty()? this->parser->readScalar<int>("numberOfSides",true) : 0)
    /* init */,targetSegmentLength(this->parser->readStringVector("targetSegmentLength").empty()? 0.0 : this->parser->readScalar<double>("targetSegmentLength",true))
    /* init */,radiusDistributionMean(targetDensity>0.0? this->parser->readScalar<double>("radiusDistributionMean_SI",true) : 0.0)
    /* init */,radiusDistributionStd(targetDensity>0.0? this->parser->readScalar<double>("radiusDistributionStd_SI",true) : 0.0)
    /* init */,areVacancyLoops(targetDensity>0.0? this->parser->readScalar<int>("areVacancyLoops",true) : 1)
    /* init */,allowedGrainIDs(this->parser->readArray<int>("allowedGrainIDs",true))
    /* init */,allowedPlaneIDs(this->parser->readArray<int>("allowedPlaneIDs",true))
    {
        
    }
    void FrankLoopsDensitySpecification::validateDiscretization() const
    {
        if(!std::isfinite(targetSegmentLength) || targetSegmentLength<0.0 ||
           (targetSegmentLength==0.0 && numberOfSides<3))
        {
            throw std::invalid_argument("Frank loops require a positive targetSegmentLength or numberOfSides >= 3");
        }
    }
    size_t FrankLoopsDensitySpecification::sidesForRadius(double radius) const
    {
        validateDiscretization();
        if(!std::isfinite(radius) || radius<=0.0)
        {
            throw std::invalid_argument("Frank loop radius must be finite and positive");
        }
        if(targetSegmentLength==0.0) { return static_cast<size_t>(numberOfSides); }
        const double count=std::max(8.0,std::ceil(2.0*std::numbers::pi*radius/targetSegmentLength));
        if(!std::isfinite(count) || count>std::numeric_limits<int>::max())
        {
            throw std::overflow_error("Frank loop discretization exceeds supported node count");
        }
        return static_cast<size_t>(count);
    }
}
#endif
