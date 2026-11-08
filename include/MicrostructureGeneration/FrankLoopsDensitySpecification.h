/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FrankLoopsDensitySpecification_H_
#define model_FrankLoopsDensitySpecification_H_


#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct FrankLoopsDensitySpecification : public MicrostructureSpecificationBase
    {
        double targetDensity;
        int numberOfSides;
        // Positive: choose max(8,ceil(2*pi*radius/length)); zero: fixed numberOfSides.
        // Both radius and length are in Burgers-vector units.
        double targetSegmentLength;
        void validateDiscretization() const;
        size_t sidesForRadius(double radius) const;
        double radiusDistributionMean;
        double radiusDistributionStd;
        bool areVacancyLoops;
        std::vector<int> allowedGrainIDs;
        std::vector<int> allowedPlaneIDs;

        
        FrankLoopsDensitySpecification();
        FrankLoopsDensitySpecification(const std::string& fileName);
    };
}
#endif
