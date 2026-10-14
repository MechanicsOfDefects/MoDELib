/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_EllipsoidalInclusionIndividualSpecification_H_
#define model_EllipsoidalInclusionIndividualSpecification_H_

#include <Eigen/Dense>
#include <MicrostructureSpecificationBase.h>

namespace model
{
    struct EllipsoidalInclusionIndividualSpecification : public MicrostructureSpecificationBase
    {
        static constexpr int dim=3;

        Eigen::Matrix<double,dim,1> center; // Simulation length units (Burgers vector).
        Eigen::Matrix<double,dim,dim> axes_SI; // Columns are orthogonal semiaxis vectors [m].
        // Dimensionless components 11,12,13,21,22,23,31,32,33, as in other specifications.
        Eigen::Matrix<double,1,dim*dim> eigenDistortion;
        double velocityReductionFactor;
        int phaseID;
        bool allowOverlap;
        bool allowOutside;

        // Zero geometry must be filled in before generating the inclusion.
        EllipsoidalInclusionIndividualSpecification();
        EllipsoidalInclusionIndividualSpecification(const std::string& fileName);
    };
}
#endif
