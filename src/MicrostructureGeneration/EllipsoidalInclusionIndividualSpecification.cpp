/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#include <EllipsoidalInclusionIndividualSpecification.h>

namespace model
{
    EllipsoidalInclusionIndividualSpecification::EllipsoidalInclusionIndividualSpecification() :
    MicrostructureSpecificationBase("EllipsoidalInclusion","Individual"),
    center(Eigen::Matrix<double,dim,1>::Zero()),
    crystalAxes_SI(Eigen::Matrix<double,dim,dim>::Zero()),
    crystalEigenDistortion(Eigen::Matrix<double,1,dim*dim>::Zero()),
    velocityReductionFactor(1.0), phaseID(0), allowOverlap(false), allowOutside(false)
    {
    }

    EllipsoidalInclusionIndividualSpecification::EllipsoidalInclusionIndividualSpecification(const std::string& fileName) :
    MicrostructureSpecificationBase("EllipsoidalInclusion","Individual",fileName),
    center(this->parser->readMatrix<double>("center",dim,1,true)),
    crystalAxes_SI(this->parser->readMatrix<double>("crystalAxes_SI",dim,dim,true)),
    crystalEigenDistortion(this->parser->readMatrix<double>("crystalEigenDistortion",1,dim*dim,true)),
    velocityReductionFactor(this->parser->readScalar<double>("velocityReductionFactor",true)),
    phaseID(this->parser->readScalar<int>("phaseID",true)),
    allowOverlap(this->parser->readScalar<int>("allowOverlap",true)),
    allowOutside(this->parser->readScalar<int>("allowOutside",true))
    {
    }
}
