/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#include <InclusionMicrostructure/EllipsoidalInclusionIO.h>

namespace model
{
    template <int dim>
    EllipsoidalInclusionIO<dim>::EllipsoidalInclusionIO(
        const std::size_t& inclusionID_in, const VectorDim& C_in,
        const MatrixDim& axes_in, const MatrixDim& eT_in,
        const double& mobilityReduction_in, const int& phaseID_in) :
    inclusionID(inclusionID_in), C(C_in), axes(axes_in), eT(eT_in),
    mobilityReduction(mobilityReduction_in), phaseID(phaseID_in)
    {
    }

    template <int dim>
    EllipsoidalInclusionIO<dim>::EllipsoidalInclusionIO() :
    inclusionID(0), C(VectorDim::Zero()), axes(MatrixDim::Zero()),
    eT(MatrixDim::Zero()), mobilityReduction(1.0), phaseID(0)
    {
    }

    template <int dim>
    EllipsoidalInclusionIO<dim>::EllipsoidalInclusionIO(std::stringstream& ss) :
    EllipsoidalInclusionIO()
    {
        // Preserve stream failure flags for the caller, matching the other IO records.
        ss >> inclusionID;
        for(int d=0; d<dim; ++d)
        {
            ss >> C(d);
        }
        for(int i=0; i<dim; ++i)
        {
            for(int j=0; j<dim; ++j)
            {
                ss >> axes(i,j);
            }
        }
        for(int i=0; i<dim; ++i)
        {
            for(int j=0; j<dim; ++j)
            {
                ss >> eT(i,j);
            }
        }
        ss >> mobilityReduction >> phaseID;
    }

    template struct EllipsoidalInclusionIO<3>;
}
