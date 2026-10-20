/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_EllipsoidalInclusionIO_H_
#define model_EllipsoidalInclusionIO_H_

#include <cstddef>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

namespace model
{
    template <int dim>
    struct EllipsoidalInclusionIO
    {
        static_assert(dim == 3, "EllipsoidalInclusionIO requires dim=3");
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;

        std::size_t inclusionID;
        VectorDim C;
        MatrixDim axes; // Columns are mutually orthogonal semiaxis vectors.
        MatrixDim eT;
        double mobilityReduction;
        int phaseID;

        template <typename EllipsoidalInclusionType>
        EllipsoidalInclusionIO(const EllipsoidalInclusionType& ei) :
        inclusionID(ei.sID), C(ei.C), axes(ei.axes), eT(ei.eT),
        mobilityReduction(ei.mobilityReduction), phaseID(ei.phaseID)
        {
        }

        EllipsoidalInclusionIO(const std::size_t& inclusionID_in,
                              const VectorDim& C_in,
                              const MatrixDim& axes_in,
                              const MatrixDim& eT_in,
                              const double& mobilityReduction_in,
                              const int& phaseID_in);
        EllipsoidalInclusionIO();
        EllipsoidalInclusionIO(std::stringstream& ss);

        // One record: ID, C, axes (row-major), eT (row-major), mobility, phase.
        // Material constants and numerical settings are supplied when reconstructing
        // the inclusion, as material constants are for SphericalInclusionIO.
        template <class T>
        friend T& operator << (T& os, const EllipsoidalInclusionIO<dim>& ei)
        {
            os << ei.inclusionID << "\t"
               << std::setprecision(15) << std::scientific << ei.C.transpose() << "\t";
            for(int d=0; d<dim; ++d)
            {
                os << ei.axes.row(d) << "\t";
            }
            for(int d=0; d<dim; ++d)
            {
                os << ei.eT.row(d) << "\t";
            }
            os << ei.mobilityReduction << "\t" << ei.phaseID;
            return os;
        }
    };
}
#endif
