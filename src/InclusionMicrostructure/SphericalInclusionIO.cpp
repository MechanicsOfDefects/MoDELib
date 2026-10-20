/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <SphericalInclusionIO.h>
#include <SphericalInclusion.h>

namespace model
{

    template <int dim>
    SphericalInclusionIO<dim>::SphericalInclusionIO(const SphericalInclusion<dim> &ei)
        : inclusionID(ei.sID), C(ei.C), a(ei.a), eT(ei.eT), mobilityReduction(ei.mobilityReduction), phaseID(ei.phaseID)
    {
    }

    template <int dim>
    SphericalInclusionIO<dim>::SphericalInclusionIO(const size_t &inclusionID_in, const VectorDim &C_in,
                                                    const double &a_in, const MatrixDim &eT_in,
                                                    const double &mobilityReduction_in, const int &phaseID_in)
        : inclusionID(inclusionID_in), C(C_in), a(a_in), eT(eT_in), mobilityReduction(mobilityReduction_in),
          phaseID(phaseID_in)
    {
    }

    template <int dim>
    SphericalInclusionIO<dim>::SphericalInclusionIO()
        : inclusionID(0), C(VectorDim::Zero()), a(0), eT(MatrixDim::Zero()), mobilityReduction(1.0), phaseID(0)
    {
    }

    template <int dim>
    SphericalInclusionIO<dim>::SphericalInclusionIO(std::stringstream &ss)
        : inclusionID(0), C(VectorDim::Zero()), a(0), eT(MatrixDim::Zero()), mobilityReduction(1.0), phaseID(0)
    {
        ss >> inclusionID;
        for (int d = 0; d < dim; ++d)
        {
            ss >> C(d);
        }
        ss >> a;
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                double temp;
                ss >> temp;
                eT(i, j) = temp;
            }
        }
        ss >> mobilityReduction;
        ss >> phaseID;
    }

    template <int dim> void SphericalInclusionIO<dim>::write(std::ostream &os) const
    {
        os << inclusionID
           << "\t"
           /**/
           << std::setprecision(15) << std::scientific << C.transpose()
           << "\t"
           /**/
           << std::setprecision(15) << std::scientific << a << "\t";
        for (int d = 0; d < dim; ++d)
        {
            os << std::setprecision(15) << std::scientific << eT.row(d) << "\t";
        }
        os << std::setprecision(15) << std::scientific << mobilityReduction
           << "\t"
           /**/
           << phaseID;
        return;
    }
    template struct SphericalInclusionIO<3>;

} // namespace model
