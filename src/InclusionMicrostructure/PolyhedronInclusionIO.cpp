/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <PolyhedronInclusionIO.h>
#include <PolyhedronInclusion.h>

namespace model
{

    template <short unsigned int dim>
    PolyhedronInclusionIO<dim>::PolyhedronInclusionIO(const PolyhedronInclusion<dim> &ei)
        : inclusionID(ei.sID), eT(ei.eT), mobilityReduction(ei.mobilityReduction), phaseID(ei.phaseID)
    {
    }

    template <short unsigned int dim>
    PolyhedronInclusionIO<dim>::PolyhedronInclusionIO(const size_t &inclusionID_in, const MatrixDim &eT_in,
                                                      const double &mobilityReduction_in, const int &phaseID_in)
        : inclusionID(inclusionID_in), eT(eT_in), mobilityReduction(mobilityReduction_in), phaseID(phaseID_in)
    {
    }

    template <short unsigned int dim>
    PolyhedronInclusionIO<dim>::PolyhedronInclusionIO()
        : inclusionID(0), eT(MatrixDim::Zero()), mobilityReduction(1.0), phaseID(0)
    {
    }

    template <short unsigned int dim>
    PolyhedronInclusionIO<dim>::PolyhedronInclusionIO(std::stringstream &ss)
        : inclusionID(0), eT(MatrixDim::Zero()), mobilityReduction(1.0), phaseID(0)
    {
        ss >> inclusionID;
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

    template <short unsigned int dim> void PolyhedronInclusionIO<dim>::write(std::ostream &os) const
    {
        os << inclusionID << "\t";
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
    template struct PolyhedronInclusionIO<3>;

} // namespace model
