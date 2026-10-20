/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <PolyhedronInclusionNodeIO.h>

namespace model
{

    template <short unsigned int dim>
    PolyhedronInclusionNodeIO<dim>::PolyhedronInclusionNodeIO(

        const size_t &nodeID_in, const VectorDim &P_in)
        :

          nodeID(nodeID_in), P(P_in)
    {

        //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
    }

    template <short unsigned int dim>
    PolyhedronInclusionNodeIO<dim>::PolyhedronInclusionNodeIO()
        :

          nodeID(0), P(VectorDim::Zero())
    {
    }

    template <short unsigned int dim>
    PolyhedronInclusionNodeIO<dim>::PolyhedronInclusionNodeIO(std::stringstream &ss)
        :

          nodeID(0), P(VectorDim::Zero())
    {
        //            ss>>inclusionID;
        ss >> nodeID;
        for (int d = 0; d < dim; ++d)
        {
            ss >> P(d);
        }
    }

    template <short unsigned int dim> void PolyhedronInclusionNodeIO<dim>::write(std::ostream &os) const
    {
        os //<< inclusionID<<"\t"
            << nodeID
            << "\t"
            /**/
            << P.transpose();
        return;
    }
    template struct PolyhedronInclusionNodeIO<3>;

} // namespace model
