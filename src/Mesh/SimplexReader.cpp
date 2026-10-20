/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <SimplexReader.h>
#include <IDreader.h>
#include <GmshReader.h>

namespace model
{

    template <int dim> auto SimplexReaderBase<dim>::nodes() const -> const NodeContainerType &
    {
        return *this;
    }

    template <int dim> auto SimplexReaderBase<dim>::nodes() -> NodeContainerType &
    {
        return *this;
    }

    template <int dim> auto SimplexReaderBase<dim>::elements() const -> const ElementContainerType &
    {
        return *this;
    }

    template <int dim> auto SimplexReaderBase<dim>::elements() -> ElementContainerType &
    {
        return *this;
    }

    template <int dim> auto SimplexReaderBase<dim>::clear() -> void
    {
        nodes().clear();
        elements().clear();
    }

    template <int dim> auto SimplexReaderBase<dim>::simplexReader() -> SimplexReaderBase<dim> &
    {
        return *this;
    }

    template <int dim> auto SimplexReaderBase<dim>::simplexReader() const -> const SimplexReaderBase<dim> &
    {
        return *this;
    }

    template <int dim>
    auto SimplexReaderBase<dim>::get_P0(const typename SimplexTraits<dim, 0>::SimplexIDType &xID)
        -> const Eigen::Matrix<double, dim, 1> &
    {
        //            const typename VertexReader<'N',dim+1,double>::const_iterator nIter(nodeReader.find((xID)(0)));
        const auto nIter(nodes().find((xID)(0)));
        assert((nIter != nodes().end()) && "MESH VERTEX NOT FOUND IN N/N_x.txt.");
        return nIter->second;
    }

    template <int dim> auto SimplexReader<dim>::read(const int &meshID) -> void
    {

        this->clear();
        IDreader<1, dim, double> nodeReader("N/N");
        nodeReader.read(meshID, true);
        for (const auto &node : nodeReader)
        {
            this->nodes().emplace(node.first, Eigen::Map<const Eigen::Matrix<double, dim, 1>>(node.second.data()));
        }
        assert(nodeReader.size() == this->nodes().size());

        IDreader<1, dim + 2, size_t> elementReader("T/T");
        elementReader.read(meshID, true);
        for (const auto &eIter : elementReader)
        {
            typename SimplexTraits<dim, dim>::SimplexIDType key;
            for (int d = 0; d < dim + 1; ++d)
            {
                key[d] = eIter.second[d];
            }
            size_t regionID(eIter.second[dim + 1]);
            this->elements().emplace(eIter.first, std::make_pair(key, regionID));
        }
        //            assert(elementReader.size()==this->elements().size());
    }

    auto SimplexReader<3>::read(const std::string &meshFileName, const Eigen::Matrix<double, 3, 3> &A,
                                const Eigen::Matrix<double, 3, 1> &x0) -> void
    {
        this->clear();
        GmshReader gmshReader(meshFileName, A, x0);
        this->nodes() = gmshReader.nodes();
        // assert(nodeReader.size()==this->nodes().size());

        for (const auto &eIter : gmshReader.elements())
        {
            if (eIter.second.type == 4     // 4-node (linear) tetrahedron
                || eIter.second.type == 11 // 11-node (quadratic) tetrahedron
            )
            {
                //                    typename SimplexTraits<dim,dim>::SimplexIDType key;
                std::set<size_t> set;
                for (int d = 0; d < dim + 1; ++d)
                {
                    set.insert(eIter.second.nodeIDs[d]);
                    //                        key[d]=eIter.second.nodeIDs[d]; // first dim+1 nodes are the principa
                    //                        vertex nodes. See Node order in
                    //                        http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
                }
                size_t regionID(eIter.second.tags[1]); // FINISH HERE, WHAT IS THE RIGHT TAG FOR A POLYCRYSTAL
                this->elements().emplace(
                    eIter.first, std::make_pair(typename SimplexTraits<dim, dim>::SimplexIDType(set), regionID));
            }
        }
    }
    template struct SimplexReaderBase<3>;

} // namespace model
