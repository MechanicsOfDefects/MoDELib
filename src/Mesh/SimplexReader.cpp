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
            std::set<size_t> vertices;
            for (int d = 0; d < dim + 1; ++d)
            {
                vertices.insert(eIter.second[d]);
            }
            size_t regionID(eIter.second[dim + 1]);
            this->elements().emplace(eIter.first, std::make_pair(typename SimplexTraits<dim, dim>::SimplexIDType(vertices), regionID));
        }
        //            assert(elementReader.size()==this->elements().size());
    }

    template <int dim>
    void SimplexReader<dim>::read(const std::string &meshFileName, const Eigen::Matrix<double, dim, dim> &A,
                                  const Eigen::Matrix<double, dim, 1> &x0)
    {
        this->clear();
        GmshReader reader(meshFileName, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
        for (const auto &node : reader.nodes())
        {
            this->nodes().emplace(node.first, A * (node.second.template head<dim>() - x0));
        }
        constexpr int linearType = dim == 3 ? 4 : (dim == 2 ? 2 : 1);
        constexpr int quadraticType = dim == 3 ? 11 : (dim == 2 ? 9 : 8);
        for (const auto &entry : reader.elements())
        {
            const auto &element = entry.second;
            if (element.type == linearType || element.type == quadraticType)
            {
                if (element.nodeIDs.size() < dim + 1 || element.tags.size() < 2)
                {
                    throw std::runtime_error("Mesh element is missing vertices or region tag");
                }
                std::set<size_t> vertices;
                for (int k = 0; k <= dim; ++k)
                {
                    vertices.insert(element.nodeIDs[k]);
                }
                this->elements().emplace(entry.first, std::make_pair(
                    typename SimplexTraits<dim, dim>::SimplexIDType(vertices), element.tags[1]));
            }
        }
    }
    template struct SimplexReaderBase<1>;
    template struct SimplexReaderBase<2>;
    template struct SimplexReaderBase<3>;
    template struct SimplexReader<1>;
    template struct SimplexReader<2>;
    template struct SimplexReader<3>;
} // namespace model
