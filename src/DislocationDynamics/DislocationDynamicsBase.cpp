#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <CTM.h>
#include <filesystem>
#include <sstream>
#include <stdexcept>
namespace model
{
    namespace
    {
        template<int dim>
        unsigned readFmmTreeDepth(const std::string& file)
        {
            if(file.empty() || !std::filesystem::exists(file))
            {
                throw std::runtime_error("Required fmmTreeDepth input file is missing: " + file);
            }
            TextFileParser parser(file);
            const auto values = parser.readStringVector("fmmTreeDepth");
            if(values.empty())
            {
                throw std::runtime_error("Required fmmTreeDepth is missing in " + file);
            }
            long long depth = -1;
            std::string extra;
            std::istringstream stream(values.front().first);
            if(values.size() != 1 || !(stream >> depth) || (stream >> extra) ||
               depth < 0 || depth > BoxTree<dim>::maxDepth)
            {
                throw std::runtime_error("fmmTreeDepth in " + file +
                    " must be a single integer between 0 and " + std::to_string(BoxTree<dim>::maxDepth));
            }
            return static_cast<unsigned>(depth);
        }
    }

    template<int dim>
    DislocationDynamicsBase<dim>::DislocationDynamicsBase(
        std::shared_ptr<const SimplicialMesh<dim>> mesh,
        std::shared_ptr<const Polycrystal<dim>> poly,
        const Eigen::Matrix<double,dim,Eigen::Dynamic>& basis,
        const std::string& inputFile,
        const std::string& ddInputFile) :
        meshOwner(std::move(mesh)),
        polyOwner(std::move(poly)),
        glidePlaneFactory(*polyOwner),
        periodicGlidePlaneFactory(*polyOwner,glidePlaneFactory),
        EwaldLength(basis.cols()>0 ? TextFileParser(inputFile).readScalar<double>("EwaldLengthFactor",true)
            *std::pow(std::sqrt((basis.transpose()*basis).determinant())/CTM::factorial(basis.cols()),1.0/basis.cols()) : 0.0),
        fmmTreeDepth(readFmmTreeDepth<dim>(ddInputFile)),
        boxTree(std::make_shared<BoxTree<dim>>(meshOwner->xMin(), meshOwner->xMax(), fmmTreeDepth))
    {
    }
    template struct DislocationDynamicsBase<3>;
}
