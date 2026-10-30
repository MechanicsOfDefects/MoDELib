#include <FMM/FMMContext.h>
#include <TextFileParser.h>
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
    FMMContext<dim>::FMMContext(const Vector& lower, const Vector& upper, const std::string& inputFile)
        : fmmTreeDepth(readFmmTreeDepth<dim>(inputFile)),
          boxTree(std::make_shared<BoxTree<dim>>(lower, upper, fmmTreeDepth))
    {
    }
    template struct FMMContext<1>;
    template struct FMMContext<2>;
    template struct FMMContext<3>;
}
