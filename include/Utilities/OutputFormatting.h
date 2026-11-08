#ifndef model_OutputFormatting_H_
#define model_OutputFormatting_H_
#include <sstream>
#include <iomanip>
#include <string>
namespace model
{
    inline std::string formatSeconds(double seconds)
    {
        std::ostringstream out;
        out<<std::scientific<<std::setprecision(3)<<seconds;
        return out.str();
    }
}
#endif
