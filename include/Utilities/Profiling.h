/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_Profiling_H_
#define model_Profiling_H_
#include <cstdlib>
#include <string>
namespace model
{
    inline bool profilingEnabled()
    {
        static const bool enabled=[]()
        {
            const auto on=[](const char* key)
            {
                const char* value=std::getenv(key);
                return value && (std::string(value)=="1" || std::string(value)=="true" || std::string(value)=="on");
            };
            // Retain compatibility with the existing detailed FMM profiling gate.
            return on("MODELIB_PROFILE") || on("MODELIB_PROFILE_FMM");
        }();
        return enabled;
    }
}
#endif
