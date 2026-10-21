/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <EnvVariable.h>
#include <cstdlib>
#include <string.h>
#include <assert.h>

namespace model
{

    auto EnvVariable::getVariable(const std::string &varname) -> std::string
    {
        const char *tmpvar = std::getenv(varname.c_str());
        assert(tmpvar != 0 && *tmpvar != '\0' && "FAILED TO GET ENVIRONMENT VARIABLE.");
        return tmpvar;
    }

    EnvVariable::EnvVariable(const std::string &varname) : std::string(getVariable(varname))
    {
    }
} // namespace model
