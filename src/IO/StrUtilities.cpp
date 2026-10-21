/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <StrUtilities.h>
#include <algorithm>
#include <cctype>

namespace model
{

    auto StrUtilities::removeSpaces(std::string key) -> std::string
    {
        key.erase(std::remove_if(key.begin(), key.end(), [](unsigned char x) { return std::isspace(x); }), key.end());
        return key;
    }

    auto StrUtilities::lowercase(std::string text) -> std::string
    {
        std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) { return std::tolower(c); });
        return text;
    }
} // namespace model
