/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StrUtilities_H_
#define model_StrUtilities_H_

// #include <filesystem>
#include <string>

namespace model
{

    struct StrUtilities
    {

        static std::string removeSpaces(std::string key);

        static std::string lowercase(std::string text);
    };

} // namespace model
#endif
