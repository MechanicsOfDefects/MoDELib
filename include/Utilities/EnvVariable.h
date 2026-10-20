/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ENVVARIABLE_H_
#define model_ENVVARIABLE_H_


#include <string>

namespace model
{

    /************************************************************/
    /************************************************************/
    /*! \brief A simple class used to extract environment variables
     *
     * Sample usage:
     * \code
     * #include <iostream>
     * #include <EnvVariable.h>
     *
     * int main(){
     *  model::EnvVariable ev("PATH");
     *  std::cout<<"PATH is "<<ev<<std::endl;
     *  return 0;
     * }
     * \endcode
     */

    struct EnvVariable : public std::string
    {

        /******************************/
        static std::string getVariable(const std::string &varname);

        /******************************/
        EnvVariable(const std::string &varname);
    };

    /************************************************************/
    /************************************************************/
} // namespace model
#endif
