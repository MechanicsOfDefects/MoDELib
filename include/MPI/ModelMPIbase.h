/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _MODELMPIBASE_h_
#define _MODELMPIBASE_h_

#include <string>

namespace model
{

    class ModelMPIbase
    {

        static int _mpiInitialized;
        static int _mpiRank;
        static int _mpiProcs;
        static std::string _processorName;

      public:
        /**********************************************************************/
        ModelMPIbase(int argc, char *argv[]);

        /**********************************************************************/
        ModelMPIbase();

        /**********************************************************************/
        ~ModelMPIbase();

        /**********************************************************************/
        static void initMPI(int argc, char *argv[]);

        static const int &mpiRank();
        static const int &mpiProcs();
        static const int &mpiInitialized();
        static const std::string &processorName();
    };

} // namespace model
#endif
