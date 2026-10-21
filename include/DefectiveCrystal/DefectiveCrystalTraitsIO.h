/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalTraitsIO_H_
#define model_DefectiveCrystalTraitsIO_H_

#include <string>

namespace model
{
    
    
    struct DefectiveCrystalTraitsIO
    {
        
        enum SimulationType{FINITE_NO_FEM=0,FINITE_FEM=1,PERIODIC_IMAGES=2,PERIODIC_FEM=3};

        
        const std::string simulationFolder;
        const std::string inputFilesFolder;
        const std::string configFolder;
        const std::string auxFolder;
        const std::string fFolder;
        const std::string dcFile;
        const std::string fFile;
        const std::string flabFile;
        const std::string polyFile;
        const std::string materialFile;
        const std::string microstructureFile;
        const std::string meshFile;

        DefectiveCrystalTraitsIO(const std::string& folderName);
        
        
    };
    
}
#endif

