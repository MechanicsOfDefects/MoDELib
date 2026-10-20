/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalBaseIO_cpp_
#define model_DefectiveCrystalBaseIO_cpp_


#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <string>      // std::ifstream
#include <DefectiveCrystal/DefectiveCrystalBaseIO.h>

namespace model
{
    
        
        DefectiveCrystalBaseIO::DefectiveCrystalBaseIO(const std::string& outFolderName_in,
                 const std::string& fileNamePrefix_in,
                 const std::string& outFolderSuffix_in) :
        /* init */ outFolderName(outFolderName_in)
        /* init */,fileNamePrefix(fileNamePrefix_in)
        /* init */,outFolderSuffix(outFolderSuffix_in)
        {
        }
        

        std::string DefectiveCrystalBaseIO::getBinFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+fileNamePrefix+"_"+std::to_string(runID)+".bin";
        }
        
        std::string DefectiveCrystalBaseIO::getTxtFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+fileNamePrefix+"_"+std::to_string(runID)+".txt";
        }
        
        std::string DefectiveCrystalBaseIO::getTxtSegmentFilename(const size_t& runID)
        {
            return outFolderName+outFolderSuffix+"/"+"segments"+std::to_string(runID)+".txt";
        }
        
        bool DefectiveCrystalBaseIO::isBinGood(const size_t& frameID)
        {
            return std::ifstream(getBinFilename(frameID).c_str(), std::ios::in|std::ios::binary).good();

        }

        bool DefectiveCrystalBaseIO::isTxtGood(const size_t& frameID)
        {
            return std::ifstream(getTxtFilename(frameID).c_str(), std::ios::in).good();
            
        }

    
}
#endif

