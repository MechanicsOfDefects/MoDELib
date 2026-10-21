/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalAuxIO_H_
#define model_DefectiveCrystalAuxIO_H_

#include <vector>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <cfloat>      // std::ifstream
#include <DefectiveCrystal/DefectiveCrystalBaseIO.h>
#include <PeriodicGlidePlane.h>
//#include <GlidePlaneIO.h>
#include <GlidePlaneModule.h>
#include <DislocationDynamics/DislocationQuadraturePointIO.h>
#include <FEM/MeshNodeIO.h>
#include <DefectiveCrystal/DefectiveCrystalParameters.h>

namespace model
{
        
    template <int dim>
    struct DefectiveCrystalAuxIO : public DefectiveCrystalBaseIO
    /*            */,private std::vector<MeshNodeIO<dim>>
//    /*            */,private std::vector<GlidePlaneIO<dim>>
    /*            */,private std::vector<PeriodicPlanePatchIO<dim>>
    /*            */,private std::vector<DislocationQuadraturePointIO<dim>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
//        DefectiveCrystalAuxIO(const std::string& folderName,const std::string& suffix="");
        DefectiveCrystalAuxIO(const std::string& folderName="config");

        void clear();
        void addPeriodicGlidePlane(const PeriodicGlidePlane<dim>& pgp);
        const std::vector<MeshNodeIO<dim>>& meshNodes() const;
        std::vector<MeshNodeIO<dim>>& meshNodes();
//        const std::vector<GlidePlaneIO<dim>>& glidePlanes() const;
//        std::vector<GlidePlaneIO<dim>>& glidePlanes();
        const std::vector<PeriodicPlanePatchIO<dim>>& periodicGlidePlanePatches() const;
        std::vector<PeriodicPlanePatchIO<dim>>& periodicGlidePlanePatches();
        const std::vector<DislocationQuadraturePointIO<dim>>& quadraturePoints() const;
        std::vector<DislocationQuadraturePointIO<dim>>& quadraturePoints();
        void write(const size_t& runID,const bool& outputBinary);
        void writeTxt(const long int& runID);
        void writeBin(const size_t& runID);
        void readBin(const size_t& runID);
        void readTxt(const size_t& runID);
        void printLog(const std::chrono::time_point<std::chrono::system_clock>& t0) const;
        void read(const size_t& runID);
        
    };
    
}
#endif

