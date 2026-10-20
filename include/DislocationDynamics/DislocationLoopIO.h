/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopIO_H_
#define model_DislocationLoopIO_H_

#include <iomanip>
#include <tuple>
#include <sstream>
#include <ostream>
#include <Eigen/Dense>
#include <Eigen/Dense>

namespace model
{
    template<int dim> class DislocationLoop;
    
    template<short unsigned int dim>
    struct DislocationLoopIO
    {
        
        enum DislocationLoopType{GLISSILELOOP,SESSILELOOP,VIRTUALLOOP};

        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
         size_t sID;          // sID
         VectorDim B;          // position
         VectorDim N;          // normal
         VectorDim P;          // position
         size_t grainID;          // component ID
         int loopType;
         std::tuple<double,double,double,double> loopLength;
         double slippedArea;
        
        DislocationLoopIO(const DislocationLoop<dim>& dL);
        
        DislocationLoopIO(const size_t& sID_in,         // sID
                          const VectorDim& B_in,          // position
                          const VectorDim& N_in,          // velocity
                          const VectorDim& P_in,          // velocity
                          const size_t& grainID_in,
                          const int& loopType_in
                          );
        
        DislocationLoopIO();
        
        DislocationLoopIO(std::stringstream& ss);
        
        template <class T>
        friend T& operator << (T& os, const DislocationLoopIO<dim>& ds)
        {
            ds.write(os);
            return os;
        }

        std::ostream& write(std::ostream& os) const;
        
	};
	
}
#endif

