/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityPyComponent_h_
#define _model_DislocationMobilityPyComponent_h_

#include <DislocationMobilityBase.h>

#ifdef _MODEL_PYBIND11_
    #undef slots
    #include <pybind11/embed.h>
    #include <pybind11/eigen.h>
    #include <pybind11/stl.h>
    #define slots Q_SLOTS
#endif

namespace model
{

    /*! A mobility component (e.g. the edge or screw part of an edgeScrew mobility)
     *  evaluated by a Python module. Selected in the material file with
     *      mobility_<suffix>=python <path/to/module.py> <soluteContent>;
     *  The module must define a class MobilitySolver with a method
     *      velocityPy(stress, burgers, normal, temp, dL, dt, soluteContent)
     *  where stress is in GPa and temp in K. velocityPy must return the
     *  magnitude of the dislocation velocity in m/s.
     *  A relative module path is relative to the folder of the material file.
     */
    struct DislocationMobilityPyComponent : public DislocationMobilityBase
    {

        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;

        static constexpr double kB_SI=1.38064852e-23; // [J/K]
        const double kB;
        const double mu_SI;
        const double Tm;
        const double cs;
        const std::string pyModuleName;
        const double soluteContent;

#ifdef _MODEL_PYBIND11_
        pybind11::module pyModule;
        pybind11::object mobilitySolver;
#endif

        DislocationMobilityPyComponent(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in,const double& soluteContent_in);

        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& sfg) const override;

    };

}
#endif
