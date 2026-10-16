/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * Written in 2024 by Matthew Maron <mlm335@miami.edu>
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityPyComponent_cpp_
#define _model_DislocationMobilityPyComponent_cpp_

#include <DislocationMobilityPyComponent.h>
#include <filesystem>

namespace model
{
#ifdef _MODEL_PYBIND11_ // COMPILED WITH PYBIND11
    DislocationMobilityPyComponent::DislocationMobilityPyComponent(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in,const double& soluteContent_in) :
    /* init */ DislocationMobilityBase("DislocationMobilityPyComponent for "+material.materialName)
    /* init */,kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,mu_SI(material.mu_SI)
    /* init */,Tm(material.Tm)
    /* init */,cs(material.cs_SI)
    /* init */,pyModuleName(pyModuleName_in)
    /* init */,soluteContent(soluteContent_in)
    {// Set up pyModule
        std::filesystem::path modulePath(pyModuleName);
        if(modulePath.is_relative())
        {// a relative path is relative to the folder of the material file (the inputFiles folder)
            modulePath=(std::filesystem::path(material.materialFile).parent_path()/modulePath).lexically_normal();
        }
        const std::string moduleStem(modulePath.stem().string());
        const std::string moduleDir(modulePath.parent_path().string());
        std::cout<<"moduleStem="<<moduleStem<<std::endl;
        std::cout<<"moduleDir="<<moduleDir<<std::endl;
        pybind11::module sys = pybind11::module::import("sys");
        pybind11::list path = sys.attr("path");
        path.append(moduleDir);
        pyModule = pybind11::module::import(moduleStem.c_str());
        mobilitySolver = pyModule.attr("MobilitySolver")(); // constructed once, reused for every velocity evaluation
    }

    double DislocationMobilityPyComponent::velocity(const MatrixDim& S,
                                                    const VectorDim& b,
                                                    const VectorDim& n,
                                                    const double& T,
                                                    const double& dL,
                                                    const double& dt,
                                                    const std::shared_ptr<StochasticForceGenerator>& ) const
    {
        pybind11::gil_scoped_acquire gil; // velocity may be called from OpenMP threads
        try
        {
            const Eigen::MatrixXd stress(S*mu_SI*1e-9); // GPa
            const Eigen::VectorXd burgers(b);
            const Eigen::VectorXd normal(n);
            const double velocity_SI = mobilitySolver.attr("velocityPy")(stress, burgers, normal, T, dL, dt, soluteContent).cast<double>(); // m/s
            return velocity_SI/cs; // m/s to code units
        }
        catch (const pybind11::error_already_set& e)
        {
            std::cerr << "Python error: " << e.what() << std::endl;
            std::terminate();
        }
    }

#else // COMPILED WITHOUT PYBIND11
    DislocationMobilityPyComponent::DislocationMobilityPyComponent(const PolycrystallineMaterialBase& material,const std::string& pyModuleName_in,const double& soluteContent_in) :
    /* init */ DislocationMobilityBase("DislocationMobilityPyComponent for "+material.materialName)
    /* init */,kB(kB_SI/material.mu_SI/std::pow(material.b_SI,3))
    /* init */,mu_SI(material.mu_SI)
    /* init */,Tm(material.Tm)
    /* init */,cs(material.cs_SI)
    /* init */,pyModuleName(pyModuleName_in)
    /* init */,soluteContent(soluteContent_in)
    {
    }

    double DislocationMobilityPyComponent::velocity(const MatrixDim& ,
                                                    const VectorDim& ,
                                                    const VectorDim& ,
                                                    const double& ,
                                                    const double& ,
                                                    const double& ,
                                                    const std::shared_ptr<StochasticForceGenerator>& ) const
    {
        throw std::runtime_error("DislocationMobilityPyComponent used without pybind11");
        return 0.0;
    }
#endif

}
#endif
