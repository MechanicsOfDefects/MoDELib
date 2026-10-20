/* This file is part of MoDELib, distributed under the GNU GPL v2. */
#ifndef model_FEMLoadController_H_
#define model_FEMLoadController_H_
#include <BilinearWeakExpression.h>
#include <Quadrature.h>
#include <BarycentricTraits.h>
#include <Eigen/QR>
#include <Eigen/SparseCore>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <vector>
namespace model
{
    // Mixed loading acts on volume-average engineering strain. Its variational
    // force supplies uniform traction on physical faces and work on free periodic jumps.
    template<class Trial>
    class FEMLoadController : public BilinearWeakExpression<FEMLoadController<Trial>>
    {
        std::shared_ptr<Trial> trial;
        Eigen::Matrix<double,6,1> alpha,D,strain0,strainRate,stress0,stressRate;
        std::vector<int> fixed;
        Eigen::MatrixXd meanStrain,periodicDirections;
        double volume=0.0;
        double length=1.0;
    public:
        using TrialFunctionType=Trial;
        using Vector6=Eigen::Matrix<double,6,1>;
        // Values such as 1e20 conventionally denote strain control. Enforce
        // these exactly rather than placing enormous coefficients in the matrix.
        static constexpr double rigidDeviceThreshold=1.e15;
        struct Load
        {
            std::shared_ptr<Trial> trial;
            Eigen::VectorXd force;
            const auto& field() const
            {
                return trial->field();
            }
            Eigen::VectorXd globalVector() const
            {
                return force;
            }
        };
        template<class Settings>
        FEMLoadController(std::shared_ptr<Trial> u,const Settings& settings):
            trial(std::move(u)),alpha(settings.alpha),D(Vector6::Zero()),
            strain0(settings.g0),strainRate(settings.g0Dot),stress0(settings.f0),stressRate(settings.f0Dot)
        {
            if(!alpha.allFinite() || (alpha.array()<0.0).any() || !strain0.allFinite() ||
               !strainRate.allFinite() || !stress0.allFinite() || !stressRate.allFinite())
            {
                throw std::invalid_argument("Invalid FEM loading inputs; stiffnessRatio must be finite and nonnegative.");
            }
            const auto& fe=trial->fe();
            volume=fe.mesh.volume();
            length=(fe.mesh.xMax()-fe.mesh.xMin()).norm();
            meanStrain=Eigen::MatrixXd::Zero(6,trial->gSize());
            using Q=Quadrature<3,4,GaussLegendre>;
            for(const auto& entry : fe.elements())
            {
                const auto& e=entry.second;
                for(size_t k=0;k<Q::quadratureOrder;++k)
                {
                    const auto bary=BarycentricTraits<3>::x2l(Q::abscissa(k));
                    const auto local=(trial->sfmDef(e,bary)*(e.absJ(bary)*Q::weight(k)/volume)).eval();
                    for(int n=0;n<Trial::nodesPerElement;++n)
                    {
                        meanStrain.middleCols(3*e.node(n).gID,3)+=local.middleCols(3*n,3);
                    }
                }
            }
            for(int i=0;i<6;++i)
            {
                if(alpha(i)>=rigidDeviceThreshold)
                {
                    fixed.push_back(i);
                }
                else
                {
                    D(i)=alpha(i)*settings.A(i,i);
                }
            }
            const auto& basis=fe.mesh.periodicBasis();
            periodicDirections=Eigen::MatrixXd::Zero(3,basis.cols());
            if(basis.cols()>0)
            {
                periodicDirections=basis.householderQr().householderQ()*Eigen::MatrixXd::Identity(3,basis.cols());
            }
        }
        const auto& field() const
        {
            return trial->field();
        }
        std::vector<Eigen::Triplet<double>> globalTriplets() const
        {
            return {};
        }
        Eigen::MatrixXd lowRankPenaltyRows() const
        {
            return (volume*D.array()).sqrt().matrix().asDiagonal()*meanStrain;
        }
        Eigen::MatrixXd constraintRows() const
        {
            Eigen::MatrixXd rows(fixed.size(),trial->gSize());
            for(size_t i=0;i<fixed.size();++i)
            {
                rows.row(i)=std::sqrt(volume)*meanStrain.row(fixed[i]);
            }
            return rows;
        }
        Eigen::VectorXd constraintValues(double time,const Vector6& strainOffset) const
        {
            const Vector6 target=strain0+time*strainRate-strainOffset;
            Eigen::VectorXd values(fixed.size());
            for(size_t i=0;i<fixed.size();++i)
            {
                values(i)=std::sqrt(volume)*target(fixed[i]);
            }
            return values;
        }
        std::shared_ptr<Load> rightHandSide(const Eigen::VectorXd& source,double time,const Vector6& strainOffset) const
        {
            const Vector6 force=stress0+time*stressRate+D.asDiagonal()*(strain0+time*strainRate-strainOffset);
            return std::make_shared<Load>(Load{trial,source+volume*meanStrain.transpose()*force});
        }
        Eigen::MatrixXd affineModes() const
        {
            Eigen::MatrixXd modes=Eigen::MatrixXd::Zero(trial->gSize(),3*periodicDirections.cols());
            for(size_t node=0;node<trial->fe().nodeSize();++node)
            {
                const Eigen::VectorXd shift=periodicDirections.transpose()*trial->fe().periodicShift(node)/length;
                for(Eigen::Index j=0;j<shift.size();++j)
                {
                    modes.block<3,3>(3*node,3*j)=shift(j)*Eigen::Matrix3d::Identity();
                }
            }
            return modes;
        }
        Eigen::Matrix3d gradient(const Eigen::VectorXd& coefficients) const
        {
            Eigen::Matrix3d result=Eigen::Matrix3d::Zero();
            for(Eigen::Index j=0;j<periodicDirections.cols();++j)
            {
                result+=coefficients.segment<3>(3*j)*periodicDirections.col(j).transpose()/length;
            }
            return result;
        }
        Vector6 averageStrain() const
        {
            return meanStrain*trial->dofVector();
        }
    };
}
#endif
