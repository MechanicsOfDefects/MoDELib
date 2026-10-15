/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_EllipsoidalInclusion_cpp_
#define model_EllipsoidalInclusion_cpp_

#include <EllipsoidalInclusion.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace
{
    Eigen::Matrix3d checkedStrain(const Eigen::Matrix3d &e, const double nu, const double mu)
    {
        if(!e.allFinite() || !std::isfinite(nu) || !std::isfinite(mu) || mu <= 0.0 || nu <= -1.0 || nu >= 0.5)
        {
            throw std::invalid_argument("EllipsoidalInclusion: invalid strain or elastic constants");
        }
        const double scale = e.cwiseAbs().maxCoeff();
        if(scale > 0.0 && (e / scale - e.transpose() / scale).norm() > 1.0e-12)
        {
            throw std::invalid_argument("EllipsoidalInclusion: eigenstrain must be symmetric");
        }
        return (0.5 * (e + e.transpose())).eval();
    }
} // namespace

namespace model
{
    template <int dim>
    EllipsoidalInclusion<dim>::EllipsoidalInclusion(const VectorDim &_C, const MatrixDim &_axes, const MatrixDim &_eT,
                                                    const double &_nu, const double &_mu,
                                                    const double &_mobilityReduction, const int &_phaseID,
                                                    const std::shared_ptr<SecondPhase<dim>> &sph,
                                                    const double quadratureTolerance,
                                                    const unsigned int quadratureMaxDepth)
        : EshelbyInclusionBase<dim>(checkedStrain(_eT, _nu, _mu), _nu, _mu, _mobilityReduction, _phaseID, sph), C(_C),
          axes(_axes), lengthScale(0.0), rotation(_axes), axisProduct(0.0), tolerance(quadratureTolerance),
          maxDepth(quadratureMaxDepth)
    {
        if(!C.allFinite() || !axes.allFinite() || !std::isfinite(tolerance) || tolerance < 1.0e-14 ||
           tolerance > 1.0e-3 || maxDepth == 0 || maxDepth > 50)
        {
            throw std::invalid_argument("EllipsoidalInclusion: invalid geometry or quadrature settings");
        }
        VectorDim lengths;
        for(int i = 0; i < 3; ++i)
        {
            lengths(i) = rotation.col(i).stableNorm();
            if(!(lengths(i) > 0.0) || !std::isfinite(lengths(i)))
            {
                throw std::invalid_argument("EllipsoidalInclusion: semiaxes must have positive finite lengths");
            }
            rotation.col(i) /= lengths(i);
        }
        if((rotation.transpose() * rotation - MatrixDim::Identity()).norm() > 1.0e-12)
        {
            throw std::invalid_argument("EllipsoidalInclusion: axes columns must be mutually orthogonal");
        }
        // Remove roundoff-size nonorthogonality, preserving handedness.
        rotation.col(1) =
            (rotation.col(1) - rotation.col(0) * rotation.col(0).dot(rotation.col(1))).normalized().eval();
        VectorDim third = rotation.col(0).cross(rotation.col(1));
        if(third.dot(rotation.col(2)) < 0.0)
        {
            third = -third;
        }
        rotation.col(2) = third;
        lengthScale = lengths.maxCoeff();
        lengths /= lengthScale;
        // A deliberate conditioning bound; do not silently return unreliable stresses.
        if(lengths.minCoeff() < 1.0e-6)
        {
            throw std::invalid_argument("EllipsoidalInclusion: semiaxis ratio exceeds 1e6");
        }
        axisSquared = lengths.array().square().matrix();
        axisProduct = lengths.prod();
        interiorIntegrals = integrate(0.0);
        interiorStressCache = std::make_shared<const InteriorStressCache>(
            InteriorStressCache{this->eT, evaluateStress(VectorDim::Zero(), true, this->eT)});
    }

    template <int dim>
    typename EllipsoidalInclusion<dim>::VectorDim EllipsoidalInclusion<dim>::localPosition(const VectorDim &x) const
    {
        if(!x.allFinite())
        {
            throw std::invalid_argument("EllipsoidalInclusion: nonfinite observation point");
        }
        const VectorDim r = rotation.transpose() * ((x - C) / lengthScale);
        if(!r.allFinite() || !std::isfinite(r.squaredNorm()))
        {
            throw std::overflow_error("EllipsoidalInclusion: observation point exceeds numerical range");
        }
        return r;
    }

    template <int dim> bool EllipsoidalInclusion<dim>::contains(const VectorDim &x) const
    {
        const VectorDim r = localPosition(x);
        return (r.array().square() / axisSquared.array()).sum() < 1.0;
    }

    template <int dim> double EllipsoidalInclusion<dim>::ellipsoidalCoordinate(const VectorDim &r) const
    {
        // f(ell)=sum r_i^2/(a_i^2+ell)-1 is decreasing and convex.
        // max(0, |r|^2-max(a_i^2)) <= root <= |r|^2.
        const VectorDim rSquared = r.array().square().matrix();
        const double radiusSquared = rSquared.sum();
        const double smallestAxisSquared = axisSquared.minCoeff();
        const double relativeTolerance = 8.0 * std::numeric_limits<double>::epsilon();
        double lo = std::max(0.0, radiusSquared - axisSquared.maxCoeff());
        double hi = radiusSquared;
        double ell = lo;
        double previousStep = hi - lo;
        for(unsigned int k = 0; k < 256; ++k)
        {
            const double positionTolerance = relativeTolerance * (smallestAxisSquared + ell);
            if(hi - lo <= positionTolerance)
            {
                return lo + 0.5 * (hi - lo);
            }
            const VectorDim denominator = (axisSquared.array() + ell).matrix();
            const VectorDim terms = (rSquared.array() / denominator.array()).matrix();
            const double f = terms.sum() - 1.0;
            const double minusDerivative = (terms.array() / denominator.array()).sum();
            const double correction = f / minusDerivative;
            // Do not require the bracket to collapse after Newton has converged.
            if(f == 0.0 || (std::abs(f) <= relativeTolerance && std::abs(correction) <= positionTolerance))
            {
                return ell;
            }
            if(f > 0.0)
            {
                lo = ell;
            }
            else
            {
                hi = ell;
            }
            const double newton = ell + correction;
            double next = newton;
            if(!std::isfinite(newton) || newton <= lo || newton >= hi || std::abs(correction) > 0.5 * previousStep)
            {
                // Preserve a valid bracket when Newton is unsafe or too slow.
                next = lo + 0.5 * (hi - lo);
            }
            previousStep = std::abs(next - ell);
            ell = next;
        }
        throw std::runtime_error("EllipsoidalInclusion: ellipsoidal coordinate did not converge");
    }

    template <int dim>
    typename EllipsoidalInclusion<dim>::Integrals EllipsoidalInclusion<dim>::integrate(const double ell) const
    {
        // J_i=int ds/((a_i^2+s) Delta), K_ij=int s ds/((a_i^2+s)(a_j^2+s) Delta).
        // s=(1+ell)/t^2-1 maps infinity to t=0, eliminating the endpoint singularity.
        // Integrate values scaled by (1+ell)^(3/2) so far-field accuracy is relative.
        const double h = 1.0 + ell;
        const auto integrand = [&](double t) -> Integrals
        {
            const double t2 = t * t, oneMinusT2 = (1.0 - t) * (1.0 + t);
            // This form avoids cancellation for very small semiaxes at t near 1.
            const VectorDim d = (oneMinusT2 + t2 * (axisSquared.array() + ell) / h).matrix();
            const double w = 2.0 * t2 / std::sqrt(d.prod());
            Integrals f;
            for(int i = 0; i < 3; ++i)
            {
                f(i) = w / d(i);
                for(int j = 0; j < 3; ++j)
                {
                    f(3 + 3 * i + j) = w * (oneMinusT2 + t2 * ell / h) / (d(i) * d(j));
                }
            }
            return f;
        };
        // Embedded Gauss 7 / Kronrod 15 rule, all components positive.
        const double z[8] = {0.9914553711208126, 0.9491079123427585, 0.8648644233597691, 0.7415311855993945,
                             0.5860872354676911, 0.4058451513773972, 0.2077849550078985, 0.0};
        const double wk[8] = {0.02293532201052922, 0.06309209262997855, 0.1047900103222502, 0.1406532597155259,
                              0.1690047266392679,  0.1903505780647854,  0.2044329400752989, 0.2094821410847278};
        const double wg[4] = {0.1294849661688697, 0.2797053914892767, 0.3818300505051189, 0.4179591836734694};
        unsigned int evaluations = 0;
        std::function<Integrals(double, double, unsigned int)> adaptive;
        adaptive = [&](double left, double right, unsigned int depth) -> Integrals
        {
            if(++evaluations > 200000)
            {
                throw std::runtime_error("EllipsoidalInclusion: quadrature work limit exceeded");
            }
            const double mid = 0.5 * (left + right), half = 0.5 * (right - left);
            const Integrals center = integrand(mid);
            Integrals fine = wk[7] * center, coarse = wg[3] * center;
            for(int k = 0; k < 7; ++k)
            {
                const Integrals pair = integrand(mid - half * z[k]) + integrand(mid + half * z[k]);
                fine += wk[k] * pair;
                if(k % 2 == 1)
                {
                    coarse += wg[(k - 1) / 2] * pair;
                }
            }
            fine *= half;
            coarse *= half;
            if(!fine.allFinite())
            {
                throw std::runtime_error("EllipsoidalInclusion: nonfinite quadrature");
            }
            if(((fine - coarse).array().abs() <= tolerance * ((right - left) + fine.array().abs())).all())
            {
                return fine;
            }
            if(depth == maxDepth)
            {
                throw std::runtime_error("EllipsoidalInclusion: quadrature failed to converge");
            }
            const Integrals first = adaptive(left, mid, depth + 1);
            return (first + adaptive(mid, right, depth + 1)).eval();
        };
        return (adaptive(0.0, 1.0, 0) / h / std::sqrt(h)).eval();
    }

    template <int dim>
    typename EllipsoidalInclusion<dim>::MatrixDim EllipsoidalInclusion<dim>::stress(const VectorDim &x) const
    {
        const VectorDim r = localPosition(x);
        const bool inside = (r.array().square() / axisSquared.array()).sum() <= 1.0;
        if(inside)
        {
            const auto cache = std::atomic_load(&interiorStressCache);
            if((cache->eigenstrain.array() == this->eT.array()).all())
            {
                return cache->stress;
            }
            // The base exposes mutable eT. Refresh only after successful validation.
            const MatrixDim eigenstrain = this->eT;
            const MatrixDim sigma = evaluateStress(r, true, eigenstrain);
            const auto updated = std::make_shared<const InteriorStressCache>(InteriorStressCache{eigenstrain, sigma});
            std::atomic_store(&interiorStressCache, updated);
            return sigma;
        }
        return evaluateStress(r, false, this->eT);
    }

    template <int dim>
    typename EllipsoidalInclusion<dim>::MatrixDim
    EllipsoidalInclusion<dim>::evaluateStress(const VectorDim &r, const bool inside, const MatrixDim &eigenstrain) const
    {
        const MatrixDim e = checkedStrain(eigenstrain, this->nu, this->mu);
        // Derive from eT rather than the potentially stale base-class pT.
        const MatrixDim p =
            rotation.transpose() *
            (2.0 * this->mu * (e + this->nu / (1.0 - 2.0 * this->nu) * e.trace() * MatrixDim::Identity())) * rotation;
        if(!p.allFinite())
        {
            throw std::overflow_error("EllipsoidalInclusion: eigenstress overflow");
        }
        if(p.cwiseAbs().maxCoeff() == 0.0)
        {
            return MatrixDim::Zero();
        }
        const double ell = inside ? 0.0 : ellipsoidalCoordinate(r);
        const Integrals values = inside ? interiorIntegrals : integrate(ell);
        const double pi = std::acos(-1.0), factor = pi * axisProduct;
        MatrixDim H = MatrixDim::Zero(), T = MatrixDim::Zero();
        for(int i = 0; i < 3; ++i)
        {
            H(i, i) = -2.0 * factor * values(i);
            for(int j = 0; j < 3; ++j)
            {
                T(i, j) = -4.0 * factor * p(i, j) * values(3 + 3 * i + j);
            }
            for(int j = 0; j < 3; ++j)
            {
                T(i, i) -= 2.0 * factor * p(j, j) * values(3 + 3 * i + j);
            }
        }
        if(!inside)
        {
            const VectorDim b = (axisSquared.array() + ell).inverse().matrix();
            const MatrixDim B = b.asDiagonal();
            const VectorDim q = B * r;
            const double d = q.squaredNorm();
            const double invDelta = std::sqrt(b.prod());
            const VectorDim l = 2.0 * q / d;
            const VectorDim Bq = B * q;
            const MatrixDim ll = l * l.transpose();
            const MatrixDim lH = (2.0 * B - 2.0 * (Bq * l.transpose() + l * Bq.transpose()) + 2.0 * q.dot(Bq) * ll) / d;
            H += 4.0 * factor * invDelta / d * (q * q.transpose());
            // Exact moving-lower-limit terms for Hessian(p:Hessian(Psi)).
            // At s=ell, f=1-r.B.r=0, F=s/Delta*(8 q.p.q-4 f tr(pB)).
            const double trPB = (p.diagonal().array() * b.array()).sum();
            const double v = 8.0 * q.dot(p * q);
            const double w = ell * invDelta;
            const VectorDim Fx = w * (16.0 * B * p * q + 8.0 * trPB * q);
            const double Fs = invDelta * (1.0 - 0.5 * ell * b.sum()) * v + w * (-16.0 * Bq.dot(p * q) - 4.0 * d * trPB);
            T += 0.25 * factor * (Fx * l.transpose() + l * Fx.transpose() + Fs * ll + w * v * lH);
        }
        // Kelvin Green tensor: u=-[4(1-nu)p.grad(Phi)-grad(p:Hess(Psi))]/[16 pi mu(1-nu)].
        MatrixDim sigma = -(2.0 * (1.0 - this->nu) * (p * H + H * p) - T +
                            2.0 * this->nu * (p.cwiseProduct(H)).sum() * MatrixDim::Identity()) /
                          (8.0 * pi * (1.0 - this->nu));
        if(inside)
        {
            sigma -= p;
        }
        sigma = (rotation * sigma * rotation.transpose()).eval();
        if(!sigma.allFinite())
        {
            throw std::overflow_error("EllipsoidalInclusion: stress overflow");
        }
        return (0.5 * sigma + 0.5 * sigma.transpose()).eval();
    }

    template class EllipsoidalInclusion<3>;
} // namespace model
#endif
