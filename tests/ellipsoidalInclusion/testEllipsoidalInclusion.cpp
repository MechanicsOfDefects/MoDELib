#include <InclusionMicrostructure/EllipsoidalInclusion.h>
#ifdef COMPARE_SUPPLIED_SPHERE
#include <InclusionMicrostructure/SphericalInclusion.h>
#endif
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <random>
#include <limits>
#include <atomic>
#include <thread>
#include <vector>
using M = Eigen::Matrix3d;
using V = Eigen::Vector3d;
using E = model::EllipsoidalInclusion<3>;
M makeAxes(const V &a1, const V &a2, const V &a3)
{
    M axes;
    axes.col(0) = a1;
    axes.col(1) = a2;
    axes.col(2) = a3;
    return axes;
}
const double nu = 0.31, mu = 43.0;
int checks = 0;
double worstSphere = 0.0, worstTraction = 0.0;
void check(bool ok, const char *msg)
{
    ++checks;
    if(!ok)
    {
        throw std::runtime_error(msg);
    }
}
void close(const M &a, const M &b, double tol, const char *msg)
{
    const double error = (a - b).norm() / std::max(1.e-30, b.norm());
    if(error > tol)
    {
        std::cerr << msg << ": relative error=" << error << "\n";
    }
    check(error <= tol, msg);
}
M sphere(const V &r, double a, const M &e, double n = nu, double m = mu)
{
    const M I = M::Identity(), p = 2 * m * (e + n / (1 - 2 * n) * e.trace() * I);
    const double R = r.norm(), L = (5 * n - 1) / (15 * (1 - n)), B = (4 - 5 * n) / (15 * (1 - n));
    if(R <= a)
    {
        return 2 * m * ((L + n / (1 - 2 * n) * (2 * B + 3 * L)) * e.trace() * I + 2 * B * e) - p;
    }
    const double R2 = R * R, c = a * a / R2;
    const V pr = p * r;
    const double rr = r.dot(pr), tr = p.trace();
    return std::pow(a / R, 3) / (2 * (1 - n)) *
           (((10 * (1 - 2 * n) + 6 * c) / 15) * p + (2 * n - 2 * c) / R2 * (pr * r.transpose() + r * pr.transpose()) +
            ((3 * c - 5 * (1 - 2 * n)) / 15 * tr + (1 - 2 * n - c) / R2 * rr) * I +
            (-(5 - 7 * c) / (R2 * R2) * rr + (1 - c) / R2 * tr) * r * r.transpose());
}
int main()
{
    try
    {
        M e;
        e << 0.01, 0.003, -0.002, 0.003, -0.006, 0.001, -0.002, 0.001, 0.004;
        std::mt19937 gen(21);
        std::normal_distribution<double> normal;
        const M Q = Eigen::AngleAxisd(0.71, V(1, 2, 3).normalized()).toRotationMatrix();
        const V center(0.2, -0.3, 0.1);
        for(const double n : {-0.8, 0.0, 0.31, 0.49})
        {
            for(int kind = 0; kind < 8; ++kind)
            {
                M strain = M::Zero();
                if(kind < 3)
                {
                    strain(kind, kind) = 0.01;
                }
                else if(kind < 6)
                {
                    int i = kind == 3 ? 0 : kind == 4 ? 0 : 1;
                    int j = kind == 3 ? 1 : 2;
                    strain(i, j) = strain(j, i) = 0.01;
                }
                else
                {
                    strain = kind == 6 ? (0.01 * M::Identity()).eval() : e;
                }
                E ell(center, 2 * Q, strain, n, mu, 1.0, 0, nullptr);
#ifdef COMPARE_SUPPLIED_SPHERE
                model::SphericalInclusion<3> original(center, 2, strain, n, mu, 1.0, 0, nullptr);
#endif
                for(double radius : {0.0, 0.7, 1.999, 2.0000000002, 2.2, 4.0, 20.0, 2000.0, 2.e6})
                {
                    V dir(normal(gen), normal(gen), normal(gen));
                    dir.normalize();
                    const V x = center + radius * dir;
                    const M ref = sphere(x - center, 2, strain, n, mu), got = ell.stress(x);
                    worstSphere = std::max(worstSphere, (got - ref).norm() / ref.norm());
                    close(got, ref, 3.e-9, "sphere formula");
#ifdef COMPARE_SUPPLIED_SPHERE
                    close(got, original.stress(x), 3.e-9, "supplied sphere class");
#endif
                }
            }
        }
        E boundary(V::Zero(), M::Identity(), e, nu, mu, 1, 0, nullptr);
        check(!boundary.contains(V::UnitX()), "contains excludes exact boundary");
        close(boundary.stress(V::UnitX()), sphere(V::UnitX(), 1, e), 1.e-10, "interior boundary convention");
        for(const V &lengths : {V(3, 2, 1), V(10, 1, 0.1), V(1, 1, 0.01), V(100, 1, 1)})
        {
            E local(V::Zero(), lengths.asDiagonal(), e, nu, mu, 1, 0, nullptr, 1.e-10);
            E rotated(center, Q * lengths.asDiagonal(), Q * e * Q.transpose(), nu, mu, 1, 0, nullptr, 1.e-10);
            for(int k = 0; k < 12; ++k)
            {
                V dir(normal(gen), normal(gen), normal(gen));
                dir.normalize();
                const V surf = lengths.cwiseProduct(dir), inside = 0.3 * surf, outside = 1.6 * surf;
                close(local.stress(inside), local.stress(V::Zero()), 1.e-13, "uniform interior");
                close(rotated.stress(center + Q * outside), Q * local.stress(outside) * Q.transpose(), 2.e-8,
                      "rotation covariance");
                const V normalTo = (dir.array() / lengths.array()).matrix().normalized();
                const V ti = local.stress((1 - 1.e-10) * surf) * normalTo;
                const V to = local.stress((1 + 1.e-10) * surf) * normalTo;
                const double jump = (ti - to).norm() / std::max(1.e-30, local.stress(V::Zero()).norm());
                worstTraction = std::max(worstTraction, jump);
                if(jump >= 2.e-6)
                {
                    std::cerr << "lengths " << lengths.transpose() << " jump " << jump << "\n";
                }
                check(jump < 2.e-6, "traction continuity");
                // Equilibrium in the exterior, using independent central spatial differences.
                const double step = lengths.minCoeff() * 1.e-4;
                V divergence = V::Zero();
                for(int j = 0; j < 3; ++j)
                {
                    divergence += (local.stress(outside + step * V::Unit(j)).col(j) -
                                   local.stress(outside - step * V::Unit(j)).col(j)) /
                                  (2 * step);
                }
                check(divergence.norm() * lengths.minCoeff() / local.stress(V::Zero()).norm() < 2.e-6,
                      "exterior equilibrium");
            }
        }
        for(double scale : {1.e-9, 1.0, 1.e9})
        {
            E ell(V::Zero(), makeAxes(scale * V(3, 0, 0), scale * V(0, 2, 0), scale * V(0, 0, 1)), e, nu, mu, 1, 0,
                  nullptr);
            E ref(V::Zero(), makeAxes(V(3, 0, 0), V(0, 2, 0), V(0, 0, 1)), e, nu, mu, 1, 0, nullptr);
            close(ell.stress(scale * V(4, 3, 2)), ref.stress(V(4, 3, 2)), 1.e-10, "scale invariance");
        }
        // Independent dilatational sphere values, beyond the general sphere formula.
        E dil(V::Zero(), M::Identity(), 0.01 * M::Identity(), nu, mu, 1, 0, nullptr);
        close(dil.stress(V::Zero()), (-4 * mu * (1 + nu) * 0.01 / (3 * (1 - nu)) * M::Identity()).eval(), 1.e-10,
              "hydrostatic sphere interior");
        const V observation(2, 3, 4), unit = observation.normalized();
        close(dil.stress(observation),
              (2 * mu * (1 + nu) * 0.01 / (3 * (1 - nu) * std::pow(observation.norm(), 3)) *
               (M::Identity() - 3 * unit * unit.transpose()))
                  .eval(),
              1.e-10, "dilatational sphere exterior");
        E regular(V::Zero(), makeAxes(V(3, 0, 0), V(0, 2, 0), V(0, 0, 1)), e, nu, mu, 1, 0, nullptr);
        E permuted(V::Zero(), makeAxes(V(0, 0, -1), V(-3, 0, 0), V(0, 2, 0)), e, nu, mu, 1, 0, nullptr);
        close(regular.stress(V(4, 3, 2)), permuted.stress(V(4, 3, 2)), 1.e-10,
              "semiaxis sign and permutation invariance");
        E tight(V::Zero(), makeAxes(V(10, 0, 0), V(0, 1, 0), V(0, 0, 0.1)), e, nu, mu, 1, 0, nullptr, 1.e-12, 40);
        E normalTolerance(V::Zero(), makeAxes(V(10, 0, 0), V(0, 1, 0), V(0, 0, 0.1)), e, nu, mu, 1, 0, nullptr, 1.e-10);
        for(const V &point : {V::Zero().eval(), V(10.001, 0, 0), V(0, 0, 0.100001), V(12, 2, 1)})
        {
            close(normalTolerance.stress(point), tight.stress(point), 1.e-8, "quadrature refinement");
        }
        // The default must match an explicit 1e-3, while tighter overrides remain usable.
        for(const V &lengths : {V(3, 2, 1), V(10, 1, 0.1), V(100, 1, 1), V(1, 1, 0.01)})
        {
            E defaultTolerance(V::Zero(), Q * lengths.asDiagonal(), e, nu, mu, 1, 0, nullptr);
            E explicitTolerance(V::Zero(), Q * lengths.asDiagonal(), e, nu, mu, 1, 0, nullptr, 1.e-3);
            E referenceTolerance(V::Zero(), Q * lengths.asDiagonal(), e, nu, mu, 1, 0, nullptr, 1.e-12, 40);
            for(double factor : {0.0, 1.000001, 2.0, 10.0})
            {
                const V point = Q * (factor * lengths.cwiseProduct(V(1, 2, 3).normalized()));
                const M actual = defaultTolerance.stress(point);
                check((actual - explicitTolerance.stress(point)).norm() == 0.0, "default quadrature tolerance is 1e-3");
                close(actual, referenceTolerance.stress(point), 1.e-3, "default accuracy on sampled points");
            }
        }
        regular.eT = 2 * e;
        close(regular.stress(V(4, 3, 2)), 2 * permuted.stress(V(4, 3, 2)), 1.e-10, "mutable base eigenstrain honored");
        const M originalInterior = permuted.stress(V::Zero());
        close(regular.stress(V::Zero()), 2 * originalInterior, 1.e-10, "interior cache refreshed after eT update");
        close(regular.stress(V(0.1, 0.1, 0.1)), 2 * originalInterior, 1.e-10, "refreshed interior cache reused");
        regular.eT.setZero();
        check(regular.stress(V::Zero()).norm() == 0.0, "interior cache refresh to zero strain");
        regular.eT = e;
        close(regular.stress(V::Zero()), originalInterior, 1.e-10, "interior cache refresh from zero strain");
        regular.eT(0, 1) += 0.1;
        bool rejectedCachedStrain = false;
        try
        {
            regular.stress(V::Zero());
        }
        catch(const std::invalid_argument &)
        {
            rejectedCachedStrain = true;
        }
        check(rejectedCachedStrain, "invalid updated strain cannot return cached stress");
        regular.eT = e;
        close(regular.stress(V::Zero()), originalInterior, 1.e-10, "cache survives rejected strain update");
        E copied = regular;
        copied.eT = 3 * e;
        close(copied.stress(V::Zero()), 3 * originalInterior, 1.e-10, "copied inclusion refreshes its own cache");
        close(regular.stress(V::Zero()), originalInterior, 1.e-10, "copy refresh leaves original cache unchanged");
        // Concurrent readers may refresh the same stale cache; eT stays fixed.
        regular.eT = 2 * e;
        std::atomic<bool> concurrentCacheOK(true);
        std::vector<std::thread> readers;
        for(int reader = 0; reader < 4; ++reader)
        {
            readers.emplace_back(
                [&]()
                {
                    try
                    {
                        for(int sample = 0; sample < 1000; ++sample)
                        {
                            if((regular.stress(V::Zero()) - 2 * originalInterior).norm() >
                               1.e-10 * originalInterior.norm())
                            {
                                concurrentCacheOK.store(false);
                            }
                        }
                    }
                    catch(const std::exception &)
                    {
                        concurrentCacheOK.store(false);
                    }
                });
        }
        for(auto &reader : readers)
        {
            reader.join();
        }
        check(concurrentCacheOK.load(), "concurrent interior cache refresh and reads");
        bool failedQuadrature = false;
        try
        {
            E limited(V::Zero(), makeAxes(V(100, 0, 0), V::UnitY(), V::UnitZ()), e, nu, mu, 1, 0, nullptr, 1.e-14, 1);
        }
        catch(const std::runtime_error &)
        {
            failedQuadrature = true;
        }
        check(failedQuadrature, "quadrature nonconvergence reported");
        E tiny(V::Zero(), M::Identity(), 1.e-12 * e, nu, mu, 1, 0, nullptr);
        close(tiny.stress(V(2, 1, 1)), 1.e-12 * sphere(V(2, 1, 1), 1, e), 1.e-10, "small eigenstrain preserved");
        E zero(V::Zero(), M::Identity(), M::Zero(), nu, mu, 1, 0, nullptr);
        check(zero.stress(V::Zero()).norm() == 0, "zero strain");
        int rejected = 0;
        try
        {
            E bad(V::Zero(), makeAxes(V::Zero(), V::UnitY(), V::UnitZ()), e, nu, mu, 1, 0, nullptr);
        }
        catch(const std::invalid_argument &)
        {
            ++rejected;
        }
        try
        {
            E bad(V::Zero(), makeAxes(V::UnitX(), V(1, 1, 0), V::UnitZ()), e, nu, mu, 1, 0, nullptr);
        }
        catch(const std::invalid_argument &)
        {
            ++rejected;
        }
        try
        {
            E bad(V::Zero(), M::Identity(), e, 0.5, mu, 1, 0, nullptr);
        }
        catch(const std::invalid_argument &)
        {
            ++rejected;
        }
        M nonsym = e;
        nonsym(0, 1) += 0.1;
        try
        {
            E bad(V::Zero(), M::Identity(), nonsym, nu, mu, 1, 0, nullptr);
        }
        catch(const std::invalid_argument &)
        {
            ++rejected;
        }
        check(rejected == 4, "invalid inputs rejected");
        // Check every column pair, independently of their lengths.
        for(int i = 0; i < 3; ++i)
        {
            for(int j = i + 1; j < 3; ++j)
            {
                M invalidAxes = V(3.0, 2.0, 1.0).asDiagonal();
                invalidAxes.col(j) += 0.01 * invalidAxes.col(i);
                bool rejectedAxes = false;
                try
                {
                    E bad(V::Zero(), invalidAxes, e, nu, mu, 1, 0, nullptr);
                }
                catch(const std::invalid_argument &)
                {
                    rejectedAxes = true;
                }
                check(rejectedAxes, "nonorthogonal matrix columns rejected");
            }
        }
        M nonfiniteAxes = M::Identity();
        nonfiniteAxes(1, 2) = std::numeric_limits<double>::quiet_NaN();
        bool rejectedNonfinite = false;
        try
        {
            E bad(V::Zero(), nonfiniteAxes, e, nu, mu, 1, 0, nullptr);
        }
        catch(const std::invalid_argument &)
        {
            rejectedNonfinite = true;
        }
        check(rejectedNonfinite, "nonfinite axes matrix rejected");
        const M inputAxes = Q * V(3.0, 2.0, 1.0).asDiagonal();
        E matrixInput(center, inputAxes, e, nu, mu, 1, 0, nullptr);
        check((matrixInput.axes - inputAxes).norm() == 0.0, "public axes preserves semiaxis columns");
        check(matrixInput.contains(center + 0.9 * inputAxes.col(0)), "matrix column length sets semiaxis");
        check(!matrixInput.contains(center + 1.1 * inputAxes.col(0)), "point beyond matrix semiaxis excluded");
        std::cout << "PASS " << checks << " checks; worst sphere relative error " << worstSphere
                  << "; worst normalized traction jump " << worstTraction << "\n";
        return 0;
    }
    catch(const std::exception &ex)
    {
        std::cerr << "FAIL: " << ex.what() << "\n";
        return 1;
    }
}
