// Build with EllipsoidalInclusion.cpp, reference/SphericalInclusion.cpp,
// and the MoDELib base (or tests/support headers). Use -O3 -DNDEBUG.
#include <InclusionMicrostructure/EllipsoidalInclusion.h>
#include <InclusionMicrostructure/SphericalInclusion.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <string>

using V = Eigen::Vector3d;
using M = Eigen::Matrix3d;
using Base = model::EshelbyInclusionBase<3>;
using Clock = std::chrono::steady_clock;
volatile double benchmarkSink = 0.0;

// Separate translation units and no LTO keep stress() calls out of the timing loop.
double measure(const Base &inclusion, const std::vector<V> &points)
{
    const auto start = Clock::now();
    double checksum = 0.0;
    std::size_t calls = 0;
    double seconds = 0.0;
    do
    {
        for(const auto &point : points)
        {
            checksum += inclusion.stress(point).sum();
        }
        calls += points.size();
        seconds = std::chrono::duration<double>(Clock::now() - start).count();
    } while(seconds < 0.02);
    benchmarkSink = checksum;
    return seconds * 1.e9 / calls;
}

struct Errors
{
    std::vector<double> loose;
    double defaultMaximum = 0.0;
    double sphereMaximum = 0.0;
    double interiorMaximum = 0.0;
    double worstFactor = 0.0;
    double worstNu = 0.0;
    int worstMode = 0;
};
int main(int argc, char **argv)
{
    const double trialTolerance = argc > 1 ? std::stod(argv[1]) : 1.e-4;
    const M Q = Eigen::AngleAxisd(0.71, V(1, 2, 3).normalized()).toRotationMatrix();
    const std::vector<V> geometries = {V(1, 1, 1), V(3, 2, 1), V(10, 1, 0.1), V(100, 1, 1), V(1, 1, 0.01)};
    std::mt19937 rng(123);
    std::normal_distribution<double> normal;
    std::vector<V> directions;
    for(int k = 0; k < 512; ++k)
    {
        const V direction(normal(rng), normal(rng), normal(rng));
        directions.push_back(direction.normalized());
    }
    M mixed;
    mixed << 0.01, 0.003, -0.002, 0.003, -0.006, 0.001, -0.002, 0.001, 0.004;
    std::cout << "SPEED\naxes,region,sphere_ns,baseline_1e10_ns,loose_ns,speedup,loose_over_sphere\n";
    for(const V &lengths : geometries)
    {
        model::SphericalInclusion<3> sphere(V::Zero(), lengths.minCoeff(), mixed, 0.31, 43.0, 1, 0, nullptr);
        model::EllipsoidalInclusion<3> fine(V::Zero(), Q * lengths.asDiagonal(), mixed, 0.31, 43.0, 1, 0, nullptr,
                                            1.e-10);
        model::EllipsoidalInclusion<3> loose(V::Zero(), Q * lengths.asDiagonal(), mixed, 0.31, 43.0, 1, 0, nullptr,
                                             trialTolerance);
        for(double factor : {0.5, 1.000001, 2.0, 10.0})
        {
            std::vector<V> points;
            for(const V &direction : directions)
            {
                if(factor < 1)
                {
                    points.push_back(0.5 * lengths.minCoeff() * direction);
                }
                else
                {
                    points.push_back(Q * (factor * lengths.cwiseProduct(direction)));
                }
            }
            measure(sphere, points);
            measure(fine, points);
            measure(loose, points);
            std::vector<double> st, ft, lt;
            for(int trial = 0; trial < 7; ++trial)
            {
                st.push_back(measure(sphere, points));
                if(trial % 2 == 0)
                {
                    ft.push_back(measure(fine, points));
                    lt.push_back(measure(loose, points));
                }
                else
                {
                    lt.push_back(measure(loose, points));
                    ft.push_back(measure(fine, points));
                }
            }
            std::sort(st.begin(), st.end());
            std::sort(ft.begin(), ft.end());
            std::sort(lt.begin(), lt.end());
            const char *region = factor < 1     ? "interior"
                                 : factor < 1.1 ? "near_surface"
                                 : factor < 3   ? "2x_surface"
                                                : "10x_surface";
            std::cout << lengths.transpose() << ',' << region << ',' << std::setprecision(8) << st[3] << ',' << ft[3]
                      << ',' << lt[3] << ',' << ft[3] / lt[3] << ',' << lt[3] / st[3] << std::endl;
        }
    }
    std::cout
        << "ACCURACY\naxes,count,max_relative,p95_relative,baseline_1e10_max_relative,interior_max_relative,sphere_max_"
           "relative,worst_factor,worst_nu,worst_mode\n";
    for(const V &lengths : geometries)
    {
        Errors errors;
        for(double nu : {-0.8, 0.31, 0.49})
        {
            for(int mode = 0; mode < 8; ++mode)
            {
                M strain = M::Zero();
                if(mode < 3)
                {
                    strain(mode, mode) = 0.01;
                }
                else if(mode < 6)
                {
                    const int i = mode == 5 ? 1 : 0;
                    const int j = mode == 3 ? 1 : 2;
                    strain(i, j) = strain(j, i) = 0.01;
                }
                else if(mode == 6)
                {
                    strain = 0.01 * M::Identity();
                }
                else
                {
                    strain = mixed;
                }
                model::EllipsoidalInclusion<3> ref(V::Zero(), Q * lengths.asDiagonal(), strain, nu, 43.0, 1, 0, nullptr,
                                                   1.e-12, 40);
                model::EllipsoidalInclusion<3> fine(V::Zero(), Q * lengths.asDiagonal(), strain, nu, 43.0, 1, 0,
                                                    nullptr, 1.e-10);
                model::EllipsoidalInclusion<3> loose(V::Zero(), Q * lengths.asDiagonal(), strain, nu, 43.0, 1, 0,
                                                     nullptr, trialTolerance);
                model::SphericalInclusion<3> sphere(V::Zero(), 1.0, strain, nu, 43.0, 1, 0, nullptr);
                for(double factor : {0.0, 1.0 + 1.e-12, 1.0 + 1.e-8, 1.001, 1.1, 2.0, 10.0, 1.e3, 1.e6})
                {
                    const int samples = factor == 0 ? 1 : 32;
                    for(int k = 0; k < samples; ++k)
                    {
                        V direction;
                        if(k < 3)
                        {
                            direction = V::Unit(k);
                        }
                        else
                        {
                            direction = directions[k];
                        }
                        const V point = Q * (factor * lengths.cwiseProduct(direction));
                        const M reference = ref.stress(point), approx = loose.stress(point);
                        const double norm = std::max(1.e-300, reference.norm());
                        const double error = (approx - reference).norm() / norm;
                        const double defaultError = (fine.stress(point) - reference).norm() / norm;
                        if(!std::isfinite(error) || !std::isfinite(defaultError))
                        {
                            throw std::runtime_error("Nonfinite accuracy metric");
                        }
                        if(factor == 0)
                        {
                            errors.interiorMaximum = std::max(errors.interiorMaximum, error);
                        }
                        else
                        {
                            if(errors.loose.empty() ||
                               error > *std::max_element(errors.loose.begin(), errors.loose.end()))
                            {
                                errors.worstFactor = factor;
                                errors.worstNu = nu;
                                errors.worstMode = mode;
                            }
                            errors.loose.push_back(error);
                            errors.defaultMaximum = std::max(errors.defaultMaximum, defaultError);
                        }
                        if((lengths - V::Ones()).norm() == 0.0)
                        {
                            errors.sphereMaximum =
                                std::max(errors.sphereMaximum, (approx - sphere.stress(point)).norm() / norm);
                        }
                    }
                }
            }
        }
        std::sort(errors.loose.begin(), errors.loose.end());
        std::cout << lengths.transpose() << ',' << errors.loose.size() << ',' << std::setprecision(12)
                  << errors.loose.back() << ','
                  << errors.loose[static_cast<std::size_t>(0.95 * (errors.loose.size() - 1))] << ','
                  << errors.defaultMaximum << ',' << errors.interiorMaximum << ',' << errors.sphereMaximum << ','
                  << errors.worstFactor << ',' << errors.worstNu << ',' << errors.worstMode << std::endl;
    }
}
