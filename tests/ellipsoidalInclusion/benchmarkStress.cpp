// Build with EllipsoidalInclusion.cpp, reference/SphericalInclusion.cpp,
// and the MoDELib base (or tests/support headers). Use -O3 -DNDEBUG.
#include <EllipsoidalInclusion.h>
#include <SphericalInclusion.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

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

int main()
{
    M strain;
    strain << 0.01, 0.003, -0.002, 0.003, -0.006, 0.001, -0.002, 0.001, 0.004;
    const M Q = Eigen::AngleAxisd(0.71, V(1, 2, 3).normalized()).toRotationMatrix();
    std::mt19937 rng(123);
    std::normal_distribution<double> normal;
    std::vector<V> directions;
    for(int k = 0; k < 512; ++k)
    {
        const V direction(normal(rng), normal(rng), normal(rng));
        directions.push_back(direction.normalized());
    }
    std::cout << "axes,region,sphere_ns,ellipsoid_ns,ratio,ellipsoid_min_ns,ellipsoid_max_ns\n";
    for(const V &lengths : {V(1, 1, 1), V(3, 2, 1), V(10, 1, 0.1), V(100, 1, 1)})
    {
        // The spherical reference radius is the shortest semiaxis so every
        // comparison uses the same interior/exterior branch at the same point.
        model::SphericalInclusion<3> sphere(V::Zero(), lengths.minCoeff(), strain, 0.31, 43.0, 1.0, 0, nullptr);
        model::EllipsoidalInclusion<3> ellipsoid(V::Zero(), Q * lengths.asDiagonal(), strain, 0.31, 43.0, 1.0, 0,
                                                 nullptr, 1.e-10);
        for(const double radialFactor : {0.5, 1.000001, 2.0, 10.0})
        {
            std::vector<V> points;
            for(const auto &direction : directions)
            {
                if(radialFactor < 1.0)
                {
                    points.push_back(0.5 * lengths.minCoeff() * direction);
                }
                else
                {
                    points.push_back(Q * (radialFactor * lengths.cwiseProduct(direction)));
                }
            }
            measure(sphere, points);
            measure(ellipsoid, points);
            std::vector<double> sphereTimes, ellipsoidTimes;
            for(int repeat = 0; repeat < 7; ++repeat)
            {
                if(repeat % 2 == 0)
                {
                    sphereTimes.push_back(measure(sphere, points));
                    ellipsoidTimes.push_back(measure(ellipsoid, points));
                }
                else
                {
                    ellipsoidTimes.push_back(measure(ellipsoid, points));
                    sphereTimes.push_back(measure(sphere, points));
                }
            }
            std::sort(sphereTimes.begin(), sphereTimes.end());
            std::sort(ellipsoidTimes.begin(), ellipsoidTimes.end());
            const char *region = radialFactor < 1.0   ? "interior"
                                 : radialFactor < 1.1 ? "near_surface"
                                 : radialFactor < 3.0 ? "2x_surface"
                                                      : "10x_surface";
            std::cout << lengths(0) << ':' << lengths(1) << ':' << lengths(2) << ',' << region << ',' << std::fixed
                      << std::setprecision(2) << sphereTimes[3] << ',' << ellipsoidTimes[3] << ','
                      << ellipsoidTimes[3] / sphereTimes[3] << ',' << ellipsoidTimes.front() << ','
                      << ellipsoidTimes.back() << std::endl;
        }
    }
}
