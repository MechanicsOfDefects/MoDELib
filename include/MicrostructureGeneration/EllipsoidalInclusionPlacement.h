/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_EllipsoidalInclusionPlacement_H_
#define model_EllipsoidalInclusionPlacement_H_
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace model::ellipsoidalPlacement
{
    // Positive shape tensors A*A^T; touching counts as overlap.
    inline bool overlap(const Eigen::Vector3d &displacement, const Eigen::Matrix3d &first,
                        const Eigen::Matrix3d &second)
    {
        if(!displacement.allFinite() || !first.allFinite() || !second.allFinite() ||
           first.llt().info() != Eigen::Success || second.llt().info() != Eigen::Success)
        {
            throw std::invalid_argument("Ellipsoid overlap requires finite positive shape tensors");
        }
        // The ellipsoid contact function is concave on [0,1]. Its maximum
        // exceeds one exactly when the two closed ellipsoids are separated.
        const double scale = std::max(first.norm(), second.norm());
        const Eigen::Matrix3d a = first / scale, b = second / scale;
        const Eigen::Vector3d d = displacement / std::sqrt(scale);
        const auto contact = [&](double t)
        {
            const Eigen::Matrix3d shape = (1 - t) * a + t * b;
            return t * (1 - t) * d.dot(shape.ldlt().solve(d));
        };
        double lo = 0, hi = 1;
        constexpr double ratio = 0.6180339887498948482;
        double left = hi - ratio * (hi - lo), right = lo + ratio * (hi - lo);
        double fl = contact(left), fr = contact(right);
        for(int k = 0; k < 100; ++k)
        {
            if(std::max(fl, fr) > 1.0 + 1.e-12)
            {
                return false;
            }
            if(fl < fr)
            {
                lo = left;
                left = right;
                fl = fr;
                right = lo + ratio * (hi - lo);
                fr = contact(right);
            }
            else
            {
                hi = right;
                right = left;
                fr = fl;
                left = hi - ratio * (hi - lo);
                fl = contact(left);
            }
        }
        return std::max(fl, fr) <= 1.0 + 1.e-12;
    }

    // Squared distance from the origin to a triangle in unit-ellipsoid coordinates.
    inline double triangleDistanceSquared(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c)
    {
        const auto edgeDistance = [](const Eigen::Vector3d &p, const Eigen::Vector3d &q)
        {
            const Eigen::Vector3d edge = q - p;
            const double t = edge.squaredNorm() > 0 ? std::clamp(-p.dot(edge) / edge.squaredNorm(), 0.0, 1.0) : 0.0;
            return (p + t * edge).squaredNorm();
        };
        double distance = std::min({edgeDistance(a, b), edgeDistance(b, c), edgeDistance(c, a)});
        const Eigen::Vector3d u = b - a, v = c - a, n = u.cross(v);
        if(n.squaredNorm() > 0)
        {
            const double s = (-a).cross(v).dot(n) / n.squaredNorm();
            const double t = u.cross(-a).dot(n) / n.squaredNorm();
            if(s >= 0 && t >= 0 && s + t <= 1)
            {
                distance = std::min(distance, (a + s * u + t * v).squaredNorm());
            }
        }
        return distance;
    }
} // namespace model::ellipsoidalPlacement
#endif
