#include <EllipsoidalInclusionPlacement.h>
#include <Eigen/Geometry>
#include <iostream>
#include <stdexcept>
using V = Eigen::Vector3d;
using M = Eigen::Matrix3d;
void check(bool ok, const char *message)
{
    if(!ok)
    {
        throw std::runtime_error(message);
    }
}
int main()
{
    using namespace model::ellipsoidalPlacement;
    try
    {
        const M sphere = M::Identity();
        check(overlap(V(1.9, 0, 0), sphere, sphere), "Overlapping spheres missed");
        check(overlap(V(2, 0, 0), sphere, sphere), "Touching spheres missed");
        check(!overlap(V(2.01, 0, 0), sphere, sphere), "Separated spheres rejected");
        check(overlap(V::Zero(), sphere, 100 * sphere), "Contained spheres missed");
        const M shape = V(9, 4, 1).asDiagonal();
        check(!overlap(V(6.1, 0, 0), shape, shape), "Separated ellipsoids rejected");
        check(overlap(V(5.9, 0, 0), shape, shape), "Overlapping ellipsoids missed");
        check(!overlap(V(0, 0, 2.1), shape, shape), "Thin-direction gap rejected");
        check(!overlap(V(4.01, 0, 0), shape, sphere), "Separated sphere/ellipsoid rejected");
        check(overlap(V(3.99, 0, 0), shape, sphere), "Sphere/ellipsoid overlap missed");
        const M rotation = Eigen::AngleAxisd(.7, V(1, 2, 3).normalized()).toRotationMatrix();
        check(!overlap(rotation * V(6.1, 0, 0), rotation * shape * rotation.transpose(),
                       rotation * shape * rotation.transpose()),
              "Rotation changed separation");
        check(!overlap(V(6.1e-8, 0, 0), shape * 1e-16, shape * 1e-16), "Scaling changed separation");
        // All vertices outside the unit ball, but the face crosses it.
        check(triangleDistanceSquared(V(-3, -3, .5), V(3, -3, .5), V(0, 3, .5)) < 1,
              "Triangle interior intersection missed");
        check(triangleDistanceSquared(V(-3, -3, 2), V(3, -3, 2), V(0, 3, 2)) > 1, "Separated face rejected");
        check(triangleDistanceSquared(V(-2, 0, 0), V(2, 0, 0), V(2, 0, 0)) == 0, "Degenerate edge intersection missed");
        check(triangleDistanceSquared(V(2, 0, 0), V(2, 0, 0), V(2, 0, 0)) == 4, "Degenerate point distance failed");
        check(triangleDistanceSquared(V(2, 0, 0), V(0, 2, 0), V(0, 0, 2)) > 1, "Outside triangle distance failed");
        std::cout << "PASS ellipsoid contact, rotations, scales, boundary triangles and degeneracies\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
