#include <EllipsoidalInclusionPolyData.h>
#include <Eigen/Geometry>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <iostream>
#include <limits>
#include <stdexcept>

void check(bool condition, const char *message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

void verify(const Eigen::Vector3d &center, const Eigen::Matrix3d &axes)
{
    const auto mesh = model::ellipsoidalInclusionPolyData(center, axes);
    check(mesh->GetNumberOfPoints() > 0 && mesh->GetNumberOfCells() > 0, "Empty ellipsoid geometry");
    const Eigen::Matrix3d inverse = axes.inverse();
    auto *normals = mesh->GetPointData()->GetNormals();
    check(normals != nullptr, "Missing surface normals");
    Eigen::Matrix<double, 3, 2> tipDistance =
        Eigen::Matrix<double, 3, 2>::Constant(std::numeric_limits<double>::infinity());
    for(vtkIdType k = 0; k < mesh->GetNumberOfPoints(); ++k)
    {
        double p[3], n[3];
        mesh->GetPoint(k, p);
        normals->GetTuple(k, n);
        const Eigen::Vector3d position(p[0], p[1], p[2]);
        const Eigen::Vector3d unit = inverse * (position - center);
        check(std::abs(unit.squaredNorm() - 1.0) < 1.e-11, "Incorrect surface scale or orientation");
        const Eigen::Vector3d normal(n[0], n[1], n[2]);
        const Eigen::Vector3d expectedNormal = (inverse.transpose() * unit).normalized();
        check(normal.normalized().dot(expectedNormal) > 1.0 - 1.e-6, "Incorrect transformed surface normal");
        for(int j = 0; j < 3; ++j)
        {
            tipDistance(j, 0) = std::min(tipDistance(j, 0), (position - center - axes.col(j)).norm());
            tipDistance(j, 1) = std::min(tipDistance(j, 1), (position - center + axes.col(j)).norm());
        }
    }
    check(tipDistance.maxCoeff() < 1.e-10, "Semiaxis endpoints missing or misoriented");
}

int main()
{
    try
    {
        const Eigen::Vector3d center(4, -3, 2);
        verify(center, Eigen::Vector3d(3, 2, 1).asDiagonal());
        const Eigen::Matrix3d rotation =
            Eigen::AngleAxisd(.73, Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
        Eigen::Matrix3d axes = rotation * Eigen::Vector3d(7, 2, .5).asDiagonal();
        verify(center, axes);
        axes.col(1) *= -1;
        verify(center, axes);
        axes.col(0).swap(axes.col(2));
        verify(center, axes);
        verify(center, 2.0 * Eigen::Matrix3d::Identity());
        std::cout
            << "PASS rendered ellipsoid surfaces, semiaxis endpoints, rotation, reflection, permutation and normals\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
