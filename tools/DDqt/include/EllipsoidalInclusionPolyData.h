/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_EllipsoidalInclusionPolyData_H_
#define model_EllipsoidalInclusionPolyData_H_

#include <Eigen/Dense>
#include <vtkAlgorithm.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

namespace model
{
    inline vtkSmartPointer<vtkPolyData> ellipsoidalInclusionPolyData(
        const Eigen::Vector3d& center, const Eigen::Matrix3d& axes)
    {
        vtkNew<vtkSphereSource> sphere;
        sphere->SetRadius(1.0);
        sphere->SetThetaResolution(32);
        sphere->SetPhiResolution(33);
        sphere->SetOutputPointsPrecision(vtkAlgorithm::DOUBLE_PRECISION);

        // Columns of axes are the semiaxis vectors, including their lengths.
        // Transforming a unit sphere directly avoids Euler-angle conventions.
        vtkNew<vtkMatrix4x4> matrix;
        matrix->Identity();
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                matrix->SetElement(i,j,axes(i,j));
            }
            matrix->SetElement(i,3,center(i));
        }
        vtkNew<vtkTransform> transform;
        transform->SetMatrix(matrix);
        vtkNew<vtkTransformPolyDataFilter> filter;
        filter->SetTransform(transform);
        filter->SetInputConnection(sphere->GetOutputPort());
        filter->SetOutputPointsPrecision(vtkAlgorithm::DOUBLE_PRECISION);
        filter->Update();
        auto result=vtkSmartPointer<vtkPolyData>::New();
        result->ShallowCopy(filter->GetOutput());
        return result;
    }
}
#endif
