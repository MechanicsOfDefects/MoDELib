/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_FrameCameraRotation_H_
#define model_FrameCameraRotation_H_
#include <Eigen/Geometry>
#include <vtkCamera.h>
#include <vtkWeakPointer.h>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace model
{
    class FrameCameraRotation
    {
        Eigen::Matrix3d applied = Eigen::Matrix3d::Identity();
        vtkWeakPointer<vtkCamera> lastCamera;

      public:
        void apply(vtkCamera *camera, const Eigen::Vector3d &center, const Eigen::Vector3d &axis,
                   double degreesPerFrame, size_t frameID)
        {
            if(!center.allFinite() || !axis.allFinite() || axis.stableNorm() == 0.0 || !std::isfinite(degreesPerFrame))
            {
                throw std::invalid_argument("Enter a finite, nonzero rotation axis and a finite angle.");
            }
            if(lastCamera.GetPointer() != camera)
            {
                applied.setIdentity();
                lastCamera = camera;
            }
            // Reduce the absolute angle before conversion to radians. Sparse runIDs
            // and backward navigation therefore have the same orientation as direct jumps.
            const long double degrees = std::fmod(static_cast<long double>(degreesPerFrame) * frameID, 360.0L);
            const Eigen::Matrix3d target = Eigen::AngleAxisd(
                static_cast<double>(degrees) * std::acos(-1.0) / 180.0, axis.stableNormalized()).toRotationMatrix();
            if(degreesPerFrame == 0.0 && target.isApprox(applied, 1.e-14))
            {
                return;
            }
            // Remove the previous automatic rotation before applying the new one.
            // Recenter the current view on the box while retaining its viewing
            // direction and zoom. Panning must not move the automatic rotation axis.
            const Eigen::Matrix3d change = target * applied.transpose();
            const Eigen::Vector3d focal(camera->GetFocalPoint());
            const Eigen::Vector3d position = center + change * (Eigen::Vector3d(camera->GetPosition()) - focal);
            const Eigen::Vector3d up = change * Eigen::Vector3d(camera->GetViewUp());
            camera->SetFocalPoint(center.data());
            camera->SetPosition(position.data());
            camera->SetViewUp(up.data());
            camera->OrthogonalizeViewUp();
            applied = target;
        }
    };
}
#endif
