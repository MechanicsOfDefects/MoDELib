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
        bool advanced = false;

      public:
        void apply(vtkCamera *camera, const Eigen::Vector3d &center, const Eigen::Vector3d &axis,
                   double degreesPerFrame, size_t frameID, double advancePerFrame = 0.0)
        {
            if(!center.allFinite() || !axis.allFinite() || axis.stableNorm() == 0.0 || !std::isfinite(degreesPerFrame) || !std::isfinite(advancePerFrame))
            {
                throw std::invalid_argument("Enter a finite, nonzero rotation axis and finite angle and advance.");
            }
            if(lastCamera.GetPointer() != camera)
            {
                applied.setIdentity();
                lastCamera = camera;
                advanced = false;
            }
            // Reduce the absolute angle before conversion to radians. Sparse runIDs
            // and backward navigation therefore have the same orientation as direct jumps.
            const long double degrees = std::fmod(static_cast<long double>(degreesPerFrame) * frameID, 360.0L);
            const Eigen::Matrix3d target = Eigen::AngleAxisd(
                static_cast<double>(degrees) * std::acos(-1.0) / 180.0, axis.stableNormalized()).toRotationMatrix();
            if(degreesPerFrame == 0.0 && advancePerFrame == 0.0 && !advanced && target.isApprox(applied, 1.e-14))
            {
                return;
            }
            // Remove the previous automatic rotation before applying the new one.
            // Recenter the current view on the box while retaining its viewing
            // direction and zoom. Panning must not move the automatic rotation axis.
            const Eigen::Matrix3d change = target * applied.transpose();
            const Eigen::Vector3d focal(camera->GetFocalPoint());
            const Eigen::Vector3d offset = change * (Eigen::Vector3d(camera->GetPosition()) - focal);
            const double travel = advancePerFrame * static_cast<double>(frameID);
            if(!std::isfinite(travel))
            {
                throw std::invalid_argument("Camera advance times frameID must be finite.");
            }
            Eigen::Vector3d targetFocal = center;
            if(offset.stableNorm()>0.0)
            {
                targetFocal -= travel * offset.stableNormalized();
            }
            // Translate eye and focal point together so crossing the center neither
            // collapses the viewing direction nor turns the camera backward.
            const Eigen::Vector3d position = targetFocal + offset;
            if(advancePerFrame != 0.0)
            {
                camera->ParallelProjectionOff();
            }
            const Eigen::Vector3d up = change * Eigen::Vector3d(camera->GetViewUp());
            camera->SetFocalPoint(targetFocal.data());
            camera->SetPosition(position.data());
            camera->SetViewUp(up.data());
            camera->OrthogonalizeViewUp();
            applied = target;
            advanced = advancePerFrame != 0.0;
        }
    };
}
#endif
