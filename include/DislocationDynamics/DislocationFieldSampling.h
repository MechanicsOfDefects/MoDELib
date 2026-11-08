/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_DislocationFieldSampling_H_
#define model_DislocationFieldSampling_H_
#include <DislocationDynamics/DislocationFieldBuffers.h>
#include <DislocationDynamics/DislocationSegment.h>
#include <DislocationDynamics/DislocationQuadraturePoint.h>

namespace model
{
    // Explicit sampling adapter: a caller iterates unique network links and owns
    // the snapshot. It does not depend on the segment's mobility quadrature cache.
    template<int dim>
    struct DislocationFieldSampling
    {
        using Geometry=DislocationQuadratureGeometry<dim>;
        using Segment=DislocationSegment<dim>;
        using Point=DislocationQuadraturePoint<dim>;
        using Rule=typename Point::QuadratureDynamicType;
        using PowerRule=typename Point::QuadPowDynamicType;

        static std::vector<Geometry> sample(const Segment& segment, double pointsPerLength)
        {
            if(!std::isfinite(pointsPerLength) || pointsPerLength<=0.0)
            {
                throw std::invalid_argument("Source/target quadrature density must be positive and finite");
            }
            const double requested=pointsPerLength*segment.chord().norm();
            if(!std::isfinite(requested) || requested>1024.0)
            {
                throw std::invalid_argument("Requested dislocation quadrature exceeds supported order");
            }
            const int order=PowerRule::lowerOrder(requested);
            const auto coefficients=segment.sfCoeffs();
            const auto dofs=segment.hermiteDofs();
            std::vector<Geometry> points;
            points.reserve(order);
            for(int q=0;q<order;++q)
            {
                const typename Point::VectorDim position=(PowerRule::uPow(order).row(q)*coefficients*dofs).transpose();
                const typename Point::VectorDim derivative=(PowerRule::duPow(order).row(q)*
                    coefficients.template block<Point::Ncoeff-1,Point::Ncoeff>(1,0)*dofs).transpose();
                const double jacobian=derivative.norm();
                if(!std::isfinite(jacobian) || jacobian<=0.0)
                {
                    throw std::invalid_argument("Degenerate segment quadrature geometry");
                }
                Geometry point;
                point.parameter=Rule::abscissa(order,q);
                point.lineWeight=Rule::weight(order,q)*jacobian;
                for(int axis=0;axis<dim;++axis)
                {
                    point.position[axis]=position[axis];
                    point.tangent[axis]=derivative[axis]/jacobian;
                }
                points.push_back(point);
            }
            return points;
        }

        static void append(DislocationFieldBuffers<dim>& buffers, const Segment& segment,
                           double sourcePointsPerLength, double targetPointsPerLength,
                           DislocationTargetPurpose purpose)
        {
            const DislocationSegmentFieldSelection selection{segment.isGlissile(),
                segment.isBoundarySegment(),segment.isGrainBoundarySegment()};
            typename DislocationFieldBuffers<dim>::Vector burgers{};
            std::vector<Geometry> sources;
            std::vector<Geometry> targets;
            if(!segment.hasZeroBurgers())
            {
                for(int axis=0;axis<dim;++axis)
                {
                    burgers[axis]=segment.burgers()[axis];
                }
                sources=sample(segment,sourcePointsPerLength);
                if(selection.accepts(purpose))
                {
                    targets=sample(segment,targetPointsPerLength);
                }
            }
            buffers.addSegment({segment.source->sID,segment.sink->sID},burgers,
                               selection,sources,targets,purpose);
        }
    };
}
#endif
