/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <EllipsoidalInclusionsGenerator.h>
#include <EllipsoidalInclusionPlacement.h>
#include <MicrostructureGenerator.h>
#include <limits>
#include <random>

namespace model
{
    namespace
    {
        void validateEllipsoid(const EllipsoidalInclusionIndividualSpecification &spec, const double b)
        {
            using V = Eigen::Vector3d;
            using M = Eigen::Matrix3d;

            if(!std::isfinite(b) || b <= 0)
            {
                throw std::invalid_argument("EllipsoidalInclusionsGenerator: invalid Burgers-vector scale");
            }
            const M crystalAxes = spec.crystalAxes_SI / b;
            const V lengths = crystalAxes.colwise().stableNorm();
            const M crystalStrain = Eigen::Map<const M>(spec.crystalEigenDistortion.data()).transpose();
            if(!spec.center.allFinite() || !crystalAxes.allFinite() || !lengths.allFinite() ||
               lengths.minCoeff() <= 0 || lengths.minCoeff() / lengths.maxCoeff() < 1.e-6 ||
               !crystalStrain.allFinite() || !std::isfinite(spec.velocityReductionFactor))
            {
                throw std::invalid_argument("EllipsoidalInclusionsGenerator: invalid geometry or material input");
            }
            const M crystalRotation = crystalAxes * lengths.cwiseInverse().asDiagonal();
            if((crystalRotation.transpose() * crystalRotation - M::Identity()).cwiseAbs().maxCoeff() > 1.e-12)
            {
                throw std::invalid_argument("EllipsoidalInclusionsGenerator: semiaxes must be orthogonal");
            }
            const double strainScale = crystalStrain.cwiseAbs().maxCoeff();
            if((crystalStrain - crystalStrain.transpose()).cwiseAbs().maxCoeff() > 1.e-12 * strainScale)
            {
                throw std::invalid_argument("EllipsoidalInclusionsGenerator: eigen-distortion must be symmetric");
            }
        }
    } // namespace

    EllipsoidalInclusionsGenerator::EllipsoidalInclusionsGenerator(const EllipsoidalInclusionDensitySpecification &spec,
                                                                   MicrostructureGenerator &mg)
    {
        if(!std::isfinite(spec.targetDensity) || spec.targetDensity < 0.0)
        {
            throw std::invalid_argument("Ellipsoidal density must be finite and nonnegative");
        }
        if(spec.targetDensity == 0.0)
        {
            return;
        }
        if(!spec.axesLognormalDistributionMean_SI.allFinite() ||
           (spec.axesLognormalDistributionMean_SI.array() <= 0.0).any() ||
           !spec.axesLognormalDistributionStd_SI.allFinite() ||
           (spec.axesLognormalDistributionStd_SI.array() < 0.0).any() || !spec.patternVector_SI.allFinite() ||
           spec.maxAttempts == 0)
        {
            throw std::invalid_argument("Invalid ellipsoidal size distribution, pattern, or attempt limit");
        }
        if(!spec.crystalAxes.allFinite() ||
           (spec.crystalAxes.transpose() * spec.crystalAxes - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() >
               1.e-12)
        {
            throw std::invalid_argument("Ellipsoidal crystalAxes must contain orthogonal unit vectors");
        }
        Eigen::Vector3d logMean;
        Eigen::Vector3d logStd;
        for(int k = 0; k < 3; ++k)
        {
            const double mean = spec.axesLognormalDistributionMean_SI(k);
            const double stdSI = spec.axesLognormalDistributionStd_SI(k);
            // log(1 + (std/mean)^2), evaluated without squaring an extreme ratio.
            double logVariance = 0.0;
            if(stdSI > 0.0)
            {
                const double logRatio = std::log(stdSI) - std::log(mean);
                logVariance = logRatio > 0.0 ? 2.0 * logRatio + std::log1p(std::exp(-2.0 * logRatio))
                                             : std::log1p(std::exp(2.0 * logRatio));
            }
            logStd(k) = std::sqrt(logVariance);
            logMean(k) = std::log(mean) - 0.5 * logVariance;
        }
        const double b = mg.ddBase.poly.b_SI;
        EllipsoidalInclusionIndividualSpecification candidate;
        candidate.crystalAxes_SI = spec.crystalAxes * spec.axesLognormalDistributionMean_SI.asDiagonal();
        candidate.crystalEigenDistortion = spec.crystalEigenDistortion;
        candidate.allowOverlap = spec.allowOverlap;
        candidate.allowOutside = spec.allowOutside;
        candidate.velocityReductionFactor = spec.velocityReductionFactor;
        candidate.phaseID = spec.phaseID;
        validateEllipsoid(candidate, b);
        const double volumeSI = mg.ddBase.mesh.volume() * std::pow(b, 3);
        const double requested = std::ceil(spec.targetDensity * volumeSI);
        if(!std::isfinite(volumeSI) || volumeSI <= 0 || !std::isfinite(requested) ||
           requested >= static_cast<double>(std::numeric_limits<size_t>::max()))
        {
            throw std::invalid_argument("Ellipsoidal density exceeds numerical range");
        }
        const size_t targetCount = static_cast<size_t>(requested);
        const Eigen::Vector3d pattern = spec.patternVector_SI.transpose() / b;
        const double spacing = pattern.stableNorm();
        if(!pattern.allFinite() || !std::isfinite(spacing))
        {
            throw std::invalid_argument("Ellipsoidal pattern exceeds numerical range");
        }
        std::mt19937 random;
        std::normal_distribution<double> normal(0.0, 1.0);
        size_t generated = 0;
        size_t attempts = 0;
        while(generated < targetCount)
        {
            if(attempts == spec.maxAttempts)
            {
                throw std::runtime_error("Ellipsoidal density placement exhausted maxAttempts after generating " +
                                         std::to_string(generated) + " of " + std::to_string(targetCount));
            }
            ++attempts;
            const auto point = mg.ddBase.poly.randomLatticePointInMesh();
            candidate.center = point.first.cartesian();
            if(spacing > 0.0)
            {
                const Eigen::Vector3d normalGlobal = mg.ddBase.poly.grain(point.second)->C2G * (pattern / spacing);
                const double distance = candidate.center.dot(normalGlobal);
                candidate.center += (std::round(distance / spacing) * spacing - distance) * normalGlobal;
                // Projection must not change the grain defining the plane orientation.
                const auto projected = mg.ddBase.mesh.search(candidate.center);
                if(!projected.first || projected.second->region->regionID != static_cast<size_t>(point.second))
                {
                    continue;
                }
            }
            Eigen::Vector3d sampledLengths;
            for(int k = 0; k < 3; ++k)
            {
                if(spec.axesLognormalDistributionStd_SI(k) == 0.0)
                {
                    sampledLengths(k) = spec.axesLognormalDistributionMean_SI(k);
                }
                else
                {
                    sampledLengths(k) = std::exp(logMean(k) + logStd(k) * normal(random));
                }
            }
            candidate.crystalAxes_SI = spec.crystalAxes * sampledLengths.asDiagonal();
            const Eigen::Vector3d lengths = candidate.crystalAxes_SI.colwise().stableNorm();
            if(!candidate.crystalAxes_SI.allFinite() || !lengths.allFinite() || lengths.minCoeff() <= 0 ||
               lengths.minCoeff() / lengths.maxCoeff() < 1.e-6)
            {
                continue;
            }
            const auto before = mg.config().ellipsoidalInclusions().size();
            EllipsoidalInclusionsGenerator single(candidate, mg);
            if(mg.config().ellipsoidalInclusions().size() > before)
            {
                ++generated;
                attempts = 0;
            }
        }
        std::cout << "Generated ellipsoidal number density=" << generated / volumeSI << " m^-3\n";
    }

    EllipsoidalInclusionsGenerator::EllipsoidalInclusionsGenerator(
        const EllipsoidalInclusionIndividualSpecification &spec, MicrostructureGenerator &mg)
    {
        using V = Eigen::Vector3d;
        using M = Eigen::Matrix3d;
        const double b = mg.ddBase.poly.b_SI;
        validateEllipsoid(spec, b);
        const M crystalAxes = spec.crystalAxes_SI / b;
        const V lengths = crystalAxes.colwise().stableNorm();
        const M crystalRotation = crystalAxes * lengths.cwiseInverse().asDiagonal();
        const M crystalStrain = Eigen::Map<const M>(spec.crystalEigenDistortion.data()).transpose();
        const auto searchPair = mg.ddBase.mesh.search(spec.center);
        if(!searchPair.first)
        {
            std::cout << "Ellipsoidal inclusion center outside mesh; cannot determine crystal frame; not generated\n";
            return;
        }
        const auto &grain = mg.ddBase.poly.grain(searchPair.second->region->regionID);
        const M strain = grain->C2G * crystalStrain * grain->C2G.transpose();
        const M axes = grain->C2G * crystalAxes;
        const M rotation = grain->C2G * crystalRotation;
        const M inverse = lengths.cwiseInverse().asDiagonal() * rotation.transpose();
        const M shape = axes * axes.transpose();
        if(!shape.allFinite() || !inverse.allFinite() || shape.diagonal().minCoeff() <= 0)
        {
            throw std::invalid_argument("EllipsoidalInclusionsGenerator: geometry exceeds numerical range");
        }

        if(!spec.allowOutside)
        {
            // A connected ellipsoid whose center is inside remains inside unless
            // it reaches a boundary triangle. This also handles nonconvex meshes.
            for(const auto &face : mg.ddBase.mesh.template observer<2>())
            {
                if(face.second->isBoundarySimplex())
                {
                    const auto vertices = face.second->vertices();
                    if(ellipsoidalPlacement::triangleDistanceSquared(
                           inverse * (vertices[0]->P0 - spec.center), inverse * (vertices[1]->P0 - spec.center),
                           inverse * (vertices[2]->P0 - spec.center)) <= 1.0 + 1.e-12)
                    {
                        std::cout << "Ellipsoidal inclusion reaches mesh boundary; not generated\n";
                        return;
                    }
                }
            }
        }
        if(!spec.allowOverlap)
        {
            for(const auto &inclusion : mg.config().sphericalInclusions())
            {
                if(ellipsoidalPlacement::overlap(inclusion.C - spec.center, shape,
                                                 M::Identity() * inclusion.a * inclusion.a))
                {
                    std::cout << "Ellipsoidal inclusion overlaps a sphere; not generated\n";
                    return;
                }
            }
            for(const auto &inclusion : mg.config().ellipsoidalInclusions())
            {
                if(ellipsoidalPlacement::overlap(inclusion.C - spec.center, shape,
                                                 inclusion.axes * inclusion.axes.transpose()))
                {
                    std::cout << "Ellipsoidal inclusion overlaps an ellipsoid; not generated\n";
                    return;
                }
            }
            // Conservative test for arbitrary polyhedra: disjoint bounding boxes
            // certify separation. Intersecting boxes may reject a valid placement.
            const V extent = shape.diagonal().cwiseSqrt();
            for(const auto &inclusion : mg.config().polyhedronInclusions())
            {
                V lower = V::Constant(std::numeric_limits<double>::infinity()), upper = -lower;
                bool found = false;
                for(const auto &edge : mg.config().polyhedronInclusionEdges())
                {
                    if(edge.inclusionID != inclusion.inclusionID)
                    {
                        continue;
                    }
                    for(const auto nodeID : {edge.sourceID, edge.sinkID})
                    {
                        const auto &nodes = mg.config().polyhedronInclusionNodes();
                        const auto node =
                            std::find_if(nodes.begin(), nodes.end(), [&](const auto &n) { return n.nodeID == nodeID; });
                        if(node == nodes.end() || !node->P.allFinite())
                        {
                            throw std::runtime_error("EllipsoidalInclusionsGenerator: invalid polyhedron node");
                        }
                        lower = lower.cwiseMin(node->P);
                        upper = upper.cwiseMax(node->P);
                        found = true;
                    }
                }
                if(!found)
                {
                    throw std::runtime_error("EllipsoidalInclusionsGenerator: polyhedron has no boundary edges");
                }
                if(((spec.center + extent).array() >= lower.array()).all() &&
                   ((spec.center - extent).array() <= upper.array()).all())
                {
                    std::cout << "Ellipsoidal inclusion intersects polyhedron bounding box; not generated\n";
                    return;
                }
            }
        }
        mg.insertInclusion(spec.center, axes, strain, spec.velocityReductionFactor, spec.phaseID);
        std::cout << "Generated ellipsoidal inclusion\n";
    }
} // namespace model
