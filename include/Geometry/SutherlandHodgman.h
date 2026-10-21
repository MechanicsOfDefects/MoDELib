/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SutherlandHodgman_H_
#define model_SutherlandHodgman_H_

#include <algorithm>
#include <vector>
#include <Eigen/Dense>

namespace model
{
    class SutherlandHodgman
    {
    public:
        // The clipping polygon is convex and counterclockwise. The subject
        // may be concave or self-intersecting; preserve its directed contour.
        template<typename PointType>
        constexpr static std::vector<PointType> clip(const std::vector<PointType>& subject,
                                                     const std::vector<PointType>& boundary)
        {
            if(boundary.empty() || subject.empty())
            {
                return {};
            }
            std::vector<PointType> ring=subject;
            PointType a=boundary.back();
            for(const auto& b : boundary)
            {
                if(ring.empty())
                {
                    return {};
                }
                const PointType edge=b-a;
                if(edge.squaredNorm()==0.0)
                {
                    a=b;
                    continue;
                }
                const auto distance=[&](const PointType& point)
                {
                    const PointType relative=point-a;
                    return edge(0)*relative(1)-edge(1)*relative(0);
                };
                std::vector<PointType> input;
                input.swap(ring);
                PointType start=input.back();
                double startDistance=distance(start);
                for(const auto& end : input)
                {
                    const double endDistance=distance(end);
                    const bool startInside=startDistance>=0.0;
                    const bool endInside=endDistance>=0.0;
                    if(startInside!=endInside)
                    {
                        // Use the same signed distances for classification and
                        // interpolation. Opposite signs guarantee t in [0,1],
                        // even for an edge nearly coincident with the boundary.
                        // Intersecting infinite lines using absolute-coordinate
                        // determinants can instead extrapolate far beyond it.
                        const double t=std::clamp(startDistance/(startDistance-endDistance),0.0,1.0);
                        ring.push_back(start+t*(end-start));
                    }
                    if(endInside)
                    {
                        ring.push_back(end);
                    }
                    start=end;
                    startDistance=endDistance;
                }
                a=b;
            }
            return ring;
        }
    };
}
#endif
