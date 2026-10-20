/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexReader_H_
#define model_SimplexReader_H_

#include <map>
#include <iostream>
#include <Eigen/Dense>
#include <SimplexTraits.h>
// #include <VertexReader.h>
#include <SimplexTraits.h>

namespace model
{

    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    struct SimplexReaderBase
        : private std::map<size_t, Eigen::Matrix<double, dim, 1>>
        /*                      */,
          private std::map<size_t, std::pair<typename SimplexTraits<dim, dim>::SimplexIDType, size_t>>
    {

        typedef std::map<size_t, Eigen::Matrix<double, dim, 1>> NodeContainerType;
        typedef std::map<size_t, std::pair<typename SimplexTraits<dim, dim>::SimplexIDType, size_t>>
            ElementContainerType;

        /**********************************************************************/
        const NodeContainerType &nodes() const;

        /**********************************************************************/
        NodeContainerType &nodes();

        /**********************************************************************/
        const ElementContainerType &elements() const;

        /**********************************************************************/
        ElementContainerType &elements();

        /**********************************************************************/
        void clear();

        SimplexReaderBase<dim> &simplexReader();

        const SimplexReaderBase<dim> &simplexReader() const;

        //        static IDreader<'N',1,dim,double> nodeReader;

        /**********************************************************************/
        const Eigen::Matrix<double, dim, 1> &get_P0(const typename SimplexTraits<dim, 0>::SimplexIDType &xID);
    };

    /**************************************************************************/
    /**************************************************************************/
    template <int dim> struct SimplexReader : public SimplexReaderBase<dim>
    {

        /**********************************************************************/
        void read(const int &meshID);
    };

    template <> struct SimplexReader<3> : public SimplexReaderBase<3>
    {

        static constexpr int dim = 3;

        /**********************************************************************/
        void read(const std::string &meshFileName, const Eigen::Matrix<double, 3, 3> &A,
                  const Eigen::Matrix<double, 3, 1> &x0);
    };

} // namespace model
#endif
