#ifndef model_LatticeBoxTree_H_
#define model_LatticeBoxTree_H_
#include <FMM/BoxTree.h>
#include <stdexcept>
namespace model
{
    // Complete the active periodic directions with orthogonal nonperiodic axes.
    // Bounds come from polyhedron vertices, not corners of its Cartesian AABB.
    template<int dim>
    std::shared_ptr<BoxTree<dim>> makeLatticeBoxTree(
        const Eigen::Matrix<double,dim,1>& lower,const Eigen::Matrix<double,dim,1>& upper,unsigned depth,
        const Eigen::Matrix<double,dim,Eigen::Dynamic>& periodic,
        const std::map<std::size_t,Eigen::Matrix<double,dim,1>>& points)
    {
        using Vector=Eigen::Matrix<double,dim,1>;
        using Matrix=Eigen::Matrix<double,dim,dim>;
        if(periodic.cols()==0) { return std::make_shared<BoxTree<dim>>(lower,upper,depth); }
        if(periodic.cols()>dim || !periodic.allFinite() || points.empty())
        { throw std::invalid_argument("Lattice tree requires valid periodic directions and polyhedron vertices"); }
        Matrix basis=Matrix::Zero(),orthogonal=Matrix::Zero();
        int rank=0;
        for(int c=0;c<periodic.cols();++c)
        {
            Vector v=periodic.col(c);
            for(int j=0;j<rank;++j) { v-=orthogonal.col(j).dot(v)*orthogonal.col(j); }
            if(v.norm()<=1.e-12*periodic.col(c).norm())
            { throw std::invalid_argument("Dependent periodic tree directions"); }
            orthogonal.col(rank)=v.normalized();basis.col(rank++)=periodic.col(c);
        }
        while(rank<dim)
        {
            Vector best=Vector::Zero();
            for(int c=0;c<dim;++c)
            {
                Vector v=Vector::Unit(c);
                for(int j=0;j<rank;++j) { v-=orthogonal.col(j).dot(v)*orthogonal.col(j); }
                if(v.squaredNorm()>best.squaredNorm()) { best=v; }
            }
            orthogonal.col(rank)=best.normalized();basis.col(rank)=best.normalized();++rank;
        }
        const Matrix inverse=basis.inverse();
        Vector lo=Vector::Constant(std::numeric_limits<double>::max()),hi=-lo;
        for(const auto& entry:points)
        {
            const Vector q=inverse*(entry.second-lower);
            lo=lo.cwiseMin(q);hi=hi.cwiseMax(q);
        }
        const Vector span=hi-lo;
        if(!span.allFinite() || (span.array()<=0).any())
        { throw std::invalid_argument("Degenerate lattice tree domain"); }
        // This is a covering hierarchy, not a periodic-image summation operator.
        const Vector origin=lower+basis*lo;
        const Matrix cell=basis*span.asDiagonal();
        return std::make_shared<BoxTree<dim>>(origin,cell,depth,typename BoxTree<dim>::LatticeTag{});
    }
}
#endif
