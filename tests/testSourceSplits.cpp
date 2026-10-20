#include <SphericalInclusionIO.h>
#include <PolyhedronInclusionIO.h>
#include <PolyhedronInclusionNodeIO.h>
#include <PolyhedronInclusionEdgeIO.h>
#include <EllipsoidalInclusionIO.h>
#include <SphericalInclusion.h>
#include <PolyhedronInclusion.h>
#include <EllipsoidalInclusion.h>
#include <BestRationalApproximation.h>
#include <DiophantineSolver.h>
#include <Rational.h>
#include <SchurComplementSolver.h>
#include <CatmullRomSplineSegment.h>
#include <StrUtilities.h>
#include <EnvVariable.h>
#include <GaussLegendreRecursiveRelation.h>
#include <SimplexReader.h>
#include <type_traits>
#include <stdexcept>
#include <iostream>

namespace
{
    struct Unrelated {};
    void check(bool ok, const char* message)
    {
        if(!ok)
        {
            throw std::runtime_error(message);
        }
    }
    template<class Record>
    void roundTrip(const Record& original)
    {
        std::stringstream first;
        first << original;
        Record restored(first);
        check(!first.fail(), "Record parsing failed");
        std::stringstream second;
        second << restored;
        check(first.str()==second.str(), "Record serialization changed");
    }
}

int main()
{
    using namespace model;
    static_assert(std::is_constructible_v<SphericalInclusionIO<3>,const SphericalInclusion<3>&>);
    static_assert(std::is_constructible_v<PolyhedronInclusionIO<3>,const PolyhedronInclusion<3>&>);
    static_assert(std::is_constructible_v<EllipsoidalInclusionIO<3>,const EllipsoidalInclusion<3>&>);
    static_assert(!std::is_constructible_v<SphericalInclusionIO<3>,const Unrelated&>);
    static_assert(!std::is_constructible_v<PolyhedronInclusionIO<3>,const Unrelated&>);
    static_assert(!std::is_constructible_v<EllipsoidalInclusionIO<3>,const Unrelated&>);
    try
    {
        Eigen::Matrix3d e;
        e << 1,2,3,4,5,6,7,8,9;
        const Eigen::Vector3d center(1,-2,3);
        const SphericalInclusion<3> sphere(center,2.5,e,0.3,1.0,0.75,2,nullptr);
        const SphericalInclusionIO<3> sphereIO(sphere);
        check(sphereIO.C==center && sphereIO.a==2.5 && sphereIO.eT==e && sphereIO.inclusionID==sphere.sID,
              "Sphere conversion changed fields");
        roundTrip(sphereIO);
        const std::map<size_t,PolyhedronInclusionNodeIO<3>> nodes;
        const std::map<size_t,std::vector<size_t>> faces;
        const PolyhedronInclusion<3> poly(nodes,faces,e,0.3,1.0,0.6,4,nullptr);
        const PolyhedronInclusionIO<3> polyIO(poly);
        check(polyIO.eT==e && polyIO.inclusionID==poly.sID && polyIO.phaseID==4 && polyIO.mobilityReduction==0.6,
              "Polyhedron conversion changed fields");
        roundTrip(polyIO);
        roundTrip(PolyhedronInclusionNodeIO<3>(9,center));
        roundTrip(PolyhedronInclusionEdgeIO(3,4,5,6));
        std::stringstream sphereText("7 1 -2 3 2.5 1 2 3 4 5 6 7 8 9 0.75 2");
        const SphericalInclusionIO<3> parsed(sphereText);
        check(parsed.inclusionID==7 && parsed.C==center && parsed.eT==e && parsed.phaseID==2,
              "Sphere record field order changed");
        const BestRationalApproximation best(0.125,100);
        check(best.num==1 && best.den==8,"Rational approximation changed");
        const Rational sum=Rational(1,3)+Rational(1,6);
        check(sum.n==1 && sum.d==2 && sum.asDouble()==0.5,"Rational arithmetic changed");
        long x=0,y=0;
        DiophantineSolver::solveDiophantine2vars(6,9,3,x,y);
        check(6*x+9*y==3,"Diophantine solution changed");
        const Eigen::Matrix2d A=Eigen::Matrix2d::Identity()*2;
        Eigen::Matrix<double,1,2> B;
        B << 1,1;
        const Eigen::Vector2d f(3,5);
        const Eigen::VectorXd g=Eigen::VectorXd::Constant(1,1);
        const SchurComplementSolver solve(A,B,f,g);
        check((A*solve.X+B.transpose()*solve.L-f).norm()<1e-12 && (B*solve.X-g).norm()<1e-12,
              "Schur solve changed");
        const CatmullRomSplineSegment<3> segment(Eigen::Vector3d(0,0,0),Eigen::Vector3d(1,0,0),
                                                Eigen::Vector3d(2,0,0),Eigen::Vector3d(3,0,0));
        check(segment.get_r(0.25).isApprox(Eigen::Vector3d(1.25,0,0)) &&
              segment.get_ru(0.5).isApprox(Eigen::Vector3d(1,0,0)) && segment.get_ruu(0.5).norm()<1e-12,
              "Spline evaluation changed");
        check(StrUtilities::lowercase("AbC")=="abc" && StrUtilities::removeSpaces(" a\t b ")=="ab",
              "String utility changed");
        check(!EnvVariable("PATH").empty(),"Environment variable lookup failed");
        check(GaussLegendreRecursiveRelation::a(2)==1.5 && GaussLegendreRecursiveRelation::b(2)==0 &&
              GaussLegendreRecursiveRelation::c(2)==0.5,"Quadrature recurrence changed");
        SimplexReaderBase<3> reader;
        reader.nodes().emplace(9,center);
        check(reader.simplexReader().nodes().at(9)==center,"Mesh reader accessor changed");
        reader.clear();
        check(reader.nodes().empty() && reader.elements().empty(),"Mesh reader clear failed");
        std::cout << "PASS extracted APIs, concrete inclusion constructors, serialization and numerics\n";
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
