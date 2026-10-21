#ifndef model_PeriodicMeshFacePairing_h_
#define model_PeriodicMeshFacePairing_h_

#include <Polyhedron.h>
#include <PeriodicSimplexPair.h>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace model
{
    // Compare the complete boundary discretization, not just its convex hull.
    template<int dim>
    bool isPeriodicFacePair(const Polyhedron<dim,dim-1>& a, const Polyhedron<dim,dim-1>& b,
                            Eigen::Matrix<double,dim,1>& shift, double tolerance,
                            PeriodicSimplexPairMap<dim>* correspondence=nullptr)
    {
        using Vertex = Simplex<dim,0>;
        if(&a == &b || !a.isExternal() || !b.isExternal() || a.simplices().size() != b.simplices().size()
           || (a.outNormal()+b.outNormal()).norm() > 1.e-6)
        {
            return false;
        }
        std::set<const Vertex*> av, bv;
        for(const auto* t : a.simplices())
        {
            for(const auto* v : t->vertices())
            {
                av.insert(v);
            }
        }
        for(const auto* t : b.simplices())
        {
            for(const auto* v : t->vertices())
            {
                bv.insert(v);
            }
        }
        if(av.empty() || av.size() != bv.size())
        {
            return false;
        }
        Eigen::Matrix<double,dim,1> ca = Eigen::Matrix<double,dim,1>::Zero();
        Eigen::Matrix<double,dim,1> cb = ca;
        for(const auto* v : av)
        {
            ca += (v->P0-(*av.begin())->P0)/double(av.size());
        }
        for(const auto* v : bv)
        {
            cb += (v->P0-(*bv.begin())->P0)/double(bv.size());
        }
        shift = ((*bv.begin())->P0-(*av.begin())->P0)+cb-ca;
        if(shift.norm() <= tolerance || std::abs(a.outNormal().dot(shift.normalized())) <= 1.e-6)
        {
            return false;
        }
        PeriodicSimplexPairMap<dim> matched;
        auto remaining = b.simplices();
        for(const auto* ta : a.simplices())
        {
            auto match = remaining.end();
            for(auto it = remaining.begin(); it != remaining.end(); ++it)
            {
                const auto* tb = *it;
                if((ta->outNormal()+tb->outNormal()).norm() > 1.e-6
                   || (ta->center()+shift-tb->center()).norm() > tolerance)
                {
                    continue;
                }
                std::set<const Vertex*> vertices;
                for(const auto* v : tb->vertices())
                {
                    vertices.insert(v);
                }
                for(const auto* va : ta->vertices())
                {
                    auto found = std::find_if(vertices.begin(), vertices.end(), [&](const Vertex* vb)
                    {
                        return (va->P0+shift-vb->P0).norm() <= tolerance;
                    });
                    if(found == vertices.end())
                    {
                        break;
                    }
                    vertices.erase(found);
                }
                if(vertices.empty())
                {
                    if(match != remaining.end())
                    {
                        return false; // Ambiguous triangle match.
                    }
                    match = it;
                }
            }
            if(match == remaining.end())
            {
                return false;
            }
            if(correspondence)
            {
                PeriodicSimplexPair<dim> record{*match,shift,{}};
                const auto va=ta->vertices();
                const auto vb=(*match)->vertices();
                std::set<size_t> used;
                for(size_t i=0;i<dim;++i)
                {
                    size_t count=0;
                    for(size_t j=0;j<dim;++j)
                    {
                        if((va[i]->P0+shift-vb[j]->P0).norm()<=tolerance)
                        {
                            record.vertexPermutation[i]=j;
                            ++count;
                        }
                    }
                    if(count!=1 || !used.insert(record.vertexPermutation[i]).second)
                    {
                        return false;
                    }
                }
                matched.emplace(ta,record);
            }
            remaining.erase(match);
        }
        if(!remaining.empty())
        {
            return false;
        }
        if(correspondence)
        {
            *correspondence=std::move(matched);
        }
        return true;
    }

    template<int dim>
    struct PeriodicMeshLattice
    {
        using Face = Polyhedron<dim,dim-1>;
        using Vector = Eigen::Matrix<double,dim,1>;
        using Basis = Eigen::Matrix<double,dim,Eigen::Dynamic>;
        Basis basis = Basis(dim,0);
        std::map<Face*,std::pair<Vector,Face*>> pairs;
    };

    template<int dim>
    bool periodicLatticeContains(const Eigen::Matrix<double,dim,Eigen::Dynamic>& basis,
                                 const Eigen::Matrix<double,dim,1>& shift, double tolerance)
    {
        if(basis.cols() == 0)
        {
            return shift.norm() <= tolerance;
        }
        const auto coefficients = basis.fullPivLu().solve(shift).eval();
        // Check reconstruction as well as integer coefficients: for partial
        // periodicity the translation must also lie in the basis span.
        return (basis*coefficients.array().round().matrix()-shift).norm() <= tolerance;
    }

    template<int dim>
    PeriodicMeshLattice<dim> discoverPeriodicMeshLattice(
        const std::vector<Polyhedron<dim,dim-1>*>& faces, double tolerance)
    {
        using Vector = Eigen::Matrix<double,dim,1>;
        using Lattice = PeriodicMeshLattice<dim>;
        using Basis = typename Lattice::Basis;
        struct Candidate
        {
            size_t a, b;
            Vector shift;
        };
        Lattice best;
        std::vector<Candidate> candidates;
        std::vector<Vector> generators;
        for(size_t a=0; a<faces.size(); ++a)
        {
            for(size_t b=a+1; b<faces.size(); ++b)
            {
                Vector shift;
                if(!isPeriodicFacePair(*faces[a],*faces[b],shift,tolerance))
                {
                    continue;
                }
                candidates.push_back({a,b,shift});
            }
        }
        // Restore the original disambiguation using the complete, unfiltered
        // candidate sets. Unique matches are retained, including XY diagonals.
        std::vector<std::vector<size_t>> incident(faces.size());
        for(size_t k=0; k<candidates.size(); ++k)
        {
            incident[candidates[k].a].push_back(k);
            incident[candidates[k].b].push_back(k);
        }
        std::vector<std::set<size_t>> allowed(faces.size());
        for(size_t a=0; a<faces.size(); ++a)
        {
            for(const size_t k : incident[a])
            {
                const Vector& shift = candidates[k].shift;
                bool compatible = true;
                if(incident[a].size() > 1)
                {
                    for(size_t other=0; other<faces.size(); ++other)
                    {
                        if(other == a || incident[other].empty()
                           || std::abs(faces[other]->outNormal().dot(shift.normalized())) <= 1.e-6)
                        {
                            continue;
                        }
                        const bool sameShift = std::any_of(incident[other].begin(),incident[other].end(),[&](size_t j)
                        {
                            return (candidates[j].shift-shift).norm() <= tolerance
                                || (candidates[j].shift+shift).norm() <= tolerance;
                        });
                        if(!sameShift)
                        {
                            compatible = false;
                            break;
                        }
                    }
                }
                if(compatible)
                {
                    allowed[a].insert(k);
                }
            }
            if(allowed[a].size() > 1)
            {
                throw std::runtime_error("Multiple compatible periodic partners for face " + std::to_string(faces[a]->sID));
            }
        }
        std::vector<Candidate> unambiguous;
        for(size_t k=0; k<candidates.size(); ++k)
        {
            const auto& pair = candidates[k];
            const bool forward = allowed[pair.a].count(k);
            const bool reverse = allowed[pair.b].count(k);
            if(forward != reverse)
            {
                throw std::runtime_error("Periodic face disambiguation is not reciprocal");
            }
            if(forward && reverse)
            {
                unambiguous.push_back(pair);
            }
        }
        candidates.swap(unambiguous);
        for(const auto& pair : candidates)
        {
            Vector shift = pair.shift;
            // Fix the sign for repeatable basis selection.
            for(int k=0; k<dim; ++k)
            {
                if(std::abs(shift[k]) > tolerance)
                {
                    if(shift[k] < 0.0)
                    {
                        shift = -shift;
                    }
                    break;
                }
            }
            bool duplicate = false;
            for(const auto& previous : generators)
            {
                if((shift-previous).norm() <= tolerance || (shift+previous).norm() <= tolerance)
                {
                    duplicate = true;
                    break;
                }
            }
            if(!duplicate)
            {
                generators.push_back(shift);
            }
        }
        // Prefer short independent generators, but test alternatives: independence
        // alone does not ensure that diagonal translations have integer coordinates.
        std::sort(generators.begin(),generators.end(),[](const Vector& a,const Vector& b)
        {
            return a.squaredNorm() < b.squaredNorm();
        });
        double bestScore = std::numeric_limits<double>::infinity();
        auto validate = [&](const Basis& basis)
        {
            Lattice lattice;
            lattice.basis = basis;
            for(const auto& pair : candidates)
            {
                if(periodicLatticeContains<dim>(basis,pair.shift,tolerance))
                {
                    if(!lattice.pairs.emplace(faces[pair.a],std::make_pair(pair.shift,faces[pair.b])).second
                       || !lattice.pairs.emplace(faces[pair.b],std::make_pair(Vector(-pair.shift),faces[pair.a])).second)
                    {
                        return; // More than one admissible partner for a face.
                    }
                }
            }
            for(auto* face : faces)
            {
                bool required = false;
                for(int k=0; k<basis.cols(); ++k)
                {
                    required = required || std::abs(face->outNormal().dot(basis.col(k).normalized())) > 1.e-6;
                }
                if(required && !lattice.pairs.count(face))
                {
                    return;
                }
            }
            const double score = basis.squaredNorm();
            if(basis.cols() > best.basis.cols() || (basis.cols() == best.basis.cols() && score < bestScore))
            {
                best = std::move(lattice);
                bestScore = score;
            }
        };
        std::function<void(size_t,const Basis&)> search;
        search = [&](size_t begin,const Basis& basis)
        {
            for(size_t j=begin; j<generators.size(); ++j)
            {
                Basis extended(dim,basis.cols()+1);
                extended.leftCols(basis.cols()) = basis;
                extended.col(basis.cols()) = generators[j];
                Basis normalized = extended.colwise().normalized();
                Eigen::FullPivLU<Basis> lu(normalized);
                lu.setThreshold(1.e-6);
                if(lu.rank() != extended.cols())
                {
                    continue;
                }
                validate(extended);
                // A rejected one-direction basis may become valid jointly with
                // another direction. Never prune based on per-vector coverage.
                if(extended.cols() < dim)
                {
                    search(j+1,extended);
                }
            }
        };
        search(0,Basis(dim,0));
        if(best.basis.cols() == 0)
        {
            throw std::runtime_error("No globally valid periodic lattice found");
        }
        return best;
    }

    // Selection uses partner translations from the validated full lattice. Face
    // normals and boundary-coverage requirements apply only during discovery.
    template<int dim>
    PeriodicMeshLattice<dim> selectPeriodicMeshLattice(
        const PeriodicMeshLattice<dim>& full, const std::set<int>& seeds, double tolerance)
    {
        using Lattice = PeriodicMeshLattice<dim>;
        using Basis = typename Lattice::Basis;
        using Vector = typename Lattice::Vector;
        if(seeds.count(-1))
        {
            return full;
        }
        Lattice active;
        if(seeds.empty())
        {
            return active;
        }
        std::vector<Vector> generators;
        for(const int id : seeds)
        {
            const auto pair = std::find_if(full.pairs.begin(),full.pairs.end(),[&](const auto& item)
            {
                return item.first->sID == static_cast<size_t>(id);
            });
            if(pair == full.pairs.end())
            {
                throw std::runtime_error("Face " + std::to_string(id) + " has no partner in the full periodic lattice");
            }
            const Vector shift = pair->second.first;
            const bool duplicate = std::any_of(generators.begin(),generators.end(),[&](const Vector& previous)
            {
                return (shift-previous).norm() <= tolerance || (shift+previous).norm() <= tolerance;
            });
            if(!duplicate)
            {
                generators.push_back(shift);
            }
        }
        double bestScore = std::numeric_limits<double>::infinity();
        std::function<void(size_t,const Basis&)> search;
        search = [&](size_t begin,const Basis& basis)
        {
            for(size_t j=begin; j<generators.size(); ++j)
            {
                Basis extended(dim,basis.cols()+1);
                extended.leftCols(basis.cols()) = basis;
                extended.col(basis.cols()) = generators[j];
                Eigen::FullPivLU<Basis> lu(extended.colwise().normalized());
                lu.setThreshold(1.e-6);
                if(lu.rank() != extended.cols())
                {
                    continue;
                }
                const bool representsAll = std::all_of(generators.begin(),generators.end(),[&](const Vector& shift)
                {
                    return periodicLatticeContains<dim>(extended,shift,tolerance);
                });
                if(representsAll && extended.squaredNorm() < bestScore)
                {
                    active.basis = extended;
                    bestScore = extended.squaredNorm();
                }
                if(extended.cols() < dim)
                {
                    search(j+1,extended);
                }
            }
        };
        search(0,Basis(dim,0));
        if(active.basis.cols() == 0)
        {
            throw std::runtime_error("Selected translations require a lattice basis that cannot be formed from the selected face translations");
        }
        for(const auto& pair : full.pairs)
        {
            if(periodicLatticeContains<dim>(active.basis,pair.second.first,tolerance))
            {
                active.pairs.emplace(pair);
            }
        }
        return active;
    }

}
#endif
