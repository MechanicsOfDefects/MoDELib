#ifndef model_ParallelFor_H_
#define model_ParallelFor_H_
#include <cstddef>
#include <algorithm>
#include <cstdlib>
#include <exception>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace model
{
    // Independent indexed outputs; exceptions rethrown in index order. The
    // reference gate is for regression tests and benchmarks, not physics inputs.
    template<class Function> void independentParallelFor(std::size_t count,const Function& function,std::size_t threshold=256)
    {
        bool parallel=false;
#ifdef _OPENMP
        parallel=count>=threshold && !omp_in_parallel() && !std::getenv("MODELIB_SERIAL_ASSEMBLY");
#endif
        if(!parallel)
        {
            for(std::size_t i=0;i<count;++i) { function(i); }
            return;
        }
        std::vector<std::exception_ptr> failures(count);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(count);++i)
        {
            try { function(static_cast<std::size_t>(i)); }
            catch(...) { failures[i]=std::current_exception(); }
        }
        for(const auto& failure:failures)
        {
            if(failure) { std::rethrow_exception(failure); }
        }
    }
    template<class T> void appendOrderedRows(std::vector<T>& output,const std::vector<std::vector<T>>& rows)
    {
        std::vector<std::size_t> offsets(rows.size()+1,output.size());
        for(std::size_t i=0;i<rows.size();++i) { offsets[i+1]=offsets[i]+rows[i].size(); }
        output.resize(offsets.back());
        independentParallelFor(rows.size(),[&](std::size_t i)
        {
            std::copy(rows[i].begin(),rows[i].end(),output.begin()+offsets[i]);
        });
    }
}
#endif
