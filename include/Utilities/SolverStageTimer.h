/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_SolverStageTimer_H_
#define model_SolverStageTimer_H_
#include <chrono>
#include <exception>
#include <iostream>
#include <string>
namespace model
{
    // Flushed starts identify long-running stages before they finish.
    class SolverStageTimer
    {
        using Clock=std::chrono::steady_clock;
        std::string label;
        Clock::time_point started;
        int exceptions=0;
        bool active=false;
    public:
        explicit SolverStageTimer(const std::string& name) { start(name); }
        SolverStageTimer(const SolverStageTimer&)=delete;
        SolverStageTimer& operator=(const SolverStageTimer&)=delete;
        void finish()
        {
            if(active)
            {
                const double seconds=std::chrono::duration<double>(Clock::now()-started).count();
                std::cout<<"\n  "<<label<<(std::uncaught_exceptions()>exceptions?" interrupted":" completed")
                         <<" ["<<seconds<<" sec]"<<std::flush;
                active=false;
            }
        }
        void start(const std::string& name)
        {
            finish(); label=name; exceptions=std::uncaught_exceptions(); active=true;
            std::cout<<"\n  "<<label<<" starting"<<std::flush;
            started=Clock::now();
        }
        ~SolverStageTimer() { finish(); }
    };
}
#endif
