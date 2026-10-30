#include <DefectiveCrystal/SolvedQuantityRegistry.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace model;
void require(bool condition,const std::string& message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
template<class F> void fails(F operation,const std::string& text)
{
    try
    {
        operation();
    }
    catch(const std::runtime_error& error)
    {
        require(std::string(error.what()).find(text)!=std::string::npos,error.what());
        return;
    }
    throw std::runtime_error("Expected failure: "+text);
}
SolvedQuantityContribution scalar(const std::string& name,double& value,std::function<double()> function)
{
    SolvedQuantityContribution c;
    c.provider=name;
    c.name=name;
    c.quantities={name};
    c.update=[&value,function](const SolvedQuantityContext&)
    {
        value=function();
    };
    c.state=[&value]() -> Eigen::VectorXd
    {
        return Eigen::VectorXd::Constant(1,value);
    };
    c.restore=[&value](const Eigen::VectorXd& state)
    {
        value=state(0);
    };
    c.absoluteTolerance=[]() -> Eigen::VectorXd
    {
        return Eigen::VectorXd::Constant(1,1.e-12);
    };
    return c;
}
int main()
{
    try
    {
        for(bool reversed : {false,true})
        {
            double x=0,y=0,stress=0,output=0;
            int stressCalls=0;
            SolvedQuantityRegistry registry;
            auto a=scalar("A",x,[&]() {return 1.0+0.2*y;});
            auto b=scalar("B",y,[&]() {return 2.0+0.3*x;});
            a.dependencies={{"B","",true,{}},{"Stress","",true,{}}};
            b.dependencies={{"A","",true,{}}};
            auto s=scalar("Stress",stress,[&]() {++stressCalls; return 3.0;});
            // A single operation publishes multiple outputs without duplicate execution.
            s.quantities.push_back("OtherStress");
            a.dependencies.push_back({"OtherStress","",true,{}});
            auto o=scalar("Output",output,[&]() {return x+y;});
            o.dependencies={{"A","",true,{}},{"B","",true,{}}};
            if(reversed)
            {
                registry.add(o); registry.add(b); registry.add(s); registry.add(a);
            }
            else
            {
                registry.add(a); registry.add(s); registry.add(b); registry.add(o);
            }
            registry.options.relativeTolerance=1.e-10;
            registry.finalize();
            registry.solve(12.0);
            require(std::abs(x-1.4/0.94)<1.e-10,"Coupled solution differs from analytic solution");
            require(std::abs(y-(2.0+0.3*x))<1.e-10,"Coupled second equation failed");
            require(output==x+y && stressCalls==1,"Ordering or frozen stress cache failed");
            require(registry.groups.size()==3,"Incorrect SCC groups");
        }
        double x=0,y=0;
        SolvedQuantityRegistry failure;
        auto a=scalar("A",x,[&]() {return y+1.0;});
        auto b=scalar("B",y,[&]() {return x+1.0;});
        a.dependencies={{"B","",true,{}}};
        b.dependencies={{"A","",true,{}}};
        failure.add(a); failure.add(b);
        failure.options.maxIterations=3;
        failure.options.relaxation=1.e-12;
        failure.finalize();
        fails([&]() {failure.solve(0);},"did not converge");
        require(x==0 && y==0,"Failed coupling must restore entry state");
        failure.options.relaxation=0;
        fails([&]() {failure.solve(0);},"Invalid");
        SolvedQuantityRegistry missing;
        missing.add(a);
        fails([&]() {missing.finalize();},"Missing provider");
        fails([&]() {missing.add(a);},"duplicate");
        // Conditional edge: climb feedback is absent on glide steps.
        bool feedback=false;
        SolvedQuantityRegistry conditional;
        a.dependencies={{"B","",true,[&]() {return feedback;}}};
        conditional.add(a); conditional.add(b);
        conditional.finalize();
        conditional.solve(0);
        require(conditional.groups.size()==2,"Inactive edge created a cycle");
        feedback=true;
        conditional.options.maxIterations=2;
        fails([&]() {conditional.solve(0);},"did not converge");
        require(conditional.groups.size()==1 && conditional.groups.front().cyclic,"Active edge failed to create a cycle");
        std::cout<<"PASS: quantity scheduling, convergence, caching, rollback and dynamic dependencies\n";
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
