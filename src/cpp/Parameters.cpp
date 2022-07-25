#include "Parameters.h"
#include "SSA.h"
#include "TauLeap.h"
#include "ImplicitTauLeap.h"
#include "TrapezoidalTauLeap.h"
#include "SK-tau-ROCK.h"
#include "tau-ROCK.h"
#include "R-tau-ROCK.h"
#include "SplitStepImplicitTau.h"
#include "ProblemsList.h"
#include <GetPot>

Parameters::Parameters()
{
}

Parameters::~Parameters()
{
}

void Parameters::init()
{
    ChemicalSystem* cs;
    init_ChemicalSystem(cs);
    N_species = cs->get_num_species();
    problem_name = cs->get_problem_name();
    final_time = cs->get_final_time();
    delete cs;
    
    filename  = problem_name + "/" + filename;
    
    if(refsol != string(""))
        refsol  = problem_name + "/" + refsol;
    
    if(solver_name == string("none")) //we dont simulate but just process data
    {
        MCiter =0;
        tau_max=0.;
    }
}

void Parameters::read_command_line(int argc, char** argv)
{
    GetPot cl(argc, argv);
   
    prob_num = cl.follow(prob_num, "-prob");
    MCiter = cl.follow(MCiter, 2, "-mc", "-mciter");
    tau_max = cl.follow(tau_max, 2, "-tau", "-dt");
    nout = cl.follow(nout, "-nout");
    filename  = cl.follow(filename.c_str(), "-ofile");
    post_process = cl.follow((int)post_process, 3, "-post_proc", "-pp", "-p");
    n_bins = cl.follow(n_bins, 2, "-nbins", "-bins");
    refsol = cl.follow(refsol.c_str(), 2, "-refsol", "-ref");
    Newton_tol = cl.follow(Newton_tol, 2, "-Ntol", "-ntol");
    s_add = cl.follow(s_add, "-s_add");
    s = cl.follow(s, "-s");
    damping = cl.follow(damping, 4, "-damping", "-damp", "-eps", "-eta");
    
    if(cl.search(2,"-sol", "-solver"))
    {
        string sol = cl.next(solver_name.c_str());
        if(sol=="tau-leap" || sol=="etl" || sol=="exp-tau-leap" || sol=="tl")
        {
            solver_name = "tau-leap";
            post_process = false;
        }
        else if(sol=="imp-tau-leap" || sol=="itl")
            solver_name = "imp-tau-leap";
        else if(sol=="trap-tau-leap" || sol=="ttl")
            solver_name  = "trap-tau-leap";
        else if(sol=="sk-tau-rock" || sol=="SK-tau-ROCK" || sol=="str")
            solver_name = "SK-tau-ROCK";
        else if(sol=="tau-rock" || sol=="tau-ROCK" || sol=="tr")
            solver_name = "tau-ROCK";
        else if(sol=="r-tau-rock" || sol=="R-tau-ROCK" || sol=="rtr")
            solver_name = "R-tau-ROCK";
        else if(sol=="SSA" || sol=="ssa")
            solver_name = "SSA";
        else if(sol=="SSITL" || sol=="ssitl" || sol=="ssi")//split step implicit tau leap
            solver_name = "SSITL";
        else if(sol=="none" || sol=="NONE") //in this case we must provide two sets of data to process
            solver_name = "none";
        else
        {
            solver_name = "error";
            cout<<"ERROR: solver not implemented"<<endl;
        }
    }
    
    // for these methods there is no postprocessing implemented
    if(solver_name== "trap-tau-leap" || solver_name =="SSA" || 
       solver_name =="tau-ROCK" || solver_name =="R-tau-ROCK")
        post_process = false;
    
}

void Parameters::init_ChemicalSystem(ChemicalSystem*& cs)
{
    if(prob_num==1)
        cs = new ReversibleIsomerization;
    else if(prob_num==2)
        cs = new NonlinearReversibleReaction;
    else if(prob_num==3)
        cs = new GeneticPositiveFeedbackLoop;
    else if(prob_num==4)
        cs = new MichaelisMenten;
    else if(prob_num==5)
        cs = new SchloglReaction;
    else if(prob_num==6)
        cs = new DecayingDimerizing;
    else if(prob_num==7)
        cs = new Ecoli;
    else
        cout<<"Problem not found"<<endl;
}

void Parameters::init_Solver_ChemicalSystem(Solver*& sol, ChemicalSystem*& cs)
{
    init_ChemicalSystem(cs);
    
    if(solver_name=="SSA")
        sol = new SSA(*this,cs);
    else if(solver_name=="tau-leap" || solver_name=="tl" || solver_name=="etl" || solver_name=="exp-tau-leap")
        sol = new TauLeap(*this,cs);
    else if(solver_name=="imp-tau-leap" || solver_name=="itl")
        sol = new ImplicitTauLeap(*this,cs);
    else if(solver_name=="trap-tau-leap" || solver_name=="ttl")
        sol = new TrapezoidalTauLeap(*this,cs);
    else if(solver_name=="sk-tau-rock" || solver_name=="SK-tau-ROCK" || solver_name=="str")
        sol = new SKtauROCK(*this,cs);
    else if(solver_name=="tau-rock" || solver_name=="tau-ROCK" || solver_name=="tr")
        sol = new tauROCK(*this,cs);
    else if(solver_name=="r-tau-rock" || solver_name=="R-tau-ROCK" || solver_name=="rtr")
        sol = new RtauROCK(*this,cs);
    else if(solver_name=="SSITL" || solver_name=="ssitl" || solver_name=="ssi")
        sol = new SplitStepImplicitTau(*this,cs);
    else
        cout<<"Solver not found"<<endl;
}

string Parameters::get_filename()
{
    return filename;
}

unsigned int Parameters::get_num_species()
{
    return N_species;
}

unsigned int Parameters::get_MCiter()
{
    return MCiter;
}

unsigned int Parameters::get_nout()
{
    return nout;
}

bool Parameters::get_post_proc()
{
    return post_process;
}

unsigned int Parameters::get_n_bins()
{
    return n_bins;
}

string Parameters::get_refsol()
{
    return refsol;
}

Real Parameters::get_tau_max()
{
    return tau_max;
}

string Parameters::get_solver()
{
    return solver_name;
}

Real Parameters::get_final_time()
{
    return final_time;
}

Real Parameters::get_Newton_tol()
{
    return Newton_tol;
}