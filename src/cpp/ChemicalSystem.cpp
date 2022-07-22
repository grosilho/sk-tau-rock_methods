#include "headers.h"
#include "ChemicalSystem.h"

ChemicalSystem::ChemicalSystem(const int N_, const int M_)
:N(N_), M(M_), X0(N_), nu(N_,M_), c(M_)
{
    
}

ChemicalSystem::~ChemicalSystem()
{
}

Real ChemicalSystem::get_final_time()
{
    return T;
}

VectorXd ChemicalSystem::get_X0()
{
    return X0;
}

const MatrixXd& ChemicalSystem::get_nu()
{
    return nu;
}

int ChemicalSystem::get_num_species()
{
    return N;
}

int ChemicalSystem::get_num_reactions()
{
    return M;
}

string ChemicalSystem::get_problem_name()
{
    return problem_name;
}