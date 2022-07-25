#ifndef PROBLEMSLIST_H
#define PROBLEMSLIST_H

#include "ChemicalSystem.h"

//Source: Y. Cao, L. R. Petzold, M. Rathinam, and D. T. Gillespie. The numerical stability of leaping methods for stochastic simulation of chemically reacting systems. Journal of Chemical Physics, 121(24):12169–12178, 2004.
class ReversibleIsomerization: public ChemicalSystem
{
public:
    ReversibleIsomerization();
    ~ReversibleIsomerization();
    
    void eval_propensity_fun(const VectorXd& X, VectorXd& a);
    
};


// VERY GOOD BABY PROBLEM
//Here we see that PSK-tau-ROCK gives excellent results independently of tau
//Trapezoidal is less accurate, slower and depends on tau. Implicit tau is awful.
//SK-tau-ROCK does not get the right variance unless we take more stages
// Similar to example 3 in 
// ﻿Abdulle, A., Hu, Y., & Li, T. (2010). Chebyshev Methods with Discrete Noise: the $\tau$-ROCK Methods. Journal of Computational Mathematics, 28(2), 195–217.
class NonlinearReversibleReaction: public ChemicalSystem
{
public:
    NonlinearReversibleReaction();
    ~NonlinearReversibleReaction();
    
    void eval_propensity_fun(const VectorXd& X, VectorXd& a);
    
protected:
    unsigned int XT;
    
};

//Source: Y. Yang, M. Rathinam, and J. Shen. Integral tau methods for stiff stochastic chemical systems. Journal of Chemical Physics, 134(4):044129, 2011.
class GeneticPositiveFeedbackLoop: public ChemicalSystem
{
public:
    GeneticPositiveFeedbackLoop();
    ~GeneticPositiveFeedbackLoop();
    
    void eval_propensity_fun(const VectorXd& X, VectorXd& a);
    
};

class MichaelisMenten: public ChemicalSystem
{
public:
    MichaelisMenten();
    ~MichaelisMenten();
    
    void eval_propensity_fun(const VectorXd& X, VectorXd& a);
};

// Problem with bimodal stationary distribution. But it is not stiff.
// Rathinam, M., Petzold, L. R., Cao, Y., and Gillespie, D. T. (2005). Consistency and Stability of Tau-Leaping Schemes. Multiscale Modeling and Simulation, 4(3):867–895.
class SchloglReaction: public ChemicalSystem
{
public:
    SchloglReaction();
    ~SchloglReaction();
    
    void eval_propensity_fun(const VectorXd& X, VectorXd& a);
    
protected:
    unsigned int N1, N2;
    
};

// Decaying-dimerizing (contains a nonlinear reversible reaction)
// Gillespie, D. T. (2001). Approximate accelerated stochastic simulation of chemically reacting systems. Journal of Chemical Physics, 115(4):1716-1733.
class DecayingDimerizing: public ChemicalSystem
{
public:
    DecayingDimerizing();
    ~DecayingDimerizing();
    
    void eval_propensity_fun(const VectorXd& X, VectorXd& a);
    
protected:
    
};

// Large nonlinear problem without scale separation
class Ecoli: public ChemicalSystem
{
public:
    Ecoli();
    ~Ecoli();
    
    void eval_propensity_fun(const VectorXd& X, VectorXd& a);
    
protected:
    
};

#endif /* PROBLEMSLIST_H */

