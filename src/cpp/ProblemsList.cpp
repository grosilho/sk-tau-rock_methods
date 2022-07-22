#include "ProblemsList.h"

ReversibleIsomerization::ReversibleIsomerization()
:ChemicalSystem(2,2)
{
    problem_name = "ReversibleIsomerization";
    
    T = 10.;
    
    X0(0)=1e3;
    X0(1)=0.;
    
    nu<< -1, 1,
          1,-1;
    
    c(0)=2.;
    c(1)=1.;
}

ReversibleIsomerization::~ReversibleIsomerization()
{
}

void ReversibleIsomerization::eval_propensity_fun(const VectorXd& X, VectorXd& a)
{
    a(0)= c(0)*X(0);
    a(1)= c(1)*X(1);
}


// -----------------------------------------------------------------------------

NonlinearReversibleReaction::NonlinearReversibleReaction()
:ChemicalSystem(1,2)
{
    problem_name = "NonlinearReversibleReaction";
    
    T = 0.2;
    
    XT = 8380.;
    X0(0) = 400.;
    
    nu<< -2, 2;
    
    c(0)=50.;
    c(1)=1e3;
}

NonlinearReversibleReaction::~NonlinearReversibleReaction()
{
}

void NonlinearReversibleReaction::eval_propensity_fun(const VectorXd& X, VectorXd& a)
{
    a(0)= c(0)*X(0)*(X(0)-1.)/2.;
    a(1)= c(1)*(XT-X(0))/2.;
}

// -----------------------------------------------------------------------------

GeneticPositiveFeedbackLoop::GeneticPositiveFeedbackLoop()
:ChemicalSystem(5,9)
{
    problem_name = "GeneticPositiveFeedbackLoop";
    
    T = 100.;
    
    X0(0)=10.;
    X0(1)=0.;
    X0(2)=20.;
    X0(3)=0.;
    X0(4)=0.;
    
    nu<<-2, 2, 0, 0, 0, 0, 1,-1, 0,
         1,-1,-1, 1, 0, 0, 0, 0, 0,
         0, 0,-1, 1, 0, 0, 0, 0, 0,
         0, 0, 1,-1, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 0, 0,-1;
            
    c(0)=50.;
    c(1)=1e3;
    c(2)=50;
    c(3)=1e3;
    c(4)=1.;
    c(5)=10.;
    c(6)=3.;
    c(7)=1.;
    c(8)=6.;
}

GeneticPositiveFeedbackLoop::~GeneticPositiveFeedbackLoop()
{
}

void GeneticPositiveFeedbackLoop::eval_propensity_fun(const VectorXd& X, VectorXd& a)
{
    a(0)= c(0)*X(0)*(X(0)-1)/2;
    a(1)= c(1)*X(1);
    a(2)= c(2)*X(1)*X(2);
    a(3)= c(3)*X(3);
    a(4)= c(4)*X(2);
    a(5)= c(5)*X(3);
    a(6)= c(6)*X(4);
    a(7)= c(7)*X(0);
    a(8)= c(8)*X(4);
}

// -----------------------------------------------------------------------------

MichaelisMenten::MichaelisMenten()
:ChemicalSystem(4,3)
{
    problem_name = "MichaelisMenten";
    
    T = 50.;

    X0(0) = 3000; 
    X0(1) = 120;   
    X0(2) = 0;
    X0(3) = 0;

    nu<<-1, 1, 0,
        -1, 1, 1, 
         1,-1,-1,
         0, 0, 1;
         
    c(0)=1.66e-3;
    c(1)=1e-4;
    c(2)=1e3;
}

MichaelisMenten::~MichaelisMenten()
{
}

void MichaelisMenten::eval_propensity_fun(const VectorXd& X, VectorXd& a)
{
    a(0)= c(0)*X(0)*X(1);
    a(1)= c(1)*X(2);
    a(2)= c(2)*X(2);
}

// -----------------------------------------------------------------------------

SchloglReaction::SchloglReaction()
:ChemicalSystem(1,4)
{
    problem_name = "SchloglReaction";
    
    T = 50.;

    X0(0) = 250;
    
    nu<< 1, -1, 1, -1;
    
    c(0)=3e-7;
    c(1)=1e-4;
    c(2)=1e-3;
    c(3)=3.5;
    
    N1 = 1e5;
    N2 = 2e5;
}

SchloglReaction::~SchloglReaction()
{
}

void SchloglReaction::eval_propensity_fun(const VectorXd& X, VectorXd& a)
{
    a(0)= c(0)*N1*X(0)*(X(0)-1.)/2.;
    a(1)= c(1)*X(0)*(X(0)-1.)*(X(0)-2.)/6.;
    a(2)= c(2)*N2;
    a(3)= c(3)*X(0);
}

// -----------------------------------------------------------------------------

DecayingDimerizing::DecayingDimerizing()
:ChemicalSystem(3,4)
{
    problem_name = "DecayingDimerizing";
    
    T = 0.2;

    X0(0) = 400;
    X0(1) = 798;
    X0(2) = 0;
    
    nu<<-1,-2, 2, 0,
         0, 1,-1,-1,
         0, 0, 0, 1;
    
    c(0)=1e0;
    c(1)=1e1;
    c(2)=1e3;
    c(3)=1e-1;
   
}

DecayingDimerizing::~DecayingDimerizing()
{
}

void DecayingDimerizing::eval_propensity_fun(const VectorXd& X, VectorXd& a)
{
    a(0)= c(0)*X(0);
    a(1)= c(1)*X(0)*(X(0)-1.)/2.;
    a(2)= c(2)*X(1);
    a(3)= c(3)*X(1);
}


// -----------------------------------------------------------------------------

Ecoli::Ecoli()
:ChemicalSystem(28,61)
{
    problem_name = "Ecoli";
    
    T = 80.;

    X0 << 0,0,0,0,1,4645670,1324,80,16,3413,29,584,1,22,0,171440,9150,2280,6,596,0,13,3,3,7,0,260,0;
    
    MatrixXd tmp_nu = MatrixXd::Zero(29,62);
    
    tmp_nu(1,1) =      -1;
    tmp_nu(2,1) =      -1;
    tmp_nu(3,1) =       1;
    tmp_nu(1,2) =       1;
    tmp_nu(2,2) =       1;
    tmp_nu(3,2) =      -1;
    tmp_nu(1,3) =      -1;
    tmp_nu(4,3) =      -1;
    tmp_nu(5,3) =       1;
    tmp_nu(1,4) =       1;
    tmp_nu(4,4) =       1;
    tmp_nu(5,4) =      -1;
    tmp_nu(1,5) =      -1;
    tmp_nu(6,5) =      -1;
    tmp_nu(7,5) =       1;
    tmp_nu(1,6) =       1;
    tmp_nu(6,6) =       1;
    tmp_nu(7,6) =      -1;
    tmp_nu(4,7) =      -1;
    tmp_nu(14,7) =      -1;
    tmp_nu(15,7)  =      1;
    tmp_nu(4,8)    =    1;
    tmp_nu(14,8)   =     1;
    tmp_nu(15,8)   =    -1;
    tmp_nu(14,9)   =    -1;
    tmp_nu(16,9)   =    -1;
    tmp_nu(17,9)   =     1;
    tmp_nu(14,10)  =     1;
    tmp_nu(16,10)  =     1;
    tmp_nu(17,10)  =    -1;
    tmp_nu(3,11)   =   -1;
    tmp_nu(6,11)   =   -1;
    tmp_nu(8,11)   =    1;
    tmp_nu(3,12)   =    1;
    tmp_nu(6,12)   =    1;
    tmp_nu(8,12)   =   -1;
    tmp_nu(5,13)   =   -1;
    tmp_nu(6,13)   =   -1;
    tmp_nu(9,13)   =    1;
    tmp_nu(5,14)   =    1;
    tmp_nu(6,14)   =    1;
    tmp_nu(9,14)   =   -1;
    tmp_nu(3,15)   =   -1;
    tmp_nu(10,15)  =    -1;
    tmp_nu(12,15)  =     1;
    tmp_nu(3,16)   =    1;
    tmp_nu(10,16)  =     1;
    tmp_nu(12,16)  =    -1;
    tmp_nu(5,17)   =   -1;
    tmp_nu(11,17)  =    -1;
    tmp_nu(13,17)  =     1;
    tmp_nu(5,18)   =    1;
    tmp_nu(11,18)  =     1;
    tmp_nu(13,18)  =    -1;
    tmp_nu(15,19)  =    -1;
    tmp_nu(18,19)  =    -1;
    tmp_nu(19,19)  =     1;
    tmp_nu(15,20)   =    1;
    tmp_nu(18,20)   =    1;
    tmp_nu(19,20)  =    -1;
    tmp_nu(22,21)  =     1;
    tmp_nu(22,22)  =    -1;
    tmp_nu(14,23)  =     1;
    tmp_nu(14,24)  =    -1;
    tmp_nu(16,25)  =     1;
    tmp_nu(17,25)  =    -1;
    tmp_nu(4,26)   =    1;
    tmp_nu(15,26)  =    -1;
    tmp_nu(4,27)   =    1;
    tmp_nu(18,27)  =     1;
    tmp_nu(19,27)  =    -1;
    tmp_nu(4,28)   =    1;
    tmp_nu(27,28)  =     1;
    tmp_nu(28,28)  =    -1;
    tmp_nu(23,29)  =     1;
    tmp_nu(23,30)   =   -1;
    tmp_nu(18,31)  =     1;
    tmp_nu(18,32)   =   -1;
    tmp_nu(15,33)  =     1;
    tmp_nu(19,33)   =   -1;
    tmp_nu(25,34)  =     1;
    tmp_nu(25,35)  =    -1;
    tmp_nu(4,36)   =    1;
    tmp_nu(4,37)   =   -1;
    tmp_nu(14,38)  =     1;
    tmp_nu(18,38)  =     1;
    tmp_nu(19,38)  =    -1;
    tmp_nu(20,39)  =     1;
    tmp_nu(21,39)  =    -1;
    tmp_nu(14,40)  =     1;
    tmp_nu(27,40)   =    1;
    tmp_nu(28,40)   =   -1;
    tmp_nu(24,41)   =    1;
    tmp_nu(24,42)   =   -1;
    tmp_nu(20,43)   =    1;
    tmp_nu(20,44)  =    -1;
    tmp_nu(4,45)   =    1;
    tmp_nu(21,45)  =    -1;
    tmp_nu(4,46)   =   -1;
    tmp_nu(20,46)  =    -1;
    tmp_nu(21,46)  =     1;
    tmp_nu(4,47)   =    1;
    tmp_nu(20,47)  =     1;
    tmp_nu(21,47)  =    -1;
    tmp_nu(26,48)  =     1;
    tmp_nu(26,49)  =    -1;
    tmp_nu(27,50)   =    1;
    tmp_nu(27,51)   =   -1;
    tmp_nu(15,52)   =    1;
    tmp_nu(28,52)  =    -1;
    tmp_nu(15,53)  =    -1;
    tmp_nu(27,53)  =    -1;
    tmp_nu(28,53)  =     1;
    tmp_nu(15,54)   =    1;
    tmp_nu(27,54)  =     1;
    tmp_nu(28,54)   =   -1;
    tmp_nu(1,55)    =   1;
    tmp_nu(5,55)    =  -1;
    tmp_nu(1,56)   =    1;
    tmp_nu(11,56)   =    1;
    tmp_nu(13,56)   =   -1;
    tmp_nu(7,57)   =    1;
    tmp_nu(9,57)   =   -1;
    tmp_nu(14,58)   =    1;
    tmp_nu(15,58)   =   -1;
    tmp_nu(14,59)   =    1;
    tmp_nu(18,59)  =     1;
    tmp_nu(19,59)  =    -1;
    tmp_nu(14,60)   =    1;
    tmp_nu(20,60)   =   -1;
    tmp_nu(27,60)   =    1;
    tmp_nu(20,61)   =    1;
    tmp_nu(21,61)   =   -1;
    
    nu = tmp_nu.block(1,1,28,61);
    
    c << 2.54,1,0.254,1,0.0254,10,
         254,1e4,0.000254,0.01,0.000254,1,
         0.000254,1,2.54,1,2540,1e3,
         0.0254,1,6.62,0.5,20,0.03,
         0.03,0.03,0.03,
         0.03,1.67,0.5,20,0.03,
         0.03,0.0065,0.5,7,0.03,
         3,0.7,0.5,
         1,0.5,20,0.03,0.03,
         2.54,1e4,0.43333,0.5,20,0.03,
         0.03,2.54,1e4,0.03,
         0.03,0.03,0.03,
         0.03,0.03,0.03;
   
}

Ecoli::~Ecoli()
{
}

void Ecoli::eval_propensity_fun(const VectorXd& X, VectorXd& a)
{
            
    a(0) = c(0)*X(0)*X(1);
    a(1) = c(1)*X(2);
    a(2) = c(2)*X(0)*X(3);
    a(3) = c(3)*X(4);
    a(4) = c(4)*X(0)*X(5);
    a(5) = c(5)*X(6);
    a(6) = c(6)*X(3)*X(13);
    a(7) = c(7)*X(14);
    a(8) = c(8)*X(13)*X(15);
    a(9) = c(9)*X(16);
    a(10) = c(10)*X(2)*X(5);
    a(11) = c(11)*X(7);
    a(12) = c(12)*X(4)*X(5);
    a(13) = c(13)*X(8);
    a(14) = c(14)*X(2)*X(9);
    a(15) = c(15)*X(11);
    a(16) = c(16)*X(4)*X(10);
    a(17) = c(17)*X(12);
    a(18) = c(18)*X(14)*X(17);
    a(19) = c(19)*X(18);
    a(20) = c(20)*X(12);
    a(21) = c(21)*X(21);
    a(22) = c(22)*X(21);
    a(23) = c(23)*X(13);
    a(24) = c(24)*X(16);
    a(25) = c(25)*X(14);
    a(26) = c(26)*X(18);
    a(27) = c(27)*X(27);
    a(28) = c(28)*X(12);
    a(29) = c(29)*X(22);
    a(30) = c(30)*X(22);
    a(31) = c(31)*X(17);
    a(32) = c(32)*X(18);
    a(33) = c(33)*X(11);
    a(34) = c(34)*X(24);
    a(35) = c(35)*X(24);
    a(36) = c(36)*X(3);
    a(37) = c(37)*X(18);
    a(38) = c(38)*X(20);
    a(39) = c(39)*X(27);
    a(40) = c(40)*X(12);
    a(41) = c(41)*X(23);
    a(42) = c(42)*X(23);
    a(43) = c(43)*X(19);
    a(44) = c(44)*X(20);
    a(45) = c(45)*X(3)*X(19);
    a(46) = c(46)*X(20);
    a(47) = c(47)*X(12);
    a(48) = c(48)*X(25);
    a(49) = c(49)*X(25);
    a(50) = c(50)*X(26);
    a(51) = c(51)*X(27);
    a(52) = c(52)*X(14)*X(26);
    a(53) = c(53)*X(27);
    a(54) = c(54)*X(4);
    a(55) = c(55)*X(12);
    a(56) = c(56)*X(8);
    a(57) = c(57)*X(14);
    a(58) = c(58)*X(18);
    a(59) = c(59)*X(19);
    a(60) = c(60)*X(20);
}