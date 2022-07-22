#ifndef CHEBYSHEVMETHODS_H
#define CHEBYSHEVMETHODS_H

#include "headers.h"

namespace ChebyshevMethods
{    
    void CoefficientsRKC1(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa, 
                          unsigned int s, Real eps=0.05);
    void CoefficientsModRKC1(vector<Real>& mu, vector<Real>& nu, vector<Real>& mut, Real& theta,
                          vector<Real>& c, unsigned int s, Real eps=0.05);
    
    
    Real T(Real x, unsigned int s);
    Real dT(Real x, unsigned int s);
    Real ddT(Real x, unsigned int s);
    Real U(Real x, unsigned int s);
    Real dU(Real x, unsigned int s);
};

#endif /* CHEBYSHEVMETHODS_H */

