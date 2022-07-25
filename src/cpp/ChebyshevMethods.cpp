/*
    Copyright (C) 2022 Giacomo Rosilho de Souza

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "ChebyshevMethods.h"

//ChebyshevMethods::ChebyshevMethods()
//{
//    
//}
//
//ChebyshevMethods::~ChebyshevMethods()
//{
//    
//}

void ChebyshevMethods::CoefficientsRKC1(vector<Real>& mu, vector<Real>& nu, 
                       vector<Real>& kappa, unsigned int s, Real eps)
{
    mu.resize(s);
    nu.resize(s);
    kappa.resize(s);
    
    Real w0 = 1.+eps/s/s;
    Real w1 = T(w0,s)/dT(w0,s);
   
    mu[0] = w1/w0;
    nu[0] = s*w1/2./w0;
    kappa[0] = s*w1/w0;
    
    Real Tjm2, Tjm1, Tj;
    Tjm2 = 1.;
    Tjm1 = w0;
    
    for(unsigned int j=2;j<=s;j++)
    {
        Tj = 2*w0*Tjm1-Tjm2;
        
        mu[j-1] = 2*w1*Tjm1/Tj;
        nu[j-1] = 2*w0*Tjm1/Tj;
        kappa[j-1] = -Tjm2/Tj;

        Tjm2 = Tjm1;
        Tjm1 = Tj;
    }
}

Real ChebyshevMethods::StabBoundaryRKC1(unsigned int s, Real eps)
{
    
    Real w0 = 1.+eps/s/s;
    Real w1 = T(w0,s)/dT(w0,s);
   
    return 2.*w0/w1;
}


Real ChebyshevMethods::T(Real x, unsigned int s)
{
    if(s==0)
        return 1.;
    else if(s==1)
        return x;
    else
    {
        Real Tjm2 = 1.;
        Real Tjm1 = x;
        Real Tj;
        for(unsigned int j=2;j<=s;j++)
        {
            Tj = 2*x*Tjm1-Tjm2;
            Tjm2=Tjm1;
            Tjm1=Tj;
        }
        return Tj;
    }
}

Real ChebyshevMethods::U(Real x, unsigned int s)
{
    if(s==0)
        return 1.;
    else if(s==1)
        return 2.*x;
    else
    {
        Real Ujm2 = 1.;
        Real Ujm1 = 2.*x;
        Real Uj;
        for(unsigned int j=2;j<=s;j++)
        {
            Uj = 2*x*Ujm1-Ujm2;
            Ujm2=Ujm1;
            Ujm1=Uj;
        }
        return Uj;
    }
}

Real ChebyshevMethods::dT(Real x, unsigned int s)
{
    if(s==0)
        return 0.;
    else if(s==1)
        return 1.;
    else
        return s*U(x,s-1);
}

Real ChebyshevMethods::dU(Real x, unsigned int s)
{
    if(s==0)
        return 0.;
    else if(s==1)
        return 2;
    else
        return ddT(x,s+1)/(s+1);
}

Real ChebyshevMethods::ddT(Real x, unsigned int s)
{
    if(s==0)
        return 0;
    else if(s==1)
        return 0;
    else
    {
        if(x==1)
            return s*s*(s*s-1)/3.;
        else if(x==-1)
            return pow(-1,s)*s*s*(s*s-1)/3.;
        else
            return s/(x*x-1)*((s+1)*T(x,s)-U(x,s));
    }
}