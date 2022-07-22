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

void ChebyshevMethods::CoefficientsModRKC1(vector<Real>& nu, vector<Real>& kappa, 
                       vector<Real>& mu, Real& theta, vector<Real>& c,
                       unsigned int s, Real eps)
{
    if(s%2!=0)
    {
        cout<<"ERROR: CoefficientsModRKC1 is used with an even s."<<endl;
        return;
    }
    
    unsigned int r=s/2;
    
    nu.resize(r);
    kappa.resize(r);
    mu.resize(r);
    c.resize(r);
    
    Real w0 = 1.+eps/s/s;
    Real w1 = T(w0,s)/dT(w0,s);
    theta = (1./(2.*w1))*(T(w0,r)/dT(w0,r));
    
    Real bjm2,bjm1,bj;
    
    bjm1 = 1./w0; // = 1/T(w0,1)
    mu[0] = w1/w0;
    nu[0] = s*w1/2.;
    kappa[0] = s*w1/w0;
    c[0] = mu[0];
    
    for(unsigned int j=2;j<=r;j++)
    {
        bj = 1./T(w0,j);
        mu[j-1] = 2*w1*bj/bjm1;
        nu[j-1] = 2*w0*bj/bjm1;
        
        if(j>2)
        {
            kappa[j-1] = -bj/bjm2;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2]+kappa[j-1]*c[j-3];
        }
        else
        {
            kappa[j-1]=-bj;
            c[j-1] = mu[j-1]+nu[j-1]*c[j-2];
        }
        
        bjm2 = bjm1;
        bjm1 = bj;
    }
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