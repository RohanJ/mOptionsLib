/*******************************************************************************
 * @author		: Rohan Jyoti
 * @filename	: mNormalDistLib.cpp
 * @purpose		: Normal Distribution Library Implementation File
 ******************************************************************************/
/*NOTE: based off http://www.cs.technion.ac.il/Labs/Lccn/projects/spring99/mcast/public_html/PNNI/normal_dist.C
	original author: Bernt A Oedegaard
 
 */

#include <math.h>
#include "mNormalDistLib.hpp"

#ifndef PI
#define PI 3.141592653589793238462643
#endif


double mCalcNormalDistribution(double t_val)
{
	double ndf = (1.0 / sqrt(2 * PI) ) * exp(-0.5 * t_val);
	return ndf;
}


double mCalcNormalCDF(double t_val) //univariate
{
	double b1 =  0.31938153;
    double b2 = -0.356563782;
    double b3 =  1.781477937;
    double b4 = -1.821255978;
    double b5 =  1.330274429;
    double p  =  0.2316419;
    double c2 =  0.3989423;
	
	double bound = 6.0;
    if(t_val >  bound)
		return 1.0;
    if(t_val < -bound)
		return 0.0;
	
    double a=fabs(t_val);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-t_val)*(t_val/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if(t_val < 0.0 )
		n = 1.0 - n;
    
	return n;
}


inline double f(double x, double y, double aprime, double bprime, double rho)
{
    double r = aprime*(2*x-aprime) + bprime*(2*y-bprime)
	+ 2 * rho * (x-aprime) * (y-bprime);
    return exp(r);
};


inline double sgn( double x)
{
    if (x>=0.0)  return 1.0;
    return -1.0;
};


double mCalcNormalCDF(double a, double b, double rho) //bivariate
{
	if ( (a<=0.0) && (b<=0.0) && (rho<=0.0) ) {
        double aprime = a/sqrt(2.0*(1.0-rho*rho));
        double bprime = b/sqrt(2.0*(1.0-rho*rho));
        double A[4]={0.3253030, 0.4211071, 0.1334425, 0.006374323};
        double B[4]={0.1337764, 0.6243247, 1.3425378, 2.2626645  };
        double sum = 0;
        for (int i=0;i<4;i++) {
            for (int j=0; j<4; j++) {
                sum += A[i]*A[j]* f(B[i],B[j],aprime,bprime,rho);
            };
        };
        sum = sum * ( sqrt(1.0-rho*rho)/PI);
        return sum;
    }
    else  if ( a * b * rho <= 0.0 ) {
        if ( ( a<=0.0 ) && ( b>=0.0 ) && (rho>=0.0) ) {
            return mCalcNormalCDF(a) - mCalcNormalCDF(a, -b, -rho);
        }
        else if ( (a>=0.0) && (b<=0.0) && (rho>=0.0) ) {
            return mCalcNormalCDF(b) - mCalcNormalCDF(-a, b, -rho);
        }
        else if ( (a>=0.0) && (b>=0.0) && (rho<=0.0) ) {
            return mCalcNormalCDF(a) + mCalcNormalCDF(b) - 1.0 + mCalcNormalCDF(-a, -b, rho);
        };
    }
    else  if ( a * b * rho >= 0.0 ) {
        double denum = sqrt(a*a - 2*rho*a*b + b*b);
        double rho1 = ((rho * a - b) * sgn(a))/denum;
        double rho2 = ((rho * b - a) * sgn(b))/denum;
        double delta=(1.0-sgn(a)*sgn(b))/4.0;
        return mCalcNormalCDF(a,0.0,rho1) + mCalcNormalCDF(b,0.0,rho2) - delta;
    };
    return -99.9; // should never get here
}

