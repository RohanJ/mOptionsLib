/*******************************************************************************
 * @author		: Rohan Jyoti
 * @filename	: mNormalDistLib.hpp
 * @purpose		: Normal Distribution Library Header File
 ******************************************************************************/

#ifndef mNormalDistLib_hpp
#define mNormalDistLib_hpp

double mCalcNormalDistribution(double t_val);
double mCalcNormalCDF(double t_val); //univariate
double mCalcNormalCDF(double t_val1, double t_val2, double rho); //bivariate

#endif