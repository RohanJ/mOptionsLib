/*******************************************************************************
 * @author		: Rohan Jyoti
 * @filename	: mOptionsLib.hpp
 * @purpose		: Options Library Header File
 ******************************************************************************/
#ifndef mOptionsLib_hpp
#define mOptionsLib_hpp

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <vector>

#include "mNormalDistLib.hpp"





//========== IMPORTANT ==========
// For conistency sakes, all fx
// names use camelType; all vars
// user underscore_seperation



/* ========== Notes ==========
 An option is a derivative security; its value depends on the value (aka price) 
 of some underlying security. Let V denote this value; V_t will denote time this
 value is observed at. 
 A call option allows the option holder the right to buy some underlying asset
 at price P (the exercise price) on or before the maturity date M. Similarly, a
 put option allows holder the right to sell.
 European options can only be exercised (i.e. used) at the maturity data, whereas
 American options can be exercised on any date leading up to and including the
 maturity date.
	Call_M = max(0, V_t - P)
	Put_M = max(0, P - V_t)
 
 Black-Scholes formula provides analytical solution for European Put/Call options.
 Continuosly compounded, risk-free interest rate, volatility of underlying asset,
 time till maturity (T-t where T,t are the times measured in years)
 
 American Options are more difficult than European options because of the uncertainty
 of when the option is exercised. In general, there is no solution for American Options.
 There are special circumstances, however, when finding a solution is possible.
	1. An American Call option on non-dividend paying stock, American option price
		is the same as the European Call.
 
*/

class mOption
{
	//Private Variables
	int param_counter;
	int option_exercise_type;
	
	int			option_type;
	double		option_underlying_price;
	double		option_strike;
	double		option_risk_free_rate;
	double		option_volatility;
	double		option_days_till_maturity;

	#define NUM_PARAMS 6
	int assert_array[NUM_PARAMS];
	
	double		option_price;
	double		option_delta;
	double		option_gamma;
	double		option_vega;
	double		option_theta;
	double		option_rho;
	
	//Default calculation method for implied volatility is newton-raphson
	//It can be changed to bisections at run time using the appropriate function
	int implied_volatility_calc_method;
		#define NEWTON_RAPHSON 1
		#define BISECTION 2
	double current_option_price; //necessary for implied volatility calculation. Must be supplied by user
	
	
public:
	#define PNSET 0 //Parameter not set
	#define PSET 1 //Parameter is set
	
	mOption()//default constructor
	{
		param_counter = 0;
		int i;
		for(i=0; i<NUM_PARAMS; i++)
			assert_array[i] = PNSET;
		implied_volatility_calc_method = NEWTON_RAPHSON;
		current_option_price = PNSET;
		
		option_type = PNSET, option_underlying_price = PNSET, option_strike = PNSET,
		option_risk_free_rate = PNSET, option_volatility = PNSET,
		option_days_till_maturity = PNSET;
	};
	~mOption(){}; //destructor
		
	
	int mSetOptionParam(int t_param, double tValue);
		#define OPTION_TYPE_PUT 1
		#define OPTION_TYPE_CALL 2
	
		#define OPTION_TYPE 0
		#define OPTION_UNDERLYING_PRICE 1
		#define OPTION_STRIKE 2
		#define OPTION_RISK_FREE_RATE 3
		#define OPTION_VOLATILITY 4
		#define OPTION_DAYS_TILL_MATURITY 5
	
	int mSetImpliedVolatilityCalculationMethod(int t_method);
	int mSetCurrentOptionPrice(double t_val);
	
	void mOptionExecEngine(int exercise_type);
		#define AMERICAN 1
		#define EUROPEAN 2
	
	void printMOption();
	double mGetOptionMetric(int t_metric);
	
private:
	double option_years_till_maturity; //necessary for black-scholes formula
	int assertParameters();//--> will check if all parameter are set and determine if sys can calc implied volatility
		#define CALC_IMPLIED_VOL 2 //recall the assertParameters return either PSET (0), PNSET (1), or this define (2)

	void printMetric(int t_metric, std::string msg);
		#define METRIC_PRICING 0
		#define METRIC_DELTA 1
		#define METRIC_GAMMA 2
		#define METRIC_VEGA 3
		#define METRIC_THETA 4
		#define METRIC_RHO 5

	
	void mEuropeanOptionEngine(int calculation_flag);
		#define CALC_OPTION_PRICE_ONLY 1
		#define CALC_ALL_METRICS 2
	void mAmericanOptionEngine(int calculation_flag);
	void mCalcAmericanOptionPartials(double iRS, double iRS_inverse, double up_movement,
									 double up_movement_squared, double down_movement,
									 double price_up, double price_down, int num_steps );

	void mCalcImpliedVolatility();
	double mCIV_Bisection();
	double mCIV_NewtonRaphson();
	
	
	
};

#endif