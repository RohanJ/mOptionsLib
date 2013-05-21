/*******************************************************************************
 * @author		: Rohan Jyoti
 * @filename	: mOptionsLib.hpp
 * @purpose		: Options Library for handling Greek calculations
 ******************************************************************************/
#include "mOptionsLib.hpp"
#include <boost/python.hpp>

//qualifier / namespace definitions
using std::cout;
using std::endl;


//===============================================
//===============Public Functions================
//===============================================

BOOST_PYTHON_MODULE(mOptionsLib)
{
	using namespace boost::python;
	class_<mOption>("mOption")
		.def("mSetOptionParam", &mOption::mSetOptionParam)
		.def("mSetImpliedVolatilityCalculationMethod", &mOption::mSetImpliedVolatilityCalculationMethod)
		.def("mSetCurrentOptionPrice", &mOption::mSetCurrentOptionPrice)
		.def("mOptionExecEngine", &mOption::mOptionExecEngine)
		.def("printMOption", &mOption::printMOption)
		.def("mGetOptionMetric", &mOption::mGetOptionMetric)
	;
}

//To use boost python module:
//Run python interpretter or script
//>>> import mOptionsLib
//>>> mOpt = mOptionsLib.mOption()
//>>> mOpt.mSetOptionParam(0, 2) //Set Option Type to Call
//>>> mOpt.mSetOptionParam(1, 100) //Set Option Underlying Price to 100
//>>> mOpt.mSetOptionParam(2, 100) //Set Option Strike to 100
//>>> mOpt.mSetOptionParam(3, 0.10) //Set Option Risk Free Rate to 0.10
//>>> mOpt.mSetOptionParam(4, 0.25) //Set Option Volatility to 0.25
//>>> mOpt.mSetOptionParam(5, 365) //Set Option Days Till Maturity to 365

//>>> mOpt.mOptionExecEngine(2) //European Option
//>>> mOpt.printMOption()
//>>> delete mOpt

//===============EUROPEAN OPTION===============
//Option Type: CALL
//Option Underlying Price: 100
//Option Strike: 100
//Option Interest Rate: 0.1
//Option Volatility: 0.25	(User Specified)
//Option Maturity: 365
//Option Price: 14.9758
//Delta: 0.700208
//Gamma: 0.0122735
//Vega: 30.6837
//Theta: -9.33997
//Rho: 55.045
//===============|||


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mSetOptionParam
 * @class	: mOption
 * @param	: t_param (the parameter), tValue (the parameter value)
 * @return	: int (PSET or PNSET indicating whether the parameter was set or not)
 * @purpose	: Sets option parameters (must be called before execEngine)
 ******************************************************************************/
int mOption::mSetOptionParam(int t_param, double tValue)
{
	if(t_param == OPTION_TYPE)
	{
		if((int)tValue != OPTION_TYPE_PUT && (int)tValue != OPTION_TYPE_CALL)
		{
			fprintf(stderr, "Invalid Option Type\n");
			return PNSET;
		}
		
		if(option_type == PNSET)
			param_counter++;
		option_type = (int)tValue;
		assert_array[OPTION_TYPE] = PSET;
		return PSET;
	}
	else if(t_param == OPTION_UNDERLYING_PRICE)
	{
		if(tValue < 0)
		{
			fprintf(stderr, "Invalid Option Underlying Price\n");
			return PNSET;
		}
		
		if(option_underlying_price == PNSET)
			param_counter++;
		option_underlying_price = tValue;
		assert_array[OPTION_UNDERLYING_PRICE] = PSET;
		return PSET;
	}
	else if(t_param == OPTION_STRIKE)
	{
		if(tValue < 0)
		{
			fprintf(stderr, "Invalid Option Strike Price\n");
			return PNSET;
		}
		
		if(option_strike == PNSET)
			param_counter++;
		option_strike = tValue;
		assert_array[OPTION_STRIKE] = PSET;
		return PSET;
	}
	else if(t_param == OPTION_RISK_FREE_RATE)
	{
		if(tValue < 0)
		{
			fprintf(stderr, "Invalid Option Risk Free Rate\n");
			return PNSET;
		}
		
		if(option_risk_free_rate == PNSET)
			param_counter++;
		option_risk_free_rate = tValue;
		assert_array[OPTION_RISK_FREE_RATE] = PSET;
		return PSET;
	}
	else if(t_param == OPTION_VOLATILITY)
	{
		if(tValue < 0)
		{
			fprintf(stderr, "Invalid Option Volatility\n");
			return PNSET;
		}
		
		if(option_volatility == PNSET)
			param_counter++;
		option_volatility = tValue;
		assert_array[OPTION_VOLATILITY] = PSET;
		return PSET;
	}
	else if(t_param == OPTION_DAYS_TILL_MATURITY)
	{
		if(tValue < 0)
		{
			fprintf(stderr, "Invalid Option Days Till Maturity\n");
			return PNSET;
		}
		
		if(option_days_till_maturity == PNSET)
			param_counter++;
		option_days_till_maturity = tValue;
		option_years_till_maturity = tValue/(double)365.0;
		assert_array[OPTION_DAYS_TILL_MATURITY] = PSET;
		return PSET;
	}
	else
	{
		printf("Invalid paramter\n");
		return PNSET;
	}
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mSetImpliedVolatilityCalculationMethod
 * @class	: mOption
 * @param	: calculation method (Newton Raphson or Bisection)
 *				Newton Raphson is defaulted. No need to call if not Bisection
 * @return	: int (PSET or PNSET indicating whether the parameter was set or not)
 * @purpose	: Sets the implied volatility calculation method
 ******************************************************************************/
int mOption::mSetImpliedVolatilityCalculationMethod(int t_method)
{
	if(t_method == NEWTON_RAPHSON)
		implied_volatility_calc_method = NEWTON_RAPHSON;
	else if(t_method == BISECTION)
		implied_volatility_calc_method = BISECTION;
	else
	{
		fprintf(stderr, "Unknown Calculation Method. Defaulting to NEWTON RAPHSON\n");
		implied_volatility_calc_method = NEWTON_RAPHSON;
	}
	return PSET;
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mSetCurrentOptionPrice
 * @class	: mOption
 * @param	: the option price value
 * @return	: int (PSET or PNSET indicating whether the parameter was set or not)
 * @purpose	: Sets the option price used to calculate implied volatility
 ******************************************************************************/
int mOption::mSetCurrentOptionPrice(double t_val)
{
	if(option_volatility != PNSET)
	{
		fprintf(stderr, "Current Option Price discarded. Using user specified volatility\n");
		return PNSET;
	}
	if(t_val < 0)
	{
		fprintf(stderr, "Invalid Current Option Price\n");
		return PNSET;
	}
	else
	{
		current_option_price = t_val;
		return PSET;
	}
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mOptionExecEngine
 * @class	: mOption
 * @param	: exercise_type (American, European)
 * @return	: void
 * @purpose	: Calculates option metrics. Handles case for implied volatility.
 ******************************************************************************/
void mOption::mOptionExecEngine(int exercise_type)
{
	int aP_val = assertParameters();
	if(aP_val == PNSET)
	{
		fprintf(stderr, "Please check which parameters above are not set. Program terminating\n");
		exit(1);
	}
	
	if(exercise_type == AMERICAN)
	{
		option_exercise_type = AMERICAN;
		if(aP_val == CALC_IMPLIED_VOL)
			mCalcImpliedVolatility(); //note that option_exercise_type has to be set prior
		
		mAmericanOptionEngine(CALC_ALL_METRICS);
	}
	else if(exercise_type == EUROPEAN)
	{
		option_exercise_type = EUROPEAN;
		if(aP_val == CALC_IMPLIED_VOL)
			mCalcImpliedVolatility();
		
		mEuropeanOptionEngine(CALC_ALL_METRICS);
	}
	else
		fprintf(stderr, "Invalid Exercise Type\n");
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: printMOption
 * @class	: mOption
 * @param	: None
 * @return	: void
 * @purpose	: Prints appropriate mOption object parameters
 ******************************************************************************/
void mOption::printMOption()
{
	char oet[32];
	if(option_exercise_type == AMERICAN)
		strcpy(oet, "AMERICAN OPTION");
	else if(option_exercise_type == EUROPEAN)
		strcpy(oet, "EUROPEAN OPTION");
	else
		strcpy(oet, "UNRECOGNIZED EXERCISE");
	
	
	printf("===============%s===============\n", oet);
	if(option_type == OPTION_TYPE_CALL)
		cout << "Option Type: CALL" << endl;
	else
		cout << "Option Type: PUT" << endl;
	cout << "Option Underlying Price: " << option_underlying_price << endl;
	cout << "Option Strike: " << option_strike << endl;
	cout << "Option Interest Rate: " << option_risk_free_rate << endl;
	cout << "Option Volatility: " << option_volatility;
	if(current_option_price == PNSET)
		cout << "	(User Specified)" << endl;
	else
	{
		if(implied_volatility_calc_method == NEWTON_RAPHSON)
			cout << "	(Implied Volatiliy using Newton Raphson)" << endl;
		else
			cout << "	(Implied Volatiliy using Bisection)" << endl;
	}
	
	cout << "Option Maturity: " << option_days_till_maturity << endl;
	
	std::string msgs[] = {"Option Price: ", "Delta: ", "Gamma: ", "Vega: ", "Theta: ", "Rho: "};
	
	int i; //represents the definition indexes for the Greeks
	for(i=0; i<6; i++)
		printMetric(i, msgs[i]);
	
	printf("===============|||\n\n");
}

/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mGetOptionMetric
 * @class	: mOption
 * @param	: The Option Metric Value to get
 * @return	: The Option Metric Value
 * @purpose	: Returns current mOption object's specified metric value
 ******************************************************************************/
double mOption::mGetOptionMetric(int t_metric)
{
	int assertation = assertParameters();
	if(assertation == PNSET)
		return PNSET;
	
	if(t_metric == METRIC_PRICING)
		return option_price;
	else if(t_metric == METRIC_DELTA)
		return option_delta;
	else if(t_metric == METRIC_GAMMA)
		return option_gamma;
	else if(t_metric == METRIC_VEGA)
		return option_vega;
	else if(t_metric == METRIC_THETA)
		return option_theta;
	else if(t_metric == METRIC_RHO)
		return option_rho;
	else if(t_metric == OPTION_VOLATILITY) //useful for when implied vol is calc'd
		return option_volatility;
	else
	{
		fprintf(stderr, "Metric Not Found/Accessible\n");
		return PNSET;
	}
}





//===============================================
//===============Private Functions===============
//===============================================

/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: assertParameters
 * @class	: mOption
 * @param	: None
 * @return	: Int indicating whether or not all parameters have been set
 * @purpose	: Make sure all required option parameters have been set and 
 *				determine if system is ready to calculate implied volatility
 ******************************************************************************/
int mOption::assertParameters()
{
	if(param_counter == NUM_PARAMS)
		return PSET;
	else
	{
		int volatility_flag = PNSET; //reps if sys is ready to calc implied vol instead
		int i;
		for(i=0; i<NUM_PARAMS; i++)
		{
			if(assert_array[i] == PNSET)
			{
				if(i==0)
					printf("Parameter 'Option Type' not set\n");
				else if(i==1)
					printf("Parameter 'Option Underlying Price' not set\n");
				else if(i==2)
					printf("Parameter 'Option Strike Price' not set\n");
				else if(i==3)
					printf("Parameter 'Option Risk Free Rate' not set\n");
				else if(i==4)
				{
					//if there is only one missing parameter and it happens
					//to be volatility, then we can return a different value
					//to indicate that implied volatility needs to be calculated
					if(param_counter == NUM_PARAMS - 1)
						volatility_flag = PSET;
					else
						printf("Parameter 'Option Volatility' not set\n");
				}
				else if(i==5)
					printf("Parameter 'Option Days Till Maturity' not set\n");
				else
					printf("Unknown paramter not set\n");
			}
		}
		
		if(volatility_flag == PSET)
		{
			if(current_option_price != PNSET)
				return CALC_IMPLIED_VOL;
			else
			{
				fprintf(stderr, "Current Option Price needs to be set if Volatility not spcecified. Only then can Implied Volatility be accordingly calculated\n");
				return PNSET;
			}
		}
		else
		{
			printf("About to return PNSET from assertParameters LOC=1\n");
			return PNSET;
		}
	}
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mEuropeanOptionEngine
 * @class	: mOption
 * @param	: Calculation Flag (CALC_OPTION_PRICE_ONLY or CALC_ALL_METRICS)
 * @return	: void
 * @purpose	: Calculates accordingly the option price and greek metrics for
 *				an European Option.
 ******************************************************************************/
void mOption::mEuropeanOptionEngine(int calculation_flag)
{
	//========== Notes ==========
	//This function needs the following:
	//1. Underlying Price
	//2. Strike Price
	//3. Interest Rate
	//4. Volatility (you can use implied volatility engine)
	//5. Time to maturity (years to maturyity is used)
	
	//printf("All parameters were set\nOptions_years_till_maturity: %f\n", option_years_till_maturity);
	double yTMaturity_sqrt = sqrt(option_years_till_maturity);
	
	//temp1 and temp2 are components of the black-scholes equation
	double temp1 = (log(option_underlying_price / option_strike) + option_risk_free_rate * option_years_till_maturity) / (option_volatility * yTMaturity_sqrt) + 0.5*option_volatility*yTMaturity_sqrt;
	double temp2 = temp1 - (option_volatility * yTMaturity_sqrt);
	
	//Next we will compute the following option metrics:
	//1. Option Price
	//2. Delta
	//3. Gamma
	//4. Theta
	//5. Vega
	//6. Rho
	
	if(option_type == OPTION_TYPE_CALL)
		option_price = option_underlying_price * mCalcNormalCDF(temp1) - option_strike*exp(-option_risk_free_rate*option_years_till_maturity)*mCalcNormalCDF(temp2);
	else if(option_type == OPTION_TYPE_PUT)
		option_price = option_strike*exp(-option_risk_free_rate*option_years_till_maturity)*mCalcNormalCDF(-temp2) - option_underlying_price*mCalcNormalCDF(-temp1);
	else{} //never reached
	
	if(calculation_flag == CALC_OPTION_PRICE_ONLY)
		return;
	
	option_delta = mCalcNormalCDF(temp1);
	option_gamma = mCalcNormalDistribution(temp1) / (option_underlying_price * option_volatility * yTMaturity_sqrt);
	option_theta = - (option_underlying_price * option_volatility * mCalcNormalDistribution(temp1)) / (2 * yTMaturity_sqrt) - option_risk_free_rate*option_strike*exp(-option_risk_free_rate*option_years_till_maturity)*mCalcNormalCDF(temp2);
	option_vega = option_underlying_price * yTMaturity_sqrt * mCalcNormalDistribution(temp1);
	option_rho = option_strike * option_years_till_maturity*exp(-option_risk_free_rate*option_years_till_maturity)*mCalcNormalCDF(temp2);
	
	//printf("delta: %f\ngamma: %f\ntheta: %f\nvega: %f\nrho: %f\n", delta, gamma, theta, vega, rho);
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mAmericanOptionEngine
 * @class	: mOption
 * @param	: Calculation Flag (CALC_OPTION_PRICE_ONLY or CALC_ALL_METRICS)
 * @return	: void
 * @purpose	: Calculates accordingly the option price and greek metrics for
 *				an American Option.
 ******************************************************************************/
void mOption::mAmericanOptionEngine(int calculation_flag)
{	
	//========== Notes ==========
	//This function needs the following:
	//1. Underlying Price
	//2. Strike Price
	//3. Interest Rate
	//4. Volatility (you can use implied volatility engine)
	//5. Time to maturity (years to maturity is used)
	//6. Number of step for binomial search
	
	
	//First we calculate the option_price
	
	//call
	const int num_steps = 100;

	//iRS is the interest rate per step
	double iRS = exp(option_risk_free_rate * (option_years_till_maturity/num_steps));
	double iRS_inverse = 1.0/iRS;
	
	double up_movement = exp(option_volatility * sqrt(option_years_till_maturity/num_steps));
	double up_movement_squared = up_movement * up_movement;
	double down_movement = 1.0/up_movement;
	double price_up = (iRS - down_movement) / (up_movement - down_movement);
	double price_down = 1.0 - price_up;
	
	std::vector<double> prices_vector(num_steps+1);
	prices_vector[0] = option_underlying_price * pow(down_movement, num_steps);
	
	int pv_index;
	for(pv_index = 1; pv_index <= num_steps; pv_index++)
		prices_vector[pv_index] = up_movement_squared * prices_vector[pv_index-1];
	
	//Because American Options can be exercised at any point leading up to maturity:
	if(option_type == OPTION_TYPE_CALL)
	{
		std::vector<double> call_vector(num_steps+1);
		int cv_index;
		for(cv_index = 0; cv_index <= num_steps; cv_index++)
			call_vector[cv_index] = std::max(0.0, (prices_vector[cv_index] - option_strike));
	
		int curr_step;
		for(curr_step = num_steps - 1; curr_step>=0; curr_step--)
		{
			int i;
			for(i=0; i<=curr_step; i++)
			{
				call_vector[i] = (price_up * call_vector[i+1] + price_down * call_vector[i]) * iRS_inverse;
				prices_vector[i] = down_movement * prices_vector[i+1];
				call_vector[i] = std::max(call_vector[i], (prices_vector[i] - option_strike));
			}
		}
	
		option_price = call_vector[0];
	}
	else if(option_type == OPTION_TYPE_PUT)
	{
		std::vector<double> put_vector(num_steps+1);
		int putV_index;
		for(putV_index = 0; putV_index <= num_steps; putV_index++)
			put_vector[putV_index] = std::max(0.0, (option_strike - prices_vector[putV_index]));
		
		int curr_step;
		for(curr_step = num_steps - 1; curr_step>=0; curr_step--)
		{
			int i;
			for(i=0; i<=curr_step; i++)
			{
				put_vector[i] = (price_up * put_vector[i+1] + price_down * put_vector[i]) * iRS_inverse;
				prices_vector[i] = down_movement * prices_vector[i+1];
				put_vector[i] = std::max(put_vector[i], (option_strike - prices_vector[i]));
			}
		}
		
		option_price = put_vector[0];
	}
	else {}
	
	if(calculation_flag == CALC_OPTION_PRICE_ONLY)
		return;
	
	//Next we will compute the following option metrics:
	//1. Delta
	//2. Gamma
	//3. Theta
	//4. Vega
	//5. Rho
	mCalcAmericanOptionPartials(iRS, iRS_inverse, up_movement, up_movement_squared, down_movement, price_up, price_down, num_steps);
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mCalcAmericanOptionPartials
 * @class	: mOption
 * @param	: iRS (Interest rate per step), iRS_inverse (inverse of iRS), 
 *				up_movement (binomial up), up_movment_squared, down_movement,
 *				price_up, price_down, num_steps (number of step in binomial)
 * @return	: void
 * @purpose	: Calculates the partials (greek metrics); Serves as American Option
 *				Engine helper function
 ******************************************************************************/
void mOption::mCalcAmericanOptionPartials(double iRS, double iRS_inverse,
										  double up_movement, double up_movement_squared,
										  double down_movement, double price_up,
										  double price_down, int num_steps)
{
	std::vector<double> prices_vector(num_steps+1);
	std::vector<double> call_put_values(num_steps+1);
	
	
	//Currently set up for CALL
	prices_vector[0] = option_underlying_price * pow(down_movement, num_steps);
	int pv_index, cpv_index;
	
	for(pv_index = 1; pv_index <= num_steps; pv_index++)
		prices_vector[pv_index] = up_movement_squared * prices_vector[pv_index-1];
	
	if(option_type == OPTION_TYPE_CALL)
	{
		for(cpv_index = 0; cpv_index <= num_steps; cpv_index++)
			call_put_values[cpv_index] = std::max(0.0, (prices_vector[cpv_index]-option_strike));
	}
	else if(option_type == OPTION_TYPE_PUT)
	{
		for(cpv_index = 0; cpv_index <= num_steps; cpv_index++)
			call_put_values[cpv_index] = std::max(0.0, (option_strike - prices_vector[cpv_index]));
	}
	else{}
	
	int curr_step;
	for(curr_step = num_steps-1; curr_step>=2; curr_step--)
	{
		for(int i=0; i<=curr_step; i++)
		{
			prices_vector[i] = down_movement * prices_vector[i+1];
			call_put_values[i] = (price_up * call_put_values[i+1] + price_down * call_put_values[i])*iRS_inverse;
			if(option_type == OPTION_TYPE_CALL)
				call_put_values[i] = std::max(call_put_values[i], prices_vector[i] - option_strike);
			else if(option_type == OPTION_TYPE_PUT)
				call_put_values[i] = std::max(call_put_values[i], option_strike - prices_vector[i]);
			else{}
		}
	}
	
	//Next we store some values which will be used later in partial derivative calculations
	double cv2_0 = call_put_values[0], cv2_1 = call_put_values[1], cv2_2 = call_put_values[2];
	
	//compound
	int i;
	for(i=0; i<=1; i++)
	{
		prices_vector[i] = down_movement * prices_vector[i+1];
		call_put_values[i] = (price_up * call_put_values[i+1] + price_down * call_put_values[i])*iRS_inverse;
		if(option_type == OPTION_TYPE_CALL)
			call_put_values[i] = std::max(call_put_values[i], prices_vector[i] - option_strike);
		else if(option_type == OPTION_TYPE_PUT)
			call_put_values[i] = std::max(call_put_values[i], option_strike - prices_vector[i]);
		else{}
	}
	
	double cv1_0 = call_put_values[0], cv1_1 = call_put_values[1];
	prices_vector[0] = down_movement * prices_vector[1];
	
	//compound
	call_put_values[0] = (price_up * call_put_values[1] + price_down * call_put_values[0])*iRS_inverse;
	if(option_type == OPTION_TYPE_CALL)
		call_put_values[0] = std::max(call_put_values[0], option_underlying_price - option_strike);
	else if(option_type == OPTION_TYPE_PUT)
		call_put_values[0] = std::max(call_put_values[0], option_strike - prices_vector[i]);
	else{}
	
	double cv0_0 = call_put_values[0];
	
	option_delta = (cv1_1 - cv1_0) / (option_underlying_price * (up_movement - down_movement));
	
	double g_denom = 0.5 * option_underlying_price * (up_movement_squared - down_movement * down_movement);
	option_gamma = ( (cv2_2 - cv2_1) / (option_underlying_price * (up_movement_squared - 1.0)) - (cv2_1 - cv2_0) / (option_underlying_price * (1-down_movement*down_movement)) ) / g_denom;
	
	option_theta = (cv2_1 - cv0_0) / (2 * (option_years_till_maturity / num_steps));
	
	double orig_vol_place_holder = option_volatility; //becasue vol and opt price get
	double orig_opt_price_place_holder = option_price; //changed and we must reset them
	double vega_threshold = 0.02;
	option_volatility += vega_threshold;
	mAmericanOptionEngine(CALC_OPTION_PRICE_ONLY);
	double temp_price = option_price; //just calculated from mAmericanOptionEngine
	option_vega = (temp_price - cv0_0) / vega_threshold;
	//next we restore back from originals
	option_volatility = orig_vol_place_holder;
	option_price = orig_opt_price_place_holder;
	
	double orig_risk_rate_place_holder = option_risk_free_rate;
	double rho_threshold = 0.05;
	option_risk_free_rate += rho_threshold;
	mAmericanOptionEngine(CALC_OPTION_PRICE_ONLY);
	temp_price = option_price;
	option_rho = (temp_price - cv0_0) / rho_threshold;
	//next we restore back from originals
	option_risk_free_rate = orig_risk_rate_place_holder;
	option_price = orig_opt_price_place_holder;
	
	
	//Due to demand for theta and rho to be in $/day instead of $/year:
	option_theta = option_theta / 365.0;
	option_rho = option_rho / 365.0;
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mCalcImpliedVolatility
 * @class	: mOption
 * @param	: void
 * @return	: void
 * @purpose	: Wrapper function for implied volatility calculation
 ******************************************************************************/
void mOption::mCalcImpliedVolatility()
{
	//Note that this function can only be called from mOptionExecEngine
	//and only after assertParameters has confirmed the system is ready
	//to calculate implied volatility (i.e. all params except volatility are set)
	if(implied_volatility_calc_method == NEWTON_RAPHSON)
	{
		option_volatility = mCIV_NewtonRaphson();
	}
	else if(implied_volatility_calc_method == BISECTION)
	{
		option_volatility = mCIV_Bisection();
	}
	else{}
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mCIV_Bisection
 * @class	: mOption
 * @param	: void
 * @return	: The option volatility (redundant)
 * @purpose	: Calculates implied volatility using Bisection
 ******************************************************************************/
double mOption::mCIV_Bisection()
{
	//It is assumed that by this point current_option_price has been set
	//This function uses binomial search/bisection and is much faster than
	//Newton-Raphson but does not always find a solution
	
	const double max_val = 10000000000;
	const double bisection_max_num_rounds = 100;
	const double method_accuracy = 0.00001;
	
	
	if(current_option_price < 0.99*(option_underlying_price - option_strike*exp(-option_years_till_maturity*option_risk_free_rate)))
		return 0.0;
	
	//First we calculated a maximum volatility
	double vol_low = 0.00001, vol_high = 0.3;
	option_volatility = vol_high;
	mEuropeanOptionEngine(CALC_OPTION_PRICE_ONLY);
	double temp_price = option_price; //just calculated by mEuropeanOptionEngine
	
	while (temp_price < current_option_price)
	{
		vol_high*=2.0;
		option_volatility = vol_high;
		mEuropeanOptionEngine(CALC_OPTION_PRICE_ONLY);
		temp_price = option_price;
		if(vol_high > max_val)
		{
			fprintf(stderr, "FATAL ERROR in mCIV_Bisection LOC=1\n");
			exit(1);
		}
	}
	
	//Next we perform binomial search to attain implied volatility
	int b_index;
	for(b_index = 0; b_index < bisection_max_num_rounds; b_index++)
	{
		option_volatility = (vol_low + vol_high)*0.5;
		mEuropeanOptionEngine(CALC_OPTION_PRICE_ONLY);
		temp_price = option_price;
		
		double assert_accuracy = (temp_price - current_option_price);
		if(fabs(assert_accuracy) < method_accuracy)
			return option_volatility;
		else if(assert_accuracy < 0.0)
			vol_low = option_volatility;
		else
			vol_high = option_volatility;
	}
	
	//Should theoretically never reach here
	fprintf(stderr, "FATAL ERROR in mCIV_Bisection LOC=2\n");
	exit(1);
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: mCIV_NewtonRaphson
 * @class	: mOption
 * @param	: void
 * @return	: The option volatility (redundant)
 * @purpose	: Calculates implied volatility using Newton Raphson
 ******************************************************************************/
double mOption::mCIV_NewtonRaphson()
{
	//It is assumed that by this point current_option_price has been set
	//This function uses is slower than the bisection method but has a
	//much higher probability of finding a solution
	
	const double nr_max_num_rounds = 100;
	const double method_accuracy = 0.00001;

	if(current_option_price < 0.99*(option_underlying_price - option_strike*exp(-option_years_till_maturity*option_risk_free_rate)))
		return 0.0;
	
	//First we find an initial volatility value
	double yTMaturity_sqrt = sqrt(option_years_till_maturity);
	option_volatility = (current_option_price / option_underlying_price)/(0.398*yTMaturity_sqrt);
	
	int nrIndex;
	for(nrIndex = 0; nrIndex < nr_max_num_rounds; nrIndex++)
	{
		mEuropeanOptionEngine(CALC_OPTION_PRICE_ONLY);
		double temp_price = option_price;
		double assert_accuracy = (current_option_price - temp_price);
		if(fabs(assert_accuracy) < method_accuracy)
			return option_volatility;
		else
		{
			double temp1 = (log(option_underlying_price / option_strike) + option_risk_free_rate * option_years_till_maturity) / (option_volatility * yTMaturity_sqrt) + 0.5*option_volatility*yTMaturity_sqrt;
			
			option_vega = option_underlying_price * yTMaturity_sqrt * mCalcNormalDistribution(temp1);
			option_volatility+=(assert_accuracy/option_vega);
		}
	}
	
	//Should theoretically never reach here
	fprintf(stderr, "FATAL ERROR in mCIV_NewtonRaphson LOC=1\n");
	exit(1);
}


/*******************************************************************************
 * @author	: Rohan Jyoti
 * @name	: printMetric
 * @class	: mOption
 * @param	: The Metric to print, message to print along with metric
 * @return	: void
 * @purpose	: Prints specified metric and specified msg associated with metric
 ******************************************************************************/
void mOption::printMetric(int t_metric, std::string msg)
{
	double *value = NULL;
	if(t_metric == METRIC_PRICING)
		value = &option_price;
	else if(t_metric == METRIC_DELTA)
		value = &option_delta;
	else if(t_metric == METRIC_GAMMA)
		value = &option_gamma;
	else if(t_metric == METRIC_VEGA)
		value = &option_vega;
	else if(t_metric == METRIC_THETA)
		value = &option_theta;
	else if(t_metric == METRIC_RHO)
		value = &option_rho;
	else
		fprintf(stderr, "Invalid metric provided to function 'printMetric'\n");
	
	cout << msg << *value << endl;
}

