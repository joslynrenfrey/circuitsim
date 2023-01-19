/*
 * Author: Joslyn Renfrey
 * Date: 13/01/2023
 */

#include"circuitsim.h"

const int res_terminals_count = 2;
const int res_parameters_count = 1;

void res_currentCurve(const double *parameters, const double *v, double timestep, double *i){
	double res = parameters[0];
	double current = (v[0] - v[1])/res;
	i[0] =  current;
	i[1] = -current;
}

void res_jacobian(const double *parameters, const double *v, double timestep, double *j){
	double res = parameters[0];
	double slope = 1/res;
	j[0] =  slope; j[1] = -slope;
	j[2] = -slope; j[3] =  slope;
}

void res_updateState(double *parameters, const double *v, double timestep, const double *i){}



const int src_terminals_count = 2;
const int src_parameters_count = 2;

void src_currentCurve(const double *parameters, const double *v, double timestep, double *i){
	double max_v = parameters[0];
	double max_i = parameters[1];
	double res = max_v/max_i;
	double current = (v[0] - v[1] - max_v)/res;
	
	i[0] =  current;
	i[1] = -current;
}

void src_jacobian(const double *parameters, const double *v, double timestep, double *j){
	double max_v = parameters[0];
	double max_i = parameters[1];
	double res = max_v/max_i;
	
	double slope = 1/res;
	j[0] =  slope; j[1] = -slope;
	j[2] = -slope; j[3] =  slope;
}

void src_updateState(double *parameters, const double *v, double timestep, const double *i){}



const int cap_terminals_count = 2;
const int cap_parameters_count = 2;

void cap_currentCurve(const double *parameters, const double *v, double timestep, double *i){
	double cap = parameters[0];
	double past_v = parameters[1];
	double past_i = parameters[2];
	// we use a trapezoidal approximation, which has much better convergence properties
	// i = c.dv/dt
	double mean_i = ( v[0] - v[1] - past_v)*(cap/timestep);
	// mean_i = (i + past_i)/2, 2*mean_i - past_i = i
	double current = 2*mean_i - past_i;
	i[0] =  current;
	i[1] = -current;
}

void cap_jacobian(const double *parameters, const double *v, double timestep, double *j){
	double cap = parameters[0];
	double slope = 2*cap/timestep;
	j[0] =  slope; j[1] = -slope;
	j[2] = -slope; j[3] =  slope;
}

void cap_updateState(double *parameters, const double *v, double timestep, const double *i){
	parameters[1] = v[0] - v[1];
	parameters[2] = (i[0] - i[1])/2;
}



const int ind_terminals_count = 2;
const int ind_parameters_count = 2;

void ind_currentCurve(const double *parameters, const double *v, double timestep, double *i){
	double ind = parameters[0];
	double past_i = parameters[1];
	double past_v = parameters[2];
	
	// we use a trapezoidal approximation, which has much better convergence properties
	double mean_v = (v[0] - v[1] + past_v)/2;
	// v = l.di/dt
	// i + di = i + v.dt/l
	double current = mean_v*(timestep/ind) + past_i;
	i[0] =  current;
	i[1] = -current;
}

void ind_jacobian(const double *parameters, const double *v, double timestep, double *j){
	double ind = parameters[0];
	double slope = 0.5*timestep/ind;
	j[0] =  slope; j[1] = -slope;
	j[2] = -slope; j[3] =  slope;
}

void ind_updateState(double *parameters, const double *v, double timestep, const double *i){
	parameters[1] = (i[0] - i[1])/2;
	parameters[2] = v[0] - v[1];
}

