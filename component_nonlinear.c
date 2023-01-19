/*
 * Author: Joslyn Renfrey
 * Date: 13/01/2023
 */

#include"circuitsim.h"
#include<math.h>

// a saturating exponential that has a maximum derivative
#define SATEXP_THRESHOLD 21.0
double satExp(double x){
	if(x > SATEXP_THRESHOLD){
		double y_thresh = exp(SATEXP_THRESHOLD);
		return y_thresh*(x - SATEXP_THRESHOLD) + y_thresh;
	}
	return exp(x);
}
double derivSatExp(double x){
	return exp(x > SATEXP_THRESHOLD ? SATEXP_THRESHOLD : x);
}



const int dio_terminals_count = 2;
const int dio_parameters_count = 3;

void dio_currentCurve(const double *parameters, const double *v, double timestep, double *i){
	double v_on  = parameters[0];
	double i_on  = parameters[1];
	double i_leak = parameters[2];
	
	double v_th = v_on/log(1 + i_on/i_leak);
	double current = i_leak*(satExp((v[0] - v[1])/v_th) - 1);
	
	i[0] =  current;
	i[1] = -current;
}

void dio_jacobian(const double *parameters, const double *v, double timestep, double *j){
	double v_on  = parameters[0];
	double i_on  = parameters[1];
	double i_leak = parameters[2];
	
	double v_th = v_on/log(1 + i_on/i_leak);
	double slope = i_leak*derivSatExp((v[0] - v[1])/v_th)/v_th;
	
	j[0] =  slope; j[1] = -slope;
	j[2] = -slope; j[3] =  slope;
}

void dio_updateState(double *parameters, const double *v, double timestep, const double *i){}



const int bjt_terminals_count = 3;
const int bjt_parameters_count = 4;

void bjt_currentCurve(const double *parameters, const double *v, double timestep, double *i){
	double beta     = parameters[0];
	double v_be_on  = parameters[1];
	double i_c_on   = parameters[2];
	double i_c_off  = parameters[3];
	
	double alpha_fwd = beta/(1 + beta);
	double alpha_rev = (0.1*beta)/(1 + 0.1*beta);
	double v_th = v_be_on/log(1 + i_c_on/(alpha_fwd*i_c_off));
	
	double i_ediode = i_c_off*(satExp((v[1] - v[2])/v_th) - 1);
	double i_cdiode = i_c_off*(satExp((v[1] - v[0])/v_th) - 1);
	
	i[0] = alpha_fwd*i_ediode - i_cdiode;
	i[1] = i_ediode*(1 - alpha_fwd) + i_cdiode*(1 - alpha_rev);
	i[2] = alpha_rev*i_cdiode - i_ediode;
}

void bjt_jacobian(const double *parameters, const double *v, double timestep, double *j){
	double beta     = parameters[0];
	double v_be_on  = parameters[1];
	double i_c_on   = parameters[2];
	double i_c_off  = parameters[3];
	
	double alpha_fwd = beta/(1 + beta);
	double alpha_rev = (0.25*beta)/(1 + 0.25*beta);
	double v_th = v_be_on/log(1 + i_c_on/(alpha_fwd*i_c_off));
	
	double i_ediode_d_dvb = i_c_off*derivSatExp((v[1] - v[2])/v_th)/v_th;
	double i_cdiode_d_dvb = i_c_off*derivSatExp((v[1] - v[0])/v_th)/v_th;
	double i_ediode_d_dve = -i_ediode_d_dvb;
	double i_cdiode_d_dvc = -i_cdiode_d_dvb;
	double i_cdiode_d_dve = 0;
	double i_ediode_d_dvc = 0;
	
	// ic derivatives w.r.t c, b, e
	j[0] = alpha_fwd*i_ediode_d_dvc - i_cdiode_d_dvc;
	j[1] = alpha_fwd*i_ediode_d_dvb - i_cdiode_d_dvb;
	j[2] = alpha_fwd*i_ediode_d_dve - i_cdiode_d_dve;

	// ib derivatives w.r.t. c,b,e
	j[3] = i_ediode_d_dvc*(1 - alpha_fwd) + i_cdiode_d_dvc*(1 - alpha_rev);
	j[4] = i_ediode_d_dvb*(1 - alpha_fwd) + i_cdiode_d_dvb*(1 - alpha_rev);	
	j[5] = i_ediode_d_dve*(1 - alpha_fwd) + i_cdiode_d_dve*(1 - alpha_rev);
	
	// ie derivatives w.r.t. c,b,e
	j[6] = alpha_rev*i_cdiode_d_dvc - i_ediode_d_dvc;
	j[7] = alpha_rev*i_cdiode_d_dvb - i_ediode_d_dvb;
	j[8] = alpha_rev*i_cdiode_d_dve - i_ediode_d_dve;
}

void bjt_updateState(double *parameters, const double *v, double timestep, const double *i){}
