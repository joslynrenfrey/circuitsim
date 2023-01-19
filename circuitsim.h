/*
 * Author: Joslyn Renfrey
 * Date: 13/01/2023
 */

#ifndef _CIRCUITSIM_H
#define _CIRCUITSIM_H

#include<stddef.h>
#include<stdio.h>

#define max_name_len 20
#define max_terms 5
#define max_params 5

#define default_maxiter 10000
#define default_errorsq 1e-18
#define default_convrate 0.8

#define COMPONENT_LIST( X ) X(res) X(src) X(ind) X(cap) X(dio) X(bjt)

typedef struct {
	char name[max_name_len + 1];
	
	int terminals_count; int terminals[max_params];
	int parameters_count; double parameters[max_params];
	
	uint8_t is_measured;
	
	void (*currentCurve)(const double *parameters, const double *v, double timestep, double *i);
	void (*jacobian)(const double *parameters, const double *v, double timestep, double *j);
	void (*updateState)(double *parameters, const double *v, double timestep, const double *i);
} component_t;

typedef struct {
	char name[max_name_len + 1];
	
	uint8_t is_fixed;
	float fixed_voltage;
	uint8_t is_measured;
} node_t;

typedef struct {
	double errorsq, convrate;
	double timestep, endtime;
	int maxiter;
	
	int c_count, n_count;
	component_t *c;
	node_t *n;
} sim_t;

int parseFile(FILE *f, sim_t *s);
int simulate(sim_t *s, FILE *f);

#endif
