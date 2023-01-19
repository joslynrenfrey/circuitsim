/*
 * Author: Joslyn Renfrey
 * Date: 13/01/2023
 */






#include"circuitsim.h"
#include<math.h>
#include<stdio.h>
#include <string.h>
#include<stdlib.h>

void vecSub(int n, double *x, double *y, double *r){
	for(int i = 0; i < n; i++){ r[i] = x[i] - y[i]; }
}

double vecDot(int n, double *x, double *y){
	double r = 0;
	for(int i = 0; i < n; i++){ r += x[i]*y[i]; }
	return r;
}

void vecScale(int n, double *x, double s, double *y){
	for(int i = 0; i < n; i++){ y[i] = s*x[i]; }
}

int solveLinear(int n, int rowskip, int *swap_indices, double *A, double *b){
	#define swap(x, y) {double temp = x; x = y; y = temp; }
	// gaussian elimination with full pivoting
	for(int pivot = 0; pivot < n; pivot++){
		// search for the greatest element in the remaining sub matrix
		int swaprow = pivot, swapcol = pivot;
		double swap_mag = fabs(A[rowskip*pivot + pivot]);
		for(int row = pivot; row < n; row++){
			for(int col = pivot; col < n; col++){
				double candidate_mag = fabs(A[rowskip*row + col]);
				if(candidate_mag > swap_mag){
					swap_mag = candidate_mag, swaprow = row, swapcol = col;
				}
			}
		}
		
		// swap the row
		for(int col = pivot; col < n; col++){
			swap(A[rowskip*swaprow + col], A[rowskip*pivot + col]);
		}
		swap(b[swaprow], b[pivot]);
		
		// and swap the column
		for(int row = 0; row < n; row++){
			swap(A[rowskip*row + pivot], A[rowskip*row + swapcol]);
		}
		
		// store the column permutation to be undone later
		swap_indices[pivot] = swapcol;
		
		
		// eliminate pivot column of all lower rows by row scale and subtract
		for(int row = pivot + 1; row < n; row++){
			double scale = A[rowskip*row + pivot]/A[rowskip*pivot + pivot];
			for(int col = pivot; col < n; col++){
				A[rowskip*row + col] -= scale*A[rowskip*pivot + col];
			}
			b[row] -= scale*b[pivot];
		}
		
		if(A[rowskip*pivot + pivot] == 0){ return 0; }
	}
	
	// eliminate upper diagonal
	for(int pivot = n - 1; pivot >= 0; pivot--){		
		// normalize the diagonal element in the pivotvot row by row scaling
		// there are only 2 elements in the pivotvot row
		double scale = 1/A[rowskip*pivot + pivot];
		b[pivot] *= scale; A[rowskip*pivot + pivot] *= scale;

		// for each higher row, eliminate the pivotvot column by row scale and subtract
		// again there are only 2 elements in the pivotvot row.
		for(int row = pivot - 1; row >= 0; row--){
			double scale = A[rowskip*row + pivot];
			b[row] -= scale*b[pivot];
			A[rowskip*row + pivot] -= scale; // A[rowskip*pivot + pivot] = 1;
		}
	}
	
	// undo the column swaps in reverse
	for(int pivot = n - 1; pivot >= 0; pivot--){
		int swapcol = swap_indices[pivot];
		swap(b[swapcol], b[pivot]);
	}
	return 1;
}

static void evalErrorAndJacobian(sim_t *s, double *v, double *e, double *jac){	
	// calculate error vector as the sum of currents at each node,
	// and the jacobian as the rate of change of the that w.r.t node voltage
	memset(e, 0, sizeof(double)*s->n_count);
	memset(jac, 0, sizeof(double)*s->n_count*s->n_count);
	
	for(int i = 0; i < s->c_count; i++){
		// action of G on v. v_term is the fragment of v
		// that only this component's curve operates on
		double v_term[max_terms], i_term[max_terms];
		for(int j = 0; j < s->c[i].terminals_count; j++){
			v_term[j] = v[s->c[i].terminals[j]];
		}
		
		double jac_term[max_terms*max_terms];
		s->c[i].currentCurve(s->c[i].parameters, v_term, s->timestep, i_term);
		s->c[i].jacobian(s->c[i].parameters, v_term, s->timestep, jac_term);
		
		// i_term is the corresponding fraction of F(G v)
		// linearly combine GT i_term from each component
		// to get the complete e = GT F(G v)
		for(int j = 0; j < s->c[i].terminals_count; j++){
			e[s->c[i].terminals[j]] += i_term[j];
		}
		
		for(int col = 0; col < s->c[i].terminals_count; col++){
			for(int row = 0; row < s->c[i].terminals_count; row++){
				int tcol = s->c[i].terminals[col], trow = s->c[i].terminals[row];
				jac[trow*s->n_count + tcol] += jac_term[row*s->c[i].terminals_count + col];
			}
		}
	}
}

static void printLabels(sim_t *s, FILE *f){
	fprintf(f, "time(s)");
	// print labels for measured nodes
	for(int i = 0; i < s->n_count; i++){
		if(s->n[i].is_measured){
			fprintf(f, ", %s(V)", s->n[i].name);
		}
	}
	// print labels for components that are being measured (only supports 2 terminal devices)
	for(int i = 0; i < s->c_count; i++){
		if(s->c[i].is_measured && s->c[i].terminals_count == 2){
			fprintf(f, ", %s(V), %s(A)", s->c[i].name, s->c[i].name);
		}
	}
	fprintf(f, "\n");
}

static void printState(sim_t *s, double time, double *v, FILE *f){
	fprintf(f, "%.6e", time);
	
	// print voltages for measured nodes
	for(int i = 0; i < s->n_count; i++){
		if(s->n[i].is_measured){
			fprintf(f, ", %.6e", v[i]);
		}
	}
	for(int i = 0; i < s->c_count; i++){
		// time is advanced by calling updateState time-sensitive components
		// (capacitors, inductors) with the current voltage and current
		double v_term[max_terms], i_term[max_terms];
		for(int j = 0; j < s->c[i].terminals_count; j++){
			// we need to reconstruct this, since it was clobbed
			v_term[j] = v[s->c[i].terminals[j]];
		}
		s->c[i].currentCurve(s->c[i].parameters, v_term, s->timestep, i_term);
		s->c[i].updateState(s->c[i].parameters, v_term, s->timestep, i_term);
		
		// print voltages and currents for measured components
		if(s->c[i].is_measured && s->c[i].terminals_count == 2){
			fprintf(f, ", %.6e, %.6e", v_term[0] - v_term[1], (i_term[0] - i_term[1])/2);
		}
	}
	fprintf(f, "\n");
}



int simulate(sim_t *s, FILE *f){
	double e_sqmag = 0;
	double *v = malloc(sizeof(double)*s->n_count);
	double *e = malloc(sizeof(double)*s->n_count);
	double *jac = malloc(sizeof(double)*s->n_count*s->n_count);
	int *swap_indices = malloc(sizeof(int)*s->n_count);


	int stats_steps = 0;
	int stats_iters_total = 0;
	int stats_iters_max = 0;
	int stats_iters_min = s->maxiter;

	// assume that nodes are sorted by variable nodes, then fixed nodes
	int var_n_count = 0;
	for(int i = 0; i < s->n_count; i++){
		if(!s->n[i].is_fixed){ var_n_count++; v[i] = 0; }
		else { v[i] = s->n[i].fixed_voltage; }
	}
		
	// the first line will be column labels
	printLabels(s, f);	
	
	for(double time = 0 ; time < s->endtime; time += s->timestep){
		for(int iter = 0; iter < s->maxiter; iter++){
			evalErrorAndJacobian(s, v, e, jac);
			e_sqmag = vecDot(var_n_count, e, e);
			if(e_sqmag < s->errorsq){
				if(stats_iters_max < iter){
					stats_iters_max = iter;
				}
				if(iter != 0 && iter < stats_iters_min){
					stats_iters_min = iter;
				}
				break;
			}
			
			//newton's method
			// multiply inverse jacobian by error vector, result
			// is stored in the error vector...
			if(!solveLinear(var_n_count, s->n_count, swap_indices, jac, e)){
				fprintf(stderr, "error: singular jacobian on time step %.6e, iteration %i\n", time, iter);
				return 0;
			}
			
			vecScale(var_n_count, e, s->convrate, e);
			vecSub(var_n_count, v, e, v);
			stats_iters_total++;
		}
		
		if(e_sqmag < s->errorsq){
			printState(s, time, v, f);
			stats_steps++;
		} else {
			fprintf(stderr, "error: could not converge at timestep %.6e\n", time);
			fprintf(stderr, "error: minimum E^2 = %.6g\n", e_sqmag);
			
			return 0;
		}
	}
	
	fprintf(stderr, "cycles = %i, total iterations = %i\n", stats_steps, stats_iters_total);
	fprintf(stderr, "avg iterations/cycle = %.1f\n", (double) stats_iters_total/(double) stats_steps);
	fprintf(stderr, "min iterations/cycle = %i\n", stats_iters_min);
	fprintf(stderr, "max iterations/cycle = %i\n", stats_iters_max);
	
	return 1;
}

	

