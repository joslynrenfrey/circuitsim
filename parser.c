/*
 * Author: Joslyn Renfrey
 * Date: 13/01/2023
 */

#include "circuitsim.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include<stdlib.h>

static int isendline(int c){ return (c == '\n') || (c == '\r') || (c == '\v'); }
static int isblank(int c){ return (c == ' ') || (c == '\t'); }
static int isnonblank(int c){ return !isendline(c) && !isblank(c) && c >= 0; }

static int getWord(FILE *f, char *dest){
	int c, i;
	// skip blank space, but not line end
	while(isblank(c = fgetc(f)));
	if(c < 0){ return 0; }
	ungetc(c, f);
	// if we could not get any more words for this line
	if(isendline(c)){ return 0; }
	for(i = 0; isnonblank(c = fgetc(f)); i += (i <= max_name_len)? 1 : 0){
		dest[i] = c;
	}
	if(c >= 0){ ungetc(c, f); }
	dest[i] = '\0';
	return 1;
}

static double getDouble(FILE *f){
	char buffer[max_name_len + 1];
	if(!getWord(f, buffer)){ return NAN; }
	char *endptr = NULL;
	double d = strtod(buffer, &endptr);
	if(endptr == buffer){ return NAN; }
	switch(*endptr){
		case 'T': d *=  1e12 ; break;
		case 'G': d *=  1e9  ; break;
		case 'M': d *=  1e6  ; break;
		case 'k': d *=  1e3  ; break;
		case 'm': d *=  1e-3 ; break;
		case 'u': d *=  1e-6 ; break;
		case 'n': d *=  1e-9 ; break;
		case 'p': d *=  1e-12; break;
	}
	return d;
}

// finishes a line regardless of if there are any words remaining
// returns true if there are remaining lines
static uint8_t getNextLine(FILE *f){
	int c1;
	while(!isendline(c1 = fgetc(f))){
		if(c1 < 0){ return 0; }
	}
	//check next character...
	int c2 = fgetc(f);
	// EOF: no more lines
	if(c2 < 0){
		return 0;
	}
	// non-endline: more line(s)
	else if(!isendline(c2)){
		ungetc(c2, f);
		return 1;
	}
	// unequal line-enders: a line ending pair
	// peek ahead to see if EOF and hence, whether there is a line ahead
	else if(c1 != c2){
		c2 = fgetc(f);
		ungetc(c2, f);
		return c2 >= 0;
	}
	// equal line-enders: empty line ahead
	else {
		ungetc(c2, f);
		return 1;
	}
}

static component_t *componentByName(const char *s, component_t *c, int n){
	for(int i = 0; i < n; i++){
		if(strcmp(s, c[i].name) == 0){ return c + i; }
	}
	return NULL;
}

static node_t *nodeByName(const char *s, node_t *c, int n){
	for(int i = 0; i < n; i++){
		if(strcmp(s, c[i].name) == 0){ return c + i; }
	}
	return NULL;
}

static int nodeIndexByName(const char *s, node_t *c, int n){
	for(int i = 0; i < n; i++){
		if(strcmp(s, c[i].name) == 0){ return i; }
	}
	return -1;
}

#define COMPONENT_EXTERN( n ) \
	extern const int n##_terminals_count; \
	extern const int n##_parameters_count; \
	extern void n##_currentCurve(const double *parameters, const double *v, double timestep, double *i); \
	extern void n##_jacobian(const double *parameters, const double *v, double timestep, double *j); \
	extern void n##_updateState(double *parameters, const double *v, double timestep, const double *i);
COMPONENT_LIST(COMPONENT_EXTERN)

int parseFile(FILE *f, sim_t *s){
	#define COMPONENT_PROTOTYPE( n ) { \
		.name = #n, \
		.terminals_count = n##_terminals_count, \
		.parameters_count = n##_parameters_count, \
		.currentCurve = &n##_currentCurve, \
		.jacobian = &n##_jacobian, \
		.updateState = &n##_updateState \
	},
	component_t component_prototypes[] = {COMPONENT_LIST(COMPONENT_PROTOTYPE) {.name = ""}};

	s->n_count = 0; s->c_count = 0;
	int n_space = 0, c_space = 0;

	s->errorsq = default_errorsq;
	s->convrate = default_convrate;
	s->maxiter = default_maxiter;
	
	#define ERROR(condition, ...) \
	if(condition){ \
		fprintf(stderr, "error: line %i: ", line_num); \
		fprintf(stderr, __VA_ARGS__); \
		fprintf(stderr, "\r\n"); \
		return 0; \
	}

	// we will perform 3 passes
	long f_start = ftell(f);
	
	char word[max_name_len + 1];
	int line_num = 0;
	// 1st pass: determine an upper bound on the space we need to store components and nodes
	do {
		line_num++;
		// if line is empty, skip
		if(!getWord(f, word)){ continue; }
		else if (word[0] == '#' || word[0] == '/'){ continue; }
		
		// count node space
		else if(strcmp(word, "nodes") == 0){
			while(getWord(f, word)){ n_space++; }
		}
		// count component space
		else {
			c_space++;
		}
	}while(getNextLine(f));
	
	s->c = malloc(sizeof(component_t)*c_space);
	s->n = malloc(sizeof(node_t)*n_space);
	memset(s->c, 0, sizeof(component_t)*c_space);
	memset(s->n, 0, sizeof(node_t)*n_space);
	
	// 2nd pass: read all nodes, and set voltages and measure flags
	// after this we want to sort them into variable and fixed nodes
	fseek(f, f_start, SEEK_SET);
	line_num = 0;
	do {
		line_num++;
		// if line is empty, skip
		if(!getWord(f, word)){ continue; }
		else if (word[0] == '#' || word[0] == '/'){ continue; }
				
		// also handle the timestep and endtime settings
		else if(strcmp(word, "timestep") == 0){
			s->timestep = getDouble(f);
			ERROR(isnan(s->timestep) || s->timestep <= 0, "timestep invalid");
		}	
		else if(strcmp(word, "endtime") == 0){
			s->endtime = getDouble(f);
			ERROR(isnan(s->endtime) || s->endtime <= 0, "endtime invalid");
		}
		else if(strcmp(word, "convrate") == 0){
			s->convrate = getDouble(f)/100;
			ERROR(isnan(s->convrate) || s->convrate <= 0, "convrate invalid");
		}	
		else if(strcmp(word, "errorsq") == 0){
			s->errorsq = getDouble(f);
			ERROR(isnan(s->errorsq) || s->errorsq <= 0, "errorsq invalid");
		}
		else if(strcmp(word, "maxiter") == 0){
			double d = getDouble(f);
			ERROR(isnan(d) || d <= 0, "maxiter invalid");
			s->maxiter = d;
		}
		
		// nodes: add nodes
		else if(strcmp(word, "nodes") == 0){
			while(getWord(f, word)){
				strcpy((s->n + s->n_count)->name, word);
				(s->n + s->n_count)->is_fixed = 0;
				(s->n + s->n_count)->is_measured = 0;
				(s->n_count)++;
			}
		}
		
		// set: make fixed nodes. list is pairs of node names and voltages
		else if(strcmp(word, "set") == 0){
			while(getWord(f, word)){
				node_t *target = nodeByName(word, s->n, s->n_count);
				ERROR(target == NULL, "unrecognised node \"%s\"", word);
				target->is_fixed = 1;
				target->fixed_voltage = getDouble(f);
				ERROR(isnan(target->fixed_voltage), "invalid node voltage");
			}
		}
	} while(getNextLine(f));
	
	// we want to sort the nodes into variable and fixed
	// before we parse components so that the node indices
	// are properly applied
	int var_n_count = s->n_count;
	for(int i = 0; i < var_n_count;){
		if(!(s->n)[i].is_fixed){ i++; }
		else {
			node_t temp = (s->n)[i];
			(s->n)[i] = (s->n)[var_n_count - 1];
			(s->n)[var_n_count - 1] = temp;
			var_n_count--;
		}
	}
	
	// third pass: read component related nodes, 
	fseek(f, f_start, SEEK_SET);
	line_num = 0;
	do {
		line_num++;
		// if line is empty, skip
		if(!getWord(f, word)){ continue; }
		else if (word[0] == '#' || word[0] == '/'){ continue; }
		
		// nodes: ignore all commands related to nodes and time settings
		else if(strcmp(word, "timestep") == 0){}
		else if(strcmp(word, "endtime") == 0){}
		else if(strcmp(word, "nodes") == 0){}
		else if(strcmp(word, "set") == 0){}
		else if(strcmp(word, "convrate") == 0){}
		else if(strcmp(word, "errorsq") == 0){}
		else if(strcmp(word, "maxiter") == 0){}
		
		// measure: enable the is_measured flag on selected nodes or components
		else if(strcmp(word, "measure") == 0){
			while(getWord(f, word)){
				node_t *tn; component_t *tc;
				if((tn = nodeByName(word, s->n, s->n_count)) != NULL){
					tn->is_measured = 1;
				}
				else if((tc = componentByName(word, s->c, s->c_count)) != NULL){	
					tc->is_measured = 1;
				}
				else { ERROR(1, "unrecognised node or component \"%s\"", word); }
			}
		}

		// otherwise assume first word is a component type
		else {
			// search through list of prototypes...
			component_t *t = component_prototypes;
			while(t->name[0] != '\0'){
				if(strcmp(t->name, word) == 0){ break; } t++;
			}
			ERROR(t->name[0] == '\0', "unrecognised component type \"%s\"", word);
			
			memcpy(s->c + s->c_count, t, sizeof(component_t));
			t = s->c + s->c_count; (s->c_count)++;
			t->is_measured = 0;
			
			// make sure we get the components name...
			ERROR(!getWord(f, word), "expected component name");
			strcpy(t->name, word);
			
			// match terminal node names to indices
			for(int i = 0; i < t->terminals_count; i++){
				ERROR(!getWord(f, word), "expected terminal node");
				t->terminals[i] = nodeIndexByName(word, s->n, s->n_count);
				ERROR(t->terminals[i] < 0, "unrecognised node \"%s\"", word);
			}
			
			// read parameters
			for(int i = 0; i < t->parameters_count; i++){
				t->parameters[i] = getDouble(f);
				ERROR(isnan(t->parameters[i]), "expected numerical parameter");
			}
		}
	} while(getNextLine(f));
	return 1;
}
