/*
 * Author: Joslyn Renfrey
 * Date: 13/01/2023
 */

#include"circuitsim.h"
#include<string.h>
#include<stdio.h>
#include<stdlib.h>

int main(int argc, char *argv[]){
	
	if(argc < 2){ return 0; }
	
	char *spec = argv[1], *results = NULL;
	if(argc == 2){
		results = malloc(strlen(spec) + 20);
		strcpy(results, spec);
		strcpy(results + strlen(spec), "_results.csv");
	} else {
		results = argv[2];
	}
	
	FILE *spec_f = fopen(spec, "r");
	FILE *results_f = fopen(results, "w");
	
	if(spec_f == NULL){
		fprintf(stderr, "error: could not open file \"%s\"\r\n", spec);
		return -1;
	}
	if(results_f == NULL){
		fprintf(stderr, "error: could not open file \"%s\"\r\n", results);
		return -1;
	}
	
	sim_t s;
	
	if(!parseFile(spec_f, &s)){
		return -1;
	}
	//printf("file parsed\n");
	fclose(spec_f);
	if(!simulate(&s, results_f)){
		return -1;
	}
	fclose(results_f);
	return 0;
}
