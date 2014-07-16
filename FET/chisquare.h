
#ifndef INC_chisquare_H
#define INC_chisquare_H
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define ITMAX 100 // maximum number of iterations for calculating the gamma functions....
#define TINY 1e-200 // smallest p-value for chi-square
//extern int PIVOTSAMPLE; 

int chi2pvalue(double** ctable, int* strata, int size,double* chisqpvalue,double chistatistics[]);

double kf_gammaq(double DOF,double val); // for kfunc.c 

#endif
