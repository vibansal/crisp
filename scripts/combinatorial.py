#! /usr/bin/env python
# AUTHOR VIKAS BANSAL  
## code fixed for FET p-values, checked feb 6 2015 working correctly 

## functions for calculating fisher's exact test p-value

import sys, os, glob, string, subprocess,time, math, random


def ncr(n,r):
        ll = 0;
        for i in xrange(min(r,n-r)): ll += math.log(float(n-i)/(i+1),10);
        return ll;

def ncrbase2(n,r):
        ll = 0;
        for i in xrange(min(r,n-r)): ll += math.log(float(n-i)/(i+1),2);
        return ll;


def mcr(n,rlist):
        ll =0; k=n;
        for i in xrange(len(rlist)):
                for j in xrange(rlist[i]): ll += math.log(k,10)-math.log(j+1,10); k -=1;
        return ll;

#from scipy.stats import chi2  # chi square test p-value 

def fet(c1,r1,c2,r2):  # c1 <= r1 and c2 <= r2
	minc =0; maxc = c1+c2;
	if r2 < c1+c2: minc = c1+c2-r2;
	if r1 < c1+c2: maxc = r1; 
	ptable = ncr(r1,c1) + ncr(r2,c2); c11 = minc; c21 = c1+c2-minc; 
	p1 = ncr(r1,c11); p2 =  ncr(r2,c21);
	pvalue0 = 0.0; pvalue1 =0.0; pvalue = 0.0;
	flags = [0,0,0];

	## BUG FIXED 
	## sept 7 2012 bug  in code, if we don't enter the code below or pvalue0 are not updated, they remain at very low values 
	## if the <= comparison doesn't work then also there could be incorrect value
	while c11 <= maxc:

		if p1+p2 < ptable + 1.0e-10:
			if flags[2] ==0: pvalue = p1+p2; flags[2] = 1;  
			elif p1 + p2 <= pvalue: pvalue += math.log(1+math.pow(10,p1+p2-pvalue),10); 
			else: pvalue = p1+p2 + math.log(1+math.pow(10,pvalue-p1-p2),10);
		
		if c11 <= c1:
			if flags[0] ==0: pvalue0 = p1+p2; flags[0] = 1; 
			elif p1 + p2 <= pvalue0: pvalue0 += math.log(1+math.pow(10,p1+p2-pvalue0),10); 
			else: pvalue0 = p1+p2 + math.log(1+math.pow(10,pvalue0-p1-p2),10);
		
		if c11 >= c1:
			if flags[1] ==0: pvalue1 = p1+p2; flags[1] = 1; 
			elif p1 + p2 <= pvalue1: pvalue1 += math.log(1+math.pow(10,p1+p2-pvalue1),10); 
			else: pvalue1 = p1+p2 + math.log(1+math.pow(10,pvalue1-p1-p2),10);

		if r1-c11 > 0: p1 += math.log(r1-c11,10) - math.log(c11+1,10);
		if c21 > 0: p2 += math.log(c21,10) - math.log(r2-c21+1,10);
		c11 +=1; c21 -= 1;

#	if pvalue0 <= pvalue1: pvalue = pvalue1 + math.log(1+math.pow(10,pvalue0-pvalue1),10);
#	else: pvalue = pvalue0 + math.log(1+math.pow(10,pvalue1-pvalue0),10);
#	if ptable <= pvalue0: pvalue0 += math.log(1+math.pow(10,ptable-pvalue0),10);
#	else: pvalue0 = ptable + math.log(1+math.pow(10,pvalue0-ptable),10);
	if flags[2] == 0: pvalue += ptable; 
	#if flags[1] == 0: pvalue1 = ptable; 
	#if flags[0] == 0: pvalue0 = ptable; 
	sum = ncr(r1+r2,c1+c2); pvalue -= sum; pvalue0 -= sum; pvalue1 -= sum;
	if pvalue >= -0.001: pvalue =0;
	if pvalue0 >= -0.001: pvalue0 =0;
	if pvalue1 >= -0.001: pvalue1 =0;
	return [pvalue,pvalue0,pvalue1];


#P =  fet(1,1200,0,1000); print P,math.pow(10,P[0]),math.pow(10,P[1]),math.pow(10,P[2]);
