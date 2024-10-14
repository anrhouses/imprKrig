#include "fuzzyinterval.hpp"
#include "fuzzyintervalbycuts.hpp"
#include <iostream>
#include <cstring>
#include <cmath>

using namespace std;


double triangular (double x) {
	double tmp =0;
	if(x>-1 && x<1)	{
		tmp = 1.0-fabs(x);
	}
	return tmp;
}

double modifiedEpan (double x) {
	double tmp =0;
	if(x>-1 && x<1)
		tmp = 1-x*x;
	return tmp;
}




Finter::Finter(void)
{
	membership = NULL;
	mode = 0;
	delta = 1;
}

void Finter::initialize(const char * memberFuncName, double mod, double del)	{
	if(del<=0)
		cerr << "Delta must be positive.\n";
		
	if(strcmp(memberFuncName,"TRI")==0)
		membership = triangular;
	else if(strcmp(memberFuncName,"EPA")==0)
		membership = modifiedEpan;
		
	mode = mod;
	delta = del;
}

double Finter::membershipdegree(double position)
{
	return membership(fabs(position-mode)/delta);
}



/*void Finter::extract_cuts(alphaCut* alcuts, size_t n_cuts)	{
	double alp;
	gsl_vector* bounds = gsl_vector_alloc(2); 
	//mexPrintf("inside func\n");
	for (size_t i=0;i<n_cuts;i++)
	{
		alp=(double)i/(n_cuts-1);
		extract_cut(bounds, alp);		
		//mexPrintf("cut %g : (%g %g)\n",alp, gsl_vector_get(bounds,0),gsl_vector_get(bounds,1));
		alcuts[i].cut  = gsl_vector_alloc(2); 
		gsl_vector_set(alcuts[i].cut,0,gsl_vector_get(bounds,0));
		gsl_vector_set(alcuts[i].cut,1,gsl_vector_get(bounds,1));
		alcuts[i].alpha = alp;
	}
	gsl_vector_free(bounds);
}*/


void Finter::to_FinteralphaCut(FinteralphaCut* alcuts)	{
	double alp;
	size_t n_cuts = alcuts->getnbCuts();
	gsl_vector* bounds = gsl_vector_alloc(2); 
	for (size_t i=0;i<n_cuts;i++)
	{
		alp=(double)i/(n_cuts-1);
		extract_cut(bounds, alp);		
		//mexPrintf("cut %g : (%g %g)\n",alp, gsl_vector_get(bounds,0),gsl_vector_get(bounds,1));
		alcuts->set_cut(bounds,alp,i);
	}
	gsl_vector_free(bounds);
}


/*void Finter::extract_cuts_leftright(gsl_vector* kernel, alphaCut* alcutsleft,alphaCut* alcutsright, size_t n_cuts)	{
	double alpI,alpJ;
	gsl_vector* boundsI = gsl_vector_alloc(2); 
	gsl_vector* boundsJ = gsl_vector_alloc(2); 
	//kernel  = gsl_vector_alloc(2);
	alpI=1.0;
	extract_cut(boundsI, alpI);
	gsl_vector_set(kernel,0,gsl_vector_get(boundsI,0));
	gsl_vector_set(kernel,1,gsl_vector_get(boundsI,1));
	for (int i=n_cuts-2;i>-1;i--)
	{
		alpJ=(double)i/(n_cuts-1);
		extract_cut(boundsJ, alpJ);		
		alcutsleft[i].cut  = gsl_vector_alloc(2); 
		gsl_vector_set(alcutsleft[i].cut,0,gsl_vector_get(boundsJ,0));
		gsl_vector_set(alcutsleft[i].cut,1,gsl_vector_get(boundsI,0));
		alcutsleft[i].alpha = alpJ;
		alcutsright[i].cut  = gsl_vector_alloc(2); 
		gsl_vector_set(alcutsright[i].cut,0,gsl_vector_get(boundsI,1));
		gsl_vector_set(alcutsright[i].cut,1,gsl_vector_get(boundsJ,1));
		alcutsright[i].alpha = alpJ;
		
		gsl_vector_set(boundsI,0,gsl_vector_get(boundsJ,0));
		gsl_vector_set(boundsI,1,gsl_vector_get(boundsJ,1));
		alpI=alpJ;
	}
	gsl_vector_free(boundsI);
	gsl_vector_free(boundsJ);
}
*/

void Finter::to_FinteralphaCut_leftright(gsl_vector* kernel, FinteralphaCut* alcutsleft, FinteralphaCut* alcutsright)	{
	double alpI,alpJ;
	gsl_vector* boundsI = gsl_vector_alloc(2); 
	gsl_vector* boundsJ = gsl_vector_alloc(2); 
	gsl_vector* temp_vec = gsl_vector_alloc(2); 
	size_t n_cuts = alcutsleft->getnbCuts()+1;
	alpI=1.0;
	extract_cut(boundsI, alpI);
	gsl_vector_memcpy(kernel,boundsI);
	
	for (int i=n_cuts-2;i>-1;i--)
	{
		alpJ=(double)i/(n_cuts-1);
		extract_cut(boundsJ, alpJ);		
		
		gsl_vector_set(temp_vec,0,gsl_vector_get(boundsJ,0));
		gsl_vector_set(temp_vec,1,gsl_vector_get(boundsI,0));		
		alcutsleft->set_cut(temp_vec,alpJ,i);
		
		gsl_vector_set(temp_vec,0,gsl_vector_get(boundsI,1));
		gsl_vector_set(temp_vec,1,gsl_vector_get(boundsJ,1));
		alcutsright->set_cut(temp_vec,alpJ,i);
		
		gsl_vector_memcpy(boundsI,boundsJ);
		
		alpI=alpJ;
	}
	gsl_vector_free(temp_vec);
	gsl_vector_free(boundsI);
	gsl_vector_free(boundsJ);
}


void Finter::extract_cut(gsl_vector* bounds, double alpha)	{
	if(alpha<0 || alpha > 1)
		cerr << "must be between 0 and 1.\n";
	double supBound, infBound = mode-delta;
	size_t nbSteps = 10000000 , i=0;
	do	{
		infBound += delta/(nbSteps-1);
		i++;
	}
	while (i<nbSteps && membershipdegree(infBound)<alpha);
	supBound = mode+mode-infBound;
	gsl_vector_set(bounds,0,infBound);
	gsl_vector_set(bounds,1,supBound);
}

