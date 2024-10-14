#include "fuzzyintervalbycuts.hpp"
#include <iostream>
#include <string>
#include <cmath>



using namespace std;



FinteralphaCut::FinteralphaCut(void)
{
	cut=NULL;
	alpha = NULL;
}

void FinteralphaCut::dynalloc(size_t nb_cuts)	{
	cut = gsl_matrix_alloc (nb_cuts, 2);
	alpha = (double *) malloc( nb_cuts * sizeof(double) );
}


void FinteralphaCut::free_alphaCut() {

	gsl_matrix_free(cut);
	free(alpha);
}



size_t FinteralphaCut::getnbCuts()	{
	return (size_t)(cut->size1);
}


gsl_vector * FinteralphaCut::get_cut(size_t cut_index) 	{
	gsl_vector * v  = gsl_vector_alloc(2);
	gsl_matrix_get_row(v, cut, cut_index);
	return v;
}


void FinteralphaCut::set_cut(const gsl_vector * the_new_cut, size_t cut_index)	{
	gsl_matrix_set_row (cut, cut_index, the_new_cut);
}
	
void FinteralphaCut::set_alpha(double alp, size_t cut_index)	{
	alpha[cut_index]=alp;
}

void FinteralphaCut::set_cut(const gsl_vector * the_new_cut, double alp, size_t cut_index)	{
	set_cut(the_new_cut,cut_index);
	set_alpha(alp,cut_index);
}
	
void FinteralphaCut::save_in(char * filename)	{

}


void FinteralphaCut::display()	{
	size_t nb_cuts = getnbCuts();
	for(size_t i =0;i<nb_cuts;i++)	{
		cout << "cuts " <<  alpha[i] << " : (" << gsl_matrix_get(cut,i,0) << " " << gsl_matrix_get(cut,i,1) <<")" << endl;
	}
	
}
