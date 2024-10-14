#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <cmath>
#include <cfloat>
#include <sstream>


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sort.h>

#include "utils.hpp"


using namespace std;



char * convertInt(int number)	{
   	
	char * converted_char_array = (char*) malloc (100);
	stringstream ss;//create a stringstream
   	ss << setprecision(4) << number;//add number to the stream
   	string st = ss.str();
   	converted_char_array = const_cast<char*> ( st.c_str() );//return a string with the contents of the stream
	return converted_char_array;
	
}

char * convertDouble(double number)	{
   	
	char * converted_char_array = (char*) malloc (100);
	stringstream ss;//create a stringstream
   	ss << number;//add number to the stream
   	string st = ss.str();
   	converted_char_array = const_cast<char*> ( st.c_str() );//return a string with the contents of the stream
	return converted_char_array;
	
}
	
void initOutFile(ofstream * outFile, const char * currentFileName) {

	(*outFile).open(currentFileName);
	if(!(*outFile)) { // file couldn't be opened
		cout << currentFileName << ": " ;
		cout << "ERROR- file could not be opened" << endl;
		}

}

void initInFile(ifstream * inFile, const char * currentFileName) {
	
	(*inFile).open(currentFileName);
	if(!(*inFile)) { // file couldn't be opened
		cout << currentFileName << ": " ;
		cout << "ERROR- file could not be opened" << endl;
		}
	(*inFile).seekg(0);

}


size_t countData(const char * dataFileName)	{ /* Count of the sample size */
	
	ifstream dataFile;
	char str[2000];
	size_t num_samples = 0;
	initInFile(&dataFile,dataFileName);
	while(!dataFile.eof())
	{
		  dataFile.getline(str,2000);
		  num_samples++;
	}
	return num_samples;

}





double sum_points(const gsl_vector *weights, const gsl_vector *data)	{
	
	double result = gsl_vector_get(weights,0)*gsl_vector_get(data,0) ;
	size_t N = weights->size;
	
	for (size_t i=1;i<N;i++)	{
		result += gsl_vector_get(weights,i)*gsl_vector_get(data,i) ;
	}
	return result;

}

double get_sum_bound(const gsl_vector *Z,const gsl_vector *inf_weights,const gsl_vector *sup_weights, bool isInf)	{

	double sum = 0;
	double L = 1;
	size_t k=0;
	size_t N = Z->size;

	gsl_vector * L_tocompute  = gsl_vector_alloc(N);
	gsl_vector * Z_tocompute  = gsl_vector_alloc(N);

	gsl_vector * inf_L  = gsl_vector_alloc(N);
	gsl_vector * sup_L  = gsl_vector_alloc(N);

	gsl_permutation * perm = gsl_permutation_alloc(N);

	gsl_vector_memcpy(Z_tocompute,Z);
	gsl_vector_memcpy(inf_L,inf_weights);
	gsl_vector_memcpy(sup_L,sup_weights);

	/*
	cout <<  "inf data and weights non permutes" << endl ;
	for (size_t j=0 ; j<N ; j++){
		cout << "( " << gsl_vector_get(inf_L,j) << " , " << gsl_vector_get(sup_L,j) << " )  ; " << gsl_vector_get(Z_tocompute,j)  << endl;
	}
	cout << endl;
	*/
	
	gsl_sort_vector_index (perm, Z);
	gsl_permute_vector(perm,Z_tocompute);
	gsl_permute_vector(perm,inf_L);
	gsl_permute_vector(perm,sup_L);

	/*
	cout  << "inf data and weights ordered" << endl ;
	for (size_t j=0 ; j<N ; j++){
		cout << "( " << gsl_vector_get(inf_L,j) << " , " << gsl_vector_get(sup_L,j) << " )  ; " << gsl_vector_get(Z_tocompute,j)  << endl;
	}	
	cout << endl;
	*/

	if(isInf)	{
		//cout << "isInf" << endl;
		do {
			L=1;
			if( k==0)	{
				for (size_t i=1 ; i<N ; i++){
					L -= gsl_vector_get(inf_L,i);
				}
			}
			else if (k==N-1)	{
				for (size_t i=0 ; i<N-1 ; i++){
					L -= gsl_vector_get(sup_L,i);
				}
			}
			else	{
				for (size_t i=0 ; i<k ; i++){
					L -= gsl_vector_get(sup_L,i);
				}
				for (size_t i=k+1 ; i<N ; i++){
					L -= gsl_vector_get(inf_L,i);
				}
			}
			k++;
			//cout << "L" << L << endl;
		}
		while (!(L>=gsl_vector_get(inf_L,k-1) && L<=gsl_vector_get(sup_L,k-1)) && k<N);

		if(L>=gsl_vector_get(inf_L,k-1) && L<=gsl_vector_get(sup_L,k-1))	{
			gsl_vector_set(L_tocompute,k-1,L);
			for (size_t i=0 ; i<k-1 ; i++){
				gsl_vector_set(L_tocompute,i,gsl_vector_get(sup_L,i));
			}
			for (size_t i=k ; i<N ; i++){
				gsl_vector_set(L_tocompute,i,gsl_vector_get(inf_L,i));
			}
		}
	}
	else	{
		//cout << "isSup" << endl;
		do {
			L=1;
			if( k==0)	{
				for (size_t i=1 ; i<N ; i++){
					L -= gsl_vector_get(sup_L,i);
				}
			}
			else if (k==N-1)	{
				for (size_t i=0 ; i<N-1 ; i++){
					L -= gsl_vector_get(inf_L,i);
				}
			}
			else	{
				for (size_t i=0 ; i<k ; i++){
					L -= gsl_vector_get(inf_L,i);
				}
				for (size_t i=k+1 ; i<N ; i++){
					L -= gsl_vector_get(sup_L,i);
				}
			}
			k++;
			//cout << "L" << L << endl;
		}
		while (!(L>=gsl_vector_get(inf_L,k-1) && L<=gsl_vector_get(sup_L,k-1)) && k<N);

		if(L>=gsl_vector_get(inf_L,k-1) && L<=gsl_vector_get(sup_L,k-1))	{
			gsl_vector_set(L_tocompute,k-1,L);
			for (size_t i=0 ; i<k-1 ; i++){
				gsl_vector_set(L_tocompute,i,gsl_vector_get(inf_L,i));
			}
			for (size_t i=k ; i<N ; i++){
				gsl_vector_set(L_tocompute,i,gsl_vector_get(sup_L,i));
			}
		}
	}
	sum = sum_points(L_tocompute,Z_tocompute);

	gsl_vector_free(L_tocompute);
	gsl_vector_free(Z_tocompute);
	gsl_vector_free(inf_L);
	gsl_vector_free(sup_L);


	return sum;

}

void sum_intervals_under_constraint(gsl_vector * result, const gsl_vector *inf_weights, const gsl_vector *sup_weights, const gsl_vector *inf_data, const gsl_vector *sup_data)	{

	size_t N = inf_data->size;
	size_t nb_overlapped=0;

	double sum_inf = DBL_MAX ,sum_sup =0 ;
	double temp;
	gsl_vector * inf_Z  = gsl_vector_alloc(N);
	gsl_vector * sup_Z  = gsl_vector_alloc(N);


	for (size_t i=0;i<N;i++)	{
		if(gsl_vector_get(inf_weights,i)>=0)	{
			gsl_vector_set(inf_Z,i,gsl_vector_get(inf_data,i));
			gsl_vector_set(sup_Z,i,gsl_vector_get(sup_data,i));
		} 
		else if(gsl_vector_get(sup_weights,i)<=0)	{
			gsl_vector_set(inf_Z,i,gsl_vector_get(sup_data,i));
			gsl_vector_set(sup_Z,i,gsl_vector_get(inf_data,i));
		}
		else{
			nb_overlapped ++;
		}
	}

	if(nb_overlapped>0)	{
	
		gsl_combination * comb ;
		size_t i_over = 0;
		size_t * overlapped_L_indices = (size_t*) calloc (nb_overlapped,sizeof(size_t)) ;
		size_t * inverted_indices = (size_t*) calloc (nb_overlapped,sizeof(size_t)) ;

		for (size_t i=0;i<N;i++)	{
		if (gsl_vector_get(sup_weights,i)>0 && gsl_vector_get(inf_weights,i)<0) {
				overlapped_L_indices[i_over] = i;
				i_over++;
			}
		}

		for (size_t i=0;i<=nb_overlapped;i++)	{
			comb = gsl_combination_calloc (nb_overlapped, i);
			do
			{
				inverted_indices = gsl_combination_data (comb);
				for (size_t j=0 ; j<nb_overlapped ; j++)	{
					gsl_vector_set(inf_Z,overlapped_L_indices[j],gsl_vector_get(inf_data,overlapped_L_indices[j]));
					gsl_vector_set(sup_Z,overlapped_L_indices[j],gsl_vector_get(sup_data,overlapped_L_indices[j]));
				}	
				for (size_t j=0 ; j<i ; j++){
					gsl_vector_set(inf_Z,overlapped_L_indices[inverted_indices[j]],gsl_vector_get(sup_data,overlapped_L_indices[inverted_indices[j]]));
					gsl_vector_set(sup_Z,overlapped_L_indices[inverted_indices[j]],gsl_vector_get(inf_data,overlapped_L_indices[inverted_indices[j]]));
				}	
			
				temp = get_sum_bound(inf_Z,inf_weights,sup_weights,true);
				if (temp < sum_inf)	sum_inf = temp;
				temp = get_sum_bound(sup_Z,inf_weights,sup_weights,false);
				if (temp > sum_sup)	sum_sup = temp;
			
			
			}
			while (gsl_combination_next (comb) == GSL_SUCCESS);
		
			gsl_combination_free(comb);
		}
	}
	else	{
		sum_inf = get_sum_bound(inf_Z,inf_weights,sup_weights,true);
		sum_sup = get_sum_bound(sup_Z,inf_weights,sup_weights,false);
	}

	gsl_vector_set(result,0,sum_inf);
	gsl_vector_set(result,1,sum_sup);


	gsl_vector_free(inf_Z);
	gsl_vector_free(sup_Z);

}

