/**
 * @file geostat.cpp
 * @brief implementation : imprecise kriging predictor class
 * @author Kevin Loquin
 */
#include <mex.h>
#include "variogram.hpp"
#include "geostat.hpp"
#include <cmath>
#include <map>
#include <iostream>
#include <fstream>
#include <cstring>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <algorithm>
#include <ctime>
//#include "utils.hpp"



using namespace std;


bool CGeostat::allocate(size_t smpl, size_t dim)
{
	destroy();

	for (size_t s = 0;s < smpl;s++)
		m_coords.push_back(gsl_vector_alloc(dim));
		
	for (size_t s = 0;s < smpl;s++)
		m_data_bounds.push_back(gsl_vector_alloc(2));
	
	for (size_t s = 0;s < 3;s++)
		m_parameters_bounds.push_back(gsl_vector_alloc(2));
			
	for (size_t d = 0;d < dim;d++)
		m_domain.push_back( std::pair<double, double>(0, 0) );

	m_pData = gsl_vector_alloc(smpl);
	m_parameters = gsl_vector_alloc(3);
	m_pVariogram = new CVariogram;
	m_pVarioCloud = new CVariogram;
	
	m_pOKSys = gsl_matrix_alloc(samples() + 1, samples() + 1);

	return true;
}


void CGeostat::destroy(void)
{
	while (!m_coords.empty())
	{
		delete m_coords.back();
		m_coords.pop_back();
	}
	while (!m_data_bounds.empty())
	{
		delete m_data_bounds.back();
		m_data_bounds.pop_back();
	}
	while (!m_parameters_bounds.empty())
	{
		delete m_parameters_bounds.back();
		m_parameters_bounds.pop_back();
	}
	if (m_pData != NULL)
	{
		delete m_pData;
		m_pData = NULL;
	}
	if (m_parameters != NULL)
	{
	  delete m_parameters;
		m_parameters = NULL;
	}
	if (m_pVariogram != NULL)
	{
		delete m_pVariogram;
	        m_pVariogram = NULL;
	}
	if (m_pVarioCloud != NULL)
	{
		delete m_pVarioCloud;
		m_pVarioCloud = NULL;
	}
	if (m_pOKSys != NULL)
	{
		delete m_pOKSys;
		m_pOKSys = NULL;
	}
	m_domain.clear();
}




CGeostat::CGeostat(void)
{
	m_coords.clear();
	m_data_bounds.clear();
	m_parameters_bounds.clear();
	m_pData = NULL;
	m_parameters = NULL;
	m_domain.clear();
	m_pVariogram = NULL;
	m_pVarioCloud = NULL;
	m_pOKSys = NULL;
}




CGeostat::~CGeostat(void)
{
	destroy();
}

/*
* 
* This method scans the data file for learning
*   - number of data
*   - data nature (precise, interval, fuzzy)
*   - initialize the CGeostat object krig
*   - initialize the krig data
*
**/	


/*
bool CGeostat::initialize_data_file()
{

  destroy();
	
	ifstream data_file;
	size_t num_samples = 0;
	
	char str[2000];
	string input;
	
	data_file_path = (char*) malloc (100);
	
	do {
	
		strcpy(data_file_path,"data/data_irsn_imprecise");
		cout <<endl<< "Path of the data file (default : " << data_file_path << " ) ?" << endl;
		getline(cin, input);
		if (strcmp(input.c_str(),"")!=0)	strcpy(data_file_path,const_cast<char*> ( input.c_str()));
		initInFile(&data_file,data_file_path);
	}
	while(!data_file);
	
	data_nature = 0;
	while(!data_file.eof())
	{
		data_file.getline(str,2000);
		if (num_samples == 0)	{
			for (size_t k=0;k<strlen(str);k++)	{
				if(str[k]==' ') data_nature++;
			}
			data_nature = data_nature-1;
		}
		if(strcmp(str,"")!=0)
			num_samples++;
	}
	data_file.close();

	return allocate(num_samples, 2);
}
*/



/*
* 
* This method scans the data file for initializing the data
*
**/	

 /*

bool CGeostat::initialize_data()
{


	vector<CVariogram::VarioItm> vcloud;

	vcloud = computeVariogramCloud();
	m_pVarioCloud->setSample(vcloud);

	return (m_pData != NULL);

}

 */

/*
*
* define a domain of control space
*
*/


  /*
bool CGeostat::initialize_domain()	{

	
	string input;
	bool correct_answer = true;
	double inf_dom,sup_dom;
	gsl_vector *lower = gsl_vector_alloc(dimension());
	gsl_vector *upper = gsl_vector_alloc(dimension());
	
	gsl_vector *X = gsl_vector_alloc(samples());	
	gsl_vector *Y = gsl_vector_alloc(samples());	
	getCoordinates(X,Y);
	double X_min = gsl_vector_min(X);
	double X_max = gsl_vector_max(X);
	double Y_min = gsl_vector_min(Y);
	double Y_max = gsl_vector_max(Y);

	do {
		cout <<endl<< "Do you want an automatic kriging domain definition (y/n) ? (default: y)" << endl << "If automatic : the domain is an expansion of 10 p.c. in x and y directions of the smallest domain covering all the data" << endl;
		getline(cin, input);
		correct_answer = true;
		if (strcmp(input.c_str(),"y")==0 || strcmp(input.c_str(),"")==0)
		{
			
			inf_dom = X_min-5.0*(X_max-X_min)/100;
			sup_dom = X_max+5.0*(X_max-X_min)/100;
			gsl_vector_set(lower, 0, inf_dom);
			gsl_vector_set(upper, 0, sup_dom);

			inf_dom = Y_min-5.0*(Y_max-Y_min)/100;
			sup_dom = Y_max+5.0*(Y_max-Y_min)/100;
			gsl_vector_set(lower, 1, inf_dom);
			gsl_vector_set(upper, 1, sup_dom);
			
		}
		else if  (strcmp(input.c_str(),"n")==0)	
		{
			
			cout << "OK. First, the X-domain."<<endl;

			while(true) {
				cout << "Please enter the lower bound. It must be smaller than " << X_min << endl;
				getline(cin, input);
				stringstream myStream(input);
				if (myStream >> inf_dom)
				{
					if (inf_dom > X_min)
						cout << "!!!! It must be smaller than " << X_min << " !!!!" << endl;
					else break;
				}
				cout << "Invalid number, please try again" << endl;
			} 
			
			gsl_vector_set(lower, 0, inf_dom);
			
			while(true) {
				cout << "OK. Now, please enter the upper bound. It must be greater than " << X_max << endl;
				getline(cin, input);
				stringstream myStream(input);
				if (myStream >> sup_dom)
				{
					if (sup_dom < X_max)
						cout << "!!!! It must be greater than " << X_max << " !!!!" << endl;
					else break;
				}
				cout << "Invalid number, please try again" << endl;
			} 
			
			gsl_vector_set(upper, 0, sup_dom);

			cout << "OK. Second, the Y-domain."<<endl;

			while(true) {
				cout << "Please enter the lower bound. It must be smaller than " << Y_min << endl;
				getline(cin, input);
				stringstream myStream(input);
				if (myStream >> inf_dom)
				{
					if (inf_dom > Y_min)
						cout << "!!!! It must be smaller than " << Y_min << " !!!!" << endl;
					else break;
				}
				cout << "Invalid number, please try again" << endl;
			} 
			
			
			gsl_vector_set(lower, 1, inf_dom);

			while(true) {
				cout << "OK. Now, please enter the upper bound. It must be greater than " << Y_max << endl;
				getline(cin, input);
				stringstream myStream(input);
				if (myStream >> sup_dom)
				{
					if (sup_dom < Y_max)
						cout << "!!!! It must be greater than " << Y_max << " !!!!" << endl;
					else break;
				}
				cout << "Invalid number, please try again" << endl;
			} 
			gsl_vector_set(upper, 1, sup_dom);
		}
		else {correct_answer = false;}
	}
	while(!correct_answer);
	
	//
	// set domain
	correct_answer = setDomain(lower, upper);

	gsl_vector_free(upper);
	gsl_vector_free(lower);
	gsl_vector_free(X);
	gsl_vector_free(Y);

	return correct_answer;

}
  */

  /*
bool CGeostat::initialize_variogram()	{

	size_t parameter_type;
	string input = "";
	double current_param;
	gsl_vector * current_imp_param = gsl_vector_alloc(2);

	cout <<endl<< "Do you want to plot the empirical variograms in eps files (y/n) ? (default:n) " << endl;
	getline(cin, input);

	if (strcmp(input.c_str(),"y")==0)
	{
		m_pVarioCloud->draw();
	}

	cout << "Variogram specification... "<< endl ;
	while (true) {
		parameter_type = 2;
		cout << "Which kind of variogram parameters ? "<<endl <<" 1. Precise"<<endl <<" 2. Interval (default)"<<endl <<" 3. Fuzzy " << endl;
		getline(cin, input);
		// This code converts from string to number safely.
		if (strcmp(input.c_str(),"")==0)
			break;
		stringstream myStream(input);
		if (myStream >> parameter_type)
			break;
		cout << "Invalid number, please try again" << endl;
	}
	vario_nature = parameter_type ;
	switch(parameter_type)
	{
		case 1:
		while(true)	
		{
			current_param = 0.01822;
			cout << "Enter the nugget effect value (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		setVariogramParameter(0 , current_param);
		while(true)	
		{
			current_param = 0.0682;
			cout << "Enter the sill value (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		setVariogramParameter(1 , current_param);
		while(true)	{
			current_param = 11.3;
			cout << "Enter the range value (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		setVariogramParameter(2 , current_param);
				
		break;

		case 2:
		
		cout << "Imprecise nugget effect ?" << endl;
		while(true)	
		{
			current_param = 0.0;	
			
			cout << "Lower bound (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		gsl_vector_set(current_imp_param,0,current_param);
		while(true)	
		{
			current_param = 0.04;	
			
			cout << "Upper bound (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		gsl_vector_set(current_imp_param,1,current_param);
		setImpreciseVariogramParameter(0 , current_imp_param);		
		setVariogramParameter(0 , (current_param + gsl_vector_get(current_imp_param,0) )/2.0);


		cout << "Imprecise sill ?" << endl;
		while(true)	
		{
			current_param = 0.06;	
			
			cout << "Lower bound (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		gsl_vector_set(current_imp_param,0,current_param);
		while(true)	
		{
			current_param = 0.11;	
			
			cout << "Upper bound (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		gsl_vector_set(current_imp_param,1,current_param);
		setImpreciseVariogramParameter(1 , current_imp_param);		
		setVariogramParameter(1 , (current_param + gsl_vector_get(current_imp_param,0) )/2.0);
				

		cout << "Imprecise range ?" << endl;
		while(true)	
		{
			current_param = 10.0;	
			
			cout << "Lower bound (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		gsl_vector_set(current_imp_param,0,current_param);
		while(true)	
		{
			current_param = 12.0;	
			
			cout << "Upper bound (default = "<< current_param << ") ?" << endl;
			getline(cin, input);

			// This code converts from string to number safely.
			if (strcmp(input.c_str(),"")==0)
				break;

			stringstream myStream(input);
			if (myStream >> current_param)
				break;
			cout << "Invalid number, please try again" << endl;
	 	}
		gsl_vector_set(current_imp_param,1,current_param);
		setImpreciseVariogramParameter(2 , current_imp_param);		
		setVariogramParameter(2 , (current_param + gsl_vector_get(current_imp_param,0) )/2.0);

		break;

		case 3:
		cout << "fuzzy nugget effect ?" << endl;

		break;
	}
	gsl_vector_free(current_imp_param);
	return true;
}
  */

  /*
bool CGeostat::initialize_simulations_number()
{
  size_t current_param  ;
  string input = "";

  while(true)	
    {
      current_param = 10 ;	
    
      cout << "number of simulations of the simulated annealing method (default = "<< current_param << ") ?" << endl;
      getline(cin, input);
      
      // This code converts from string to number safely.
      if (strcmp(input.c_str(),"")==0)
	break;
      stringstream myStream(input);
      if (myStream >> current_param)
	break;
      cout << "Invalid number, please try again" << endl;
    }
  n_simulations = current_param ;
  return true;
}
*/

bool CGeostat::getCoordinate(gsl_vector *c, size_t s) const
{
	if (!isActive() || s >= samples())
		return false;
	gsl_vector_memcpy(c, m_coords[s]);
	
	return true;
}

bool CGeostat::getCoordinates(gsl_vector *X, gsl_vector *Y) const
{

	gsl_vector *c = gsl_vector_alloc(dimension());	
	for (size_t s=0;s<samples(); s++)
	{
		getCoordinate(c, s);
		gsl_vector_set(X,s,gsl_vector_get(c,0));
		gsl_vector_set(Y,s,gsl_vector_get(c,1));
	}
	gsl_vector_free(c);
	return ( (X!=NULL) && (Y!=NULL) );
}


double CGeostat::getData(size_t s) const
{
	if (!isActive() || s >= samples())
		return 0.0;

	return gsl_vector_get(m_pData, s);
}


bool CGeostat::getImpreciseData(gsl_vector *c, size_t s) const
{
	if (!isActive() || s >= samples())
		return false;
		
	gsl_vector_memcpy(c, m_data_bounds[s]);
	
	return true;
}

size_t CGeostat::get_data_nature() const
{
	return data_nature;
}

size_t CGeostat::get_vario_nature() const
{
	return vario_nature;
}


double CGeostat::getVariogramParameter(size_t s) const
{
	if (!isActive() || s >= samples())
		return false;
		
	return gsl_vector_get(m_parameters,s);
}

bool CGeostat::getImpreciseVariogramParameter(gsl_vector *c, size_t s) const
{
	if (!isActive() || s >= samples())
		return false;
		
	gsl_vector_memcpy(c, m_parameters_bounds[s]);
	
	return true;
}



bool CGeostat::getDomain(gsl_vector *lower, gsl_vector *upper) const
{
	if (!isActive())
		return false;

	for (size_t d = 0;d < dimension();d++)
	{
		gsl_vector_set(lower, d, m_domain[d].first);
		gsl_vector_set(upper, d, m_domain[d].second);
	}

	return true;
}


bool CGeostat::setDomain(const gsl_vector *lower, const gsl_vector *upper)
{
	if (!isActive())
		return false;

	for (size_t d = 0;d < dimension();d++)
	{
		m_domain[d].first = gsl_vector_get(lower, d);
		m_domain[d].second = gsl_vector_get(upper, d);
	}

	return true;
}


bool CGeostat::setDomain()
{
	if (!isActive())
		return false;
	gsl_vector *lower = gsl_vector_alloc(dimension());	
	gsl_vector *upper = gsl_vector_alloc(dimension());

	gsl_vector *X = gsl_vector_alloc(samples());	
	gsl_vector *Y = gsl_vector_alloc(samples());	

	getCoordinates(X,Y);
	double X_min = gsl_vector_min(X);
	double X_max = gsl_vector_max(X);
	double Y_min = gsl_vector_min(Y);
	double Y_max = gsl_vector_max(Y);
	
	gsl_vector_set(lower, 0, X_min-5.0*(X_max-X_min)/100);
	gsl_vector_set(upper, 0, X_max+5.0*(X_max-X_min)/100);
	
	gsl_vector_set(lower, 1, Y_min-5.0*(Y_max-Y_min)/100);
	gsl_vector_set(upper, 1, Y_max+5.0*(Y_max-Y_min)/100);
	
	bool correct_answer = setDomain(lower, upper);

	gsl_vector_free(upper);
	gsl_vector_free(lower);
	gsl_vector_free(X);
	gsl_vector_free(Y);

	return correct_answer;

}


bool CGeostat::setSimulation(size_t n_sims)
{
        if (!isActive())
		return false;

	n_simulations = n_sims;

	return true;
}


bool CGeostat::setDataNature(size_t data_nat)
{
        if (!isActive())
		return false;

	data_nature = data_nat;

	return true;
}


bool CGeostat::setCoordinate(size_t s, const gsl_vector *c)
{

	if (!isActive() || s >= samples())
		return false;

	gsl_vector_memcpy(m_coords[s], c);
	return true;
}


bool CGeostat::setImpreciseData(size_t s, const gsl_vector *c)
{

	if (!isActive() || s >= samples())
		return false;

	gsl_vector_memcpy(m_data_bounds[s], c);
	return true;
}


bool CGeostat::setImpreciseVariogramParameter(size_t s, const gsl_vector *c)
{

	if (!isActive() || s >= samples())
		return false;

	gsl_vector_memcpy(m_parameters_bounds[s], c);
	return true;
}


void CGeostat::display()
{
	if (!isActive())
		exit(1);
	for(size_t s=0;s<samples();s++)	{
	  mexPrintf( "Position : ( %g , %g). Precise data : %g. Imprecise data : ( %g , %g ) \n",gsl_vector_get(m_coords[s], 0), gsl_vector_get(m_coords[s], 1), getData(s),gsl_vector_get(m_data_bounds[s], 0),gsl_vector_get(m_data_bounds[s], 1));
	  // << "   Imprecise data :" << gsl_vector_get(m_data_bounds[s], 0) << " , " << gsl_vector_get(m_data_bounds[s], 1) << endl ;
	  //cout << "Position : ( " << gsl_vector_get(m_coords[s], 0) << " , " << gsl_vector_get(m_coords[s], 1) << " ). Precise data : " << getData(s) << "   Imprecise data :" << gsl_vector_get(m_data_bounds[s], 0) << " , " << gsl_vector_get(m_data_bounds[s], 1) << endl ;
	}
}


bool CGeostat::setData(size_t s, double d)
{	
	if (!isActive() || s >= samples())
		return false;

	gsl_vector_set(m_pData, s, d);

	
	return true;
}

bool CGeostat::setVariogramParameter(size_t s, double d)
{
	
	if (!isActive() || s >= 3)
		return false;
	gsl_vector_set(m_parameters, s, d);

	return true;
}	


bool CGeostat::getModelParameters(double &sill, double &range, double &nugget, double &power) const
{
	if (!isActive())
		return false;

	sill = m_pVariogram->sill();
	range = m_pVariogram->range();
	nugget = m_pVariogram->nugget();
	power = m_pVariogram->power();

	return false;
}


vector<CVariogram::VarioItm> CGeostat::computeVariogramCloud(void) const
{
	vector<CVariogram::VarioItm> vecVario;
	CVariogram::VarioItm itm;
	gsl_vector *c1, *c2;
	double r1, r2;

	if (!isActive())
		return vecVario;

	c1 = gsl_vector_alloc(dimension());
	c2 = gsl_vector_alloc(dimension());
		
	for (size_t i = 0;i < samples() - 1;i++)
	{
		getCoordinate(c1, i);
		r1 = getData(i);

		for (size_t j = i + 1;j < samples();j++)
		{
			getCoordinate(c2, j);
			r2 = getData(j);

			itm.distance = isoDist(c1, c2);
			itm.dissimilarity = (r1 - r2) * (r1 - r2) / 2.0;
		
			vecVario.push_back(itm);
		}
	}

	gsl_vector_free(c1);
	gsl_vector_free(c2);
	m_pVarioCloud->setSample(vecVario);
	return vecVario;
}


vector<CVariogram::VarioItm> CGeostat::computeExperimentalVariogram(const std::vector<CVariogram::VarioItm> &vcloud, double step) const
{
	vector<CVariogram::VarioItm> vmodel;
	map<int, CVariogram::VarioItm> variomap;
	CVariogram::VarioItm itm;

	for (size_t i = 0;i < vcloud.size();i++)
	{
		int key = static_cast<int>(vcloud[i].distance / step);
		if (variomap.find(key) == variomap.end())
		{
			variomap[key].distance = 1.0;
			variomap[key].dissimilarity = vcloud[i].dissimilarity;
		}
		else
		{
			variomap[key].distance += 1.0;
			variomap[key].dissimilarity += vcloud[i].dissimilarity;
		}
	}

	for (std::map<int, CVariogram::VarioItm>::const_iterator it = variomap.begin();it != variomap.end();++it)
	{
		itm.distance = it->first * step;
		itm.dissimilarity = it->second.dissimilarity / it->second.distance;
		vmodel.push_back(itm);
	}
	return vmodel;
}






bool CGeostat::precomputeKrigingSystemFromOutSide()
{
  return precomputeKrigingSystem(getVariogramParameter(0), getVariogramParameter(1), getVariogramParameter(2));
}


bool CGeostat::precomputeKrigingSystem(double nugget, double sill, double range)
{
	gsl_vector *c1, *c2;
	gsl_matrix *tmp;
	int signum;
	gsl_permutation *p;

	c1 = gsl_vector_alloc(dimension());
	c2 = gsl_vector_alloc(dimension());
	tmp = gsl_matrix_alloc(samples() + 1, samples() + 1);
	p = gsl_permutation_alloc(samples() + 1);
	
	m_pVariogram->setModel(CVariogram::VARIO_SPH, nugget, sill, range, 0);
	
	for (size_t s1 = 0;s1 < samples();s1++)
	{
		getCoordinate(c1, s1);
		
		for (size_t s2 = s1;s2 < samples();s2++)
		{
			getCoordinate(c2, s2);

			double dist;
			dist = isoDist(c1, c2);

			gsl_matrix_set(tmp, s1, s2, m_pVariogram->getModelData(dist));
			gsl_matrix_set(tmp, s2, s1, gsl_matrix_get(tmp, s1, s2));
		}
	}

	for (size_t s = 0;s < samples();s++)
	{
		gsl_matrix_set(tmp, s, samples(), 1.0);
		gsl_matrix_set(tmp, samples(), s, 1.0);
	}
	gsl_matrix_set(tmp, samples(), samples(), 0.0);
	
	gsl_linalg_LU_decomp(tmp, p, &signum);
	gsl_linalg_LU_invert(tmp, p, m_pOKSys);

	gsl_permutation_free(p);

	gsl_vector_free(c1);
	gsl_vector_free(c2);
	gsl_matrix_free(tmp);

	return true;
}


double CGeostat::maxDist() const
{
	double dist =0.0;

	for (size_t d = 0;d < dimension();d++)
		dist += pow(m_domain[d].second - m_domain[d].first, 2.0);
	return sqrt(dist);
}


bool CGeostat::isDomain(const gsl_vector *c) const
{
	for (size_t d = 0;d < dimension();d++)
	{
		double ct = gsl_vector_get(c, d);
		if (ct < m_domain[d].first || ct > m_domain[d].second)
			return false;
	}
	return true;
}


bool CGeostat::getWeightVector(gsl_vector *weight, const gsl_vector *c) const
{
	if (!isActive())
		return false;

	gsl_vector *spos, *vvario;
	
	spos = gsl_vector_alloc(dimension());
	vvario = gsl_vector_alloc(samples() + 1);

	for (size_t s = 0;s < samples();s++)
	{
		getCoordinate(spos, s);
		double dist = isoDist(c, spos);

		gsl_vector_set(vvario, s, m_pVariogram->getModelData(dist));
	}
	gsl_vector_set(vvario, samples(), 1.0);

	gsl_blas_dgemv(CblasNoTrans, 1.0, m_pOKSys, vvario, 0.0, weight);


	gsl_vector_free(spos);
	gsl_vector_free(vvario);

	return true;
}

double CGeostat::getWeight(size_t i, const gsl_vector *c) const
{
	if (!isActive())
		return false;
	
	double res = 0;
	gsl_vector *spos, *vvario, *row;
	
	spos = gsl_vector_alloc(dimension());
	vvario = gsl_vector_alloc(samples() + 1);
	row = gsl_vector_alloc(samples() + 1);

	for (size_t s = 0;s < samples();s++)
	{
		getCoordinate(spos, s);
		double dist = isoDist(c, spos);

		gsl_vector_set(vvario, s, m_pVariogram->getModelData(dist));
	}
	gsl_vector_set(vvario, samples(), 1.0);

	//gsl_blas_dgemv(CblasNoTrans, 1.0, m_pOKSys, vvario, 0.0, weight);
	gsl_matrix_get_row (row, m_pOKSys, i);
	for(size_t j=0;j<=samples();j++)
		res+=gsl_vector_get(vvario,j)*gsl_vector_get(row,j);

	gsl_vector_free(spos);
	gsl_vector_free(vvario);

	return res;
}



bool CGeostat::getPredictData(double &pred, double &var, const gsl_vector *c) const
{
	if (!isActive())
		return false;

	
	double condition =0;
		
		
	gsl_vector *spos;
	gsl_vector *vvario, *weight;
	spos = gsl_vector_alloc(dimension());
	vvario = gsl_vector_alloc(samples() + 1);
	weight = gsl_vector_alloc(samples() + 1);
	
	for (size_t s = 0;s < samples();s++)
	{
		getCoordinate(spos, s);
		double dist = isoDist(c, spos);

		gsl_vector_set(vvario, s, m_pVariogram->getModelData(dist));
	}
	gsl_vector_set(vvario, samples(), 1.0);

	gsl_blas_dgemv(CblasNoTrans, 1.0, m_pOKSys, vvario, 0.0, weight);

	for (size_t s = 0;s < samples();s++)
		condition += gsl_vector_get(weight, s) ;
		
	//pred = gsl_vector_get(weight, samples());		
	
	pred = 0;
	for (size_t s = 0;s < samples();s++)
		pred += gsl_vector_get(weight, s) * getData(s);

	gsl_blas_ddot(weight, vvario, &var);

	gsl_vector_free(spos);
	gsl_vector_free(vvario);
	gsl_vector_free(weight);

	return true;
}

bool CGeostat::getOptimalImprecisePrediction(gsl_vector *pred, const gsl_vector *c) const
{
	if (!isActive())
		return false;
		
	gsl_vector *spos;
	gsl_vector *vvario, *weight;
	gsl_vector *dat;
	double pred_inf = 0,pred_sup = 0;
	double tmp_weight,dist;
	
	spos = gsl_vector_alloc(dimension());
	vvario = gsl_vector_alloc(samples() + 1);
	weight = gsl_vector_alloc(samples() + 1);
	dat = gsl_vector_alloc(2);
	for (size_t s = 0;s < samples();s++)
	{
		getCoordinate(spos, s);
		dist = isoDist(c, spos);

		gsl_vector_set(vvario, s, m_pVariogram->getModelData(dist));
	}
	gsl_vector_set(vvario, samples(), 1.0);

	gsl_blas_dgemv(CblasNoTrans, 1.0, m_pOKSys, vvario, 0.0, weight);


	for (size_t s = 0;s < samples();s++)	{
		getImpreciseData(dat,s);
		tmp_weight = gsl_vector_get(weight, s) ;
		if(tmp_weight<0)
	        {
			pred_inf +=  tmp_weight * gsl_vector_get(dat, 1);
			pred_sup +=  tmp_weight * gsl_vector_get(dat, 0);
		}
		else
		{
			pred_inf +=  tmp_weight * gsl_vector_get(dat, 0);
			pred_sup +=  tmp_weight * gsl_vector_get(dat, 1);
		}
	}
	gsl_vector_set(pred,0,pred_inf);
	gsl_vector_set(pred,1,pred_sup);
	gsl_vector_free(spos);
	gsl_vector_free(vvario);
	gsl_vector_free(weight);

	return true;
}




bool CGeostat::getPredictData(double &pred, const gsl_vector *weight) const
{
	if (!isActive())
		return false;

	pred = 0;
	for (size_t s = 0;s < samples();s++)
		pred += gsl_vector_get(weight, s) * getData(s);

	return true;
}


double isoDist(const gsl_vector *c1, const gsl_vector *c2)
{
	double norm;
	gsl_vector *tmp;

	tmp = gsl_vector_alloc(c1->size);
	gsl_vector_memcpy(tmp, c1);
	gsl_vector_sub(tmp, c2);
	norm = gsl_blas_dnrm2(tmp);
	gsl_vector_free(tmp);

	return norm;
}



bool CGeostat::isOnADataPosition(const gsl_vector *c, size_t *data_index)	{
	bool test = false;
	double dist = 0;
	gsl_vector *spos = gsl_vector_alloc(dimension());
	for (size_t s = 0;s < samples();s++)
	{
		getCoordinate(spos, s);
		dist = isoDist(c, spos);
		if(dist == 0)	{
			test = true;
			*data_index = s;
			break;
		}
	}
	return test;
}


void CGeostat::getOptimumBound_simulated_annealing(gsl_vector * pred, const gsl_vector *c)
{

	size_t data_index;
	
	if(isOnADataPosition(c, &data_index))	{
		getImpreciseData(pred,data_index);
	}
	else	{

		size_t k=0;
		double Tmax = 0.5 ;
		size_t nbCoeffs = 3;
		size_t security;

		double E_inf_Temp,E_inf,E_inf_New,E_sup_Temp,E_sup,E_sup_New,chosensep;
		
		gsl_vector * params_inf  = gsl_vector_alloc(nbCoeffs);
		gsl_vector * params_inf_Temp  = gsl_vector_alloc(nbCoeffs);
		gsl_vector * params_inf_New  = gsl_vector_alloc(nbCoeffs);

		gsl_vector * params_sup  = gsl_vector_alloc(nbCoeffs);
		gsl_vector * params_sup_Temp  = gsl_vector_alloc(nbCoeffs);
		gsl_vector * params_sup_New  = gsl_vector_alloc(nbCoeffs);
	
		int chosenSign;

		gsl_vector * temp_pred = gsl_vector_alloc(2);
	

		double tmp, test;
		double temperature , prob;
	
	
	
		gsl_vector * paramsInf  = gsl_vector_alloc(nbCoeffs);
		gsl_vector * paramsSup  = gsl_vector_alloc(nbCoeffs);
		gsl_vector * lengths  = gsl_vector_alloc(nbCoeffs);
	
	
		for(size_t i=0;i<nbCoeffs;i++)	{
			gsl_vector_set(paramsInf, i, gsl_vector_get(m_parameters_bounds[i],0)); gsl_vector_set(paramsSup, i, gsl_vector_get(m_parameters_bounds[i],1)); 

			gsl_vector_set(lengths, i, fabs(gsl_vector_get(paramsSup,i)-gsl_vector_get(paramsInf,i)));
		
			tmp = (gsl_vector_get(paramsInf,i) + gsl_vector_get(paramsSup,i))/2.0;
			
			
			gsl_vector_set(params_inf, i, tmp);
			gsl_vector_set(params_inf_Temp, i, tmp);
			gsl_vector_set(params_inf_New, i, tmp);
			gsl_vector_set(params_sup, i, tmp);
			gsl_vector_set(params_sup_Temp, i, tmp);
			gsl_vector_set(params_sup_New, i, tmp); 
		}


		precomputeKrigingSystem(gsl_vector_get(params_inf,0), gsl_vector_get(params_inf,1), gsl_vector_get(params_inf,2));
		getOptimalImprecisePrediction(temp_pred, c);

		E_inf_Temp =  E_inf = gsl_vector_get(temp_pred,0);
		E_sup_Temp =  E_sup = - gsl_vector_get(temp_pred,1);

		srand ( time(NULL) );
		k=1;
		  while(k <= simulations())  {

		
			for (size_t i=0; i< nbCoeffs;i++)	{
				security = 0;
				do {
					security++;
					chosenSign =(int)rand()%2;
					if(chosenSign==0) chosenSign=-1;
					chosensep = chosenSign * gsl_vector_get(lengths,i) * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
					tmp = gsl_vector_get(params_inf_Temp,i)+chosensep;
				}
				while(( tmp < gsl_vector_get(paramsInf,i) || tmp > gsl_vector_get(paramsSup,i) ) && security < 10000);
				gsl_vector_set(params_inf_New, i, tmp);

				do {
					security++;
					chosenSign =(int)rand()%2;
					if(chosenSign==0) chosenSign=-1;
					chosensep = chosenSign * gsl_vector_get(lengths,i) * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
					tmp = gsl_vector_get(params_sup_Temp,i)+chosensep;
				}
				while(( tmp < gsl_vector_get(paramsInf,i) || tmp > gsl_vector_get(paramsSup,i) ) && security < 10000);
				gsl_vector_set(params_sup_New, i, tmp);
			}



			precomputeKrigingSystem(gsl_vector_get(params_inf_New,0), gsl_vector_get(params_inf_New,1), gsl_vector_get(params_inf_New,2));
			getOptimalImprecisePrediction(temp_pred, c);

			E_inf_New = gsl_vector_get(temp_pred,0);
			E_sup_New = -gsl_vector_get(temp_pred,1);

			  temperature = Tmax*(1.0-static_cast<double>(k) / static_cast<double>(simulations()));

			test = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
			prob = probaAcceptance(E_inf_Temp,E_inf_New,temperature);
			if(prob>test)  {
				gsl_vector_memcpy(params_inf_Temp,params_inf_New);
				E_inf_Temp = E_inf_New;
			}
			if (E_inf_New<E_inf)  {
				gsl_vector_memcpy(params_inf,params_inf_New);
				E_inf = E_inf_New;
			}

			test = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
			prob = probaAcceptance(E_sup_Temp,E_sup_New,temperature);
			if(prob>test)  {
				gsl_vector_memcpy(params_sup_Temp,params_sup_New);
				E_sup_Temp = E_sup_New;
			}
			if (E_sup_New<E_sup)  {
				gsl_vector_memcpy(params_sup,params_sup_New);
				E_sup = E_sup_New;
			}

			precomputeKrigingSystem(gsl_vector_get(params_sup_New,0), gsl_vector_get(params_sup_New,1), gsl_vector_get(params_sup_New,2));
			getOptimalImprecisePrediction(temp_pred, c);

			E_inf_New = gsl_vector_get(temp_pred,0);
			E_sup_New = -gsl_vector_get(temp_pred,1);

			  temperature = Tmax*(1.0-static_cast<double>(k) / static_cast<double>(simulations()));	       

			test = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
			prob = probaAcceptance(E_inf_Temp,E_inf_New,temperature);
			if(prob>test)  {
				gsl_vector_memcpy(params_inf_Temp,params_inf_New);
				E_inf_Temp = E_inf_New;
			}
			if (E_inf_New<E_inf)  {
				gsl_vector_memcpy(params_inf,params_inf_New);
				E_inf = E_inf_New;
			}

			test = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
			prob = probaAcceptance(E_sup_Temp,E_sup_New,temperature);
			if(prob>test)  {
				gsl_vector_memcpy(params_sup_Temp,params_sup_New);
				E_sup_Temp = E_sup_New;
			}
			if (E_sup_New<E_sup)  {
				gsl_vector_memcpy(params_sup,params_sup_New);
				E_sup = E_sup_New;
			}


			k++;
		}	
		gsl_vector_set(pred,0,E_inf);
		gsl_vector_set(pred,1,-E_sup);


		gsl_vector_free(paramsInf);
		gsl_vector_free(paramsSup);
		gsl_vector_free(params_inf);
		gsl_vector_free(params_inf_Temp);
		gsl_vector_free(params_inf_New);
		gsl_vector_free(params_sup);
		gsl_vector_free(params_sup_Temp);
		gsl_vector_free(params_sup_New);
		gsl_vector_free(lengths);
		gsl_vector_free(temp_pred);
	}
	
}

void CGeostat::getOptimumBound_bounds_combinatorial(gsl_vector * pred, const gsl_vector *c)
{
	
	size_t data_index;
	
	if(isOnADataPosition(c, &data_index))	{
		//gsl_vector *imprecise_data = gsl_vector_alloc(2);
		getImpreciseData(pred,data_index);
		/*
		if(isInf) {predictor = gsl_vector_get(imprecise_data,0);}
		else {predictor = gsl_vector_get(imprecise_data,1);}
		gsl_vector_free(imprecise_data);
		*/
	}
	else	{
	
	size_t nbCoeffs = 3;
	gsl_vector * params  = gsl_vector_alloc(nbCoeffs);
	

	double pred_inf,pred_sup;
	gsl_vector * temp_pred = gsl_vector_alloc(2);
	
	
	gsl_vector * domain_length  = gsl_vector_alloc(nbCoeffs);
	
	for(size_t i=0;i<nbCoeffs;i++)	{
		gsl_vector_set(params, i, gsl_vector_get(m_parameters_bounds[i],0));
		gsl_vector_set(domain_length, i, fabs(gsl_vector_get(m_parameters_bounds[i],1)-gsl_vector_get(m_parameters_bounds[i],0)));
	}

	precomputeKrigingSystem(gsl_vector_get(params,0), gsl_vector_get(params,1), gsl_vector_get(params,2));
	getOptimalImprecisePrediction(temp_pred, c);
	pred_inf = gsl_vector_get(temp_pred,0);
	pred_sup = gsl_vector_get(temp_pred,1);
	for (size_t nug_i=0;nug_i<2;nug_i++)	{
		for (size_t sill_i=0;sill_i<2;sill_i++)	{
			for (size_t range_i=0;range_i<2;range_i++)	{
				if(!(nug_i==0 && sill_i==0 && range_i==0))	{
					precomputeKrigingSystem(gsl_vector_get(params,0)+nug_i*gsl_vector_get(domain_length,0), gsl_vector_get(params,1)+sill_i*gsl_vector_get(domain_length,1), gsl_vector_get(params,2)+range_i*gsl_vector_get(domain_length,2));
					getOptimalImprecisePrediction(temp_pred, c);
					if(gsl_vector_get(temp_pred,0)<pred_inf) {
						pred_inf = gsl_vector_get(temp_pred,0);
					}
					if(gsl_vector_get(temp_pred,1)>pred_sup)	{
						pred_sup = gsl_vector_get(temp_pred,1);
					}
				}
			}
		}
	}
	gsl_vector_set(pred,0,pred_inf);
	gsl_vector_set(pred,1,pred_sup);	

	gsl_vector_free(temp_pred);	
	gsl_vector_free(params);
	gsl_vector_free(domain_length);
		
	}

}


void CGeostat::getOptimumBound_hybrid(gsl_vector * pred, const gsl_vector *c)
{
	
	size_t data_index;
	
	if(isOnADataPosition(c, &data_index))	{
		getImpreciseData(pred,data_index);
	}
	else	{

	double pred_inf,pred_sup;
	gsl_vector * temp_pred = gsl_vector_alloc(2);
	
	getOptimumBound_bounds_combinatorial(temp_pred,c);
	pred_inf = gsl_vector_get(temp_pred,0);
	pred_sup = gsl_vector_get(temp_pred,1);
	getOptimumBound_simulated_annealing(temp_pred,c);		
	if(pred_inf > gsl_vector_get(temp_pred,0))
		pred_inf = gsl_vector_get(temp_pred,0);
	if(pred_sup < gsl_vector_get(temp_pred,1))
		pred_sup = gsl_vector_get(temp_pred,1);
	gsl_vector_set(pred,0,pred_inf);
	gsl_vector_set(pred,1,pred_sup);

	gsl_vector_free(temp_pred);	

	}

}

