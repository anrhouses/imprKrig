/**
 * @file variogram.h
 * @brief implementation : variogram class
 * @author Tomohiko Mukai
 */

#include <mex.h>
#include <iostream>
//#include "utils.hpp"
#include "variogram.hpp"
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


using namespace std;

/**
 * @fn CVariogram::CVariogram(void)
 * @brief default constructor
 * @return none
 */
CVariogram::CVariogram(void)
{

	m_sill = 0.0;
	m_range = 0.0;
	m_nugget = 0.0;
	m_power = 0.0;
	
	m_pDistance = NULL;
	m_pVariogram = NULL;
	m_samples = 0;

	m_model = VARIO_NONE;
}

/**
 * @fn CVariogram::~CVariogram(void)
 * @brief destructor
 * @return none
 */
CVariogram::~CVariogram(void)
{
	destroy();
}

bool CVariogram::allocate(size_t smpl)
{
	destroy();      

	m_pDistance = new double[smpl];
	m_pVariogram = new double[smpl];

	if (m_pDistance == NULL || m_pVariogram == NULL)
	{
		destroy();
		return false;
	}
	m_samples = smpl;

	return true;
}


/**
 * @fn void CVariogram::destroy(void)
 * @brief memory destruction
 * @return none
 */
void CVariogram::destroy(void)
{
	if (m_pDistance != NULL)
	{
		delete[] m_pDistance;
		m_pDistance = NULL;
	}
	if (m_pVariogram != NULL)
	{
		delete[] m_pVariogram;
		m_pVariogram = NULL;
	}

	m_model = VARIO_NONE;
	m_nugget = 0.0;
	m_sill = 0.0;
	m_range = 0.0;
	m_power = 0.0;

}


bool CVariogram::setSample(const std::vector<CVariogram::VarioItm> &vecSample)
{
	if (vecSample.empty())
		return false;

	if (!allocate(vecSample.size()))
		return false;

	for (size_t i = 0;i < samples();i++)
	{
		m_pDistance[i] = vecSample[i].distance;
		m_pVariogram[i] = vecSample[i].dissimilarity;
	}

	sortByDistance();
	m_model = VARIO_NONE;
 
	return true;
}


bool CVariogram::getSample(size_t smpl, double &dist, double &vario) const
{
	if (smpl >= samples())
		return false;

	dist = m_pDistance[smpl];
	vario = m_pVariogram[smpl];
	return true;
}


/**
 * @fn bool CVariogram::setModel(int model, double nugget, double sill, double range, double power, double step, double maxdist)
 * @brief set of theoretical variogram model
 * @param model [in] theoretical variogram model
 * @param nugget [in] nugget
 * @param sill [in] sill
 * @param range [in] range
 * @param power [in] power coefficient of stable variograms
 * @param step [in] step width of reogionalization
 * @param maxdist [in] range of experimental variogram
 * @retval true success
 * @retval false fail
 */
bool CVariogram::setModel(int model, double nugget, double sill, double range, double power)
{
	if (model >= CVariogram::VARIO_NUM)
		return false;

	m_model = model;
	m_nugget = nugget;
	m_sill = sill;
	m_range = range;
	m_power = power;

	return true;
}



/**
 * @fn double CVariogram::getModelData(double dist) const
 * @brief get theoretical dissimilarity
 * @param dist [in] distance
 * @return dissimilarity
 */
double CVariogram::getModelData(double dist) const
{
	switch (m_model)
	{
	case VARIO_SPH:
		return spherical(m_nugget, m_sill, m_range, dist);
	case VARIO_STB:
		return stable(m_nugget, m_sill, m_range, m_power, dist);
	default:
		break;
	}
	return -1;
}

/**
 * @fn double CVariogram::getModelCovariance(double dist) const
 * @brief get covariance corresponding to theoretical dissimilarity
 * @param dist [in] distance
 * @return covariance
 */
double CVariogram::getModelCovariance(double dist) const
{
	switch (m_model)
	{
	case VARIO_SPH:
		return m_sill - spherical(m_nugget, m_sill, m_range, dist);
	case VARIO_STB:
		return m_sill - stable(m_nugget, m_sill, m_range, m_power, dist);
	default:
		break;
	}

	// if variogram does not has boundary ...
	return -1.0;
}

/*

bool CVariogram::draw(void)
{
	if (!isActive())
		return false;
	mglGraphPS gr;
	ofstream data_file;
	double max_distance=0;
	double max_vario=0;
	string input;
	size_t viewer_type=1;
	char *cmd  = (char*) malloc (1000);	
	char *plot_file_path = (char*) malloc (100);
	const char extension[100] = ".eps";
	
	do {
	
		strcpy(plot_file_path,"results/variogram");
		cout <<endl<< "Please enter the path of the file where to plot the cloud variogram without extension (default : " << plot_file_path << " ) ?" << endl << "The file extension will be " << extension << endl ;
		getline(cin, input);
		if (strcmp(input.c_str(),"")!=0)	
			strcpy(plot_file_path,const_cast<char*> ( input.c_str()));
		strcat (plot_file_path,extension);
		initOutFile(&data_file,plot_file_path);
	}
	while(!data_file);
	data_file.close();

	for(size_t k=0;k<m_samples;k++)
	{	
		if(m_pDistance[k]>max_distance) max_distance = m_pDistance[k];
		if(m_pVariogram[k]>max_vario) max_vario = m_pVariogram[k];
	}

	gr.SetRanges (0, 1.05*max_distance, 0, 1.05*max_vario);

	for(size_t k=0;k<m_samples;k++)
	{	
		gr.Mark(mglPoint(m_pDistance[k],m_pVariogram[k],0),'d');
	}
	gr.Light(true);
	gr.Box();
	gr.Axis("xy");
	gr.Grid();
	gr.WriteEPS(plot_file_path);    // Don't forget to save the result!
		
	while (true) {
		viewer_type = 1;
		cout << "Which viewer do you have ? "<<endl <<" 1. evince (document viewer) default"<<endl <<" 2. okular"<<endl;
		getline(cin, input);
		// This code converts from string to number safely.
		if (strcmp(input.c_str(),"")==0)
			break;
		stringstream myStream(input);
		if (myStream >> viewer_type)
			break;
		cout << "Invalid number, please try again" << endl;
	 }
	switch(viewer_type)
	{
		case 1:
		strcpy(cmd,"evince ");
		break;

		case 2:
		strcpy(cmd,"okular ");
		break;
	}
	strcat(cmd,plot_file_path);
	strcat(cmd," &");
	system(cmd);
	return true;
}

*/
bool CVariogram::sortByDistance(void)
{
	if (!isActive())
		return false;

	for (size_t i = 0;i < m_samples;i++)
	{
		for (size_t j = i + 1;j < m_samples;j++)
		{
			if (m_pDistance[i] > m_pDistance[j])
			{
				std::swap(m_pDistance[i], m_pDistance[j]);
				std::swap(m_pVariogram[i], m_pVariogram[j]);
			}
		}
	}
	return true;
}

