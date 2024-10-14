#include <mex.h>

#include <iostream>
#include <sstream>

#include <iomanip>

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <cstring>
#include <string>


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


//#include "utils.hpp"
#include "geostat.hpp"

//#include "fuzzyinterval.hpp"
//#include "fuzzyintervalbycuts.hpp"



#define MESSAGE_D_ERREUR "intervallist_kriging fonction call is : [ INF_PRED, SUP_PRED ] = precise_kriging(DATA, VARIO_PARAMS, DOMAIN, SAMPLING, METHOD, N_SIMULATIONS)\n\nINPUT PARAMETERS\n\n\tDATA (matrix nx4) is the data set: column 1 are the X locations ; column 2 are the Y locations ; column 3 are the inferior data values ; column 4 are the superior data values\n\n\tVARIO_PARAMS (vector 2x3 or 3x2) made of the intervallist nugget effect, sill and range of the variogram\n\n\tDOMAIN can be\n \t\t(i) any integer, it results in an automatic domain definition\n \t\t(ii) a matrix 2x2 - row 1: X and Y positions of one angular extreme point of the squared domain - row 2: X and Y positions of its opposite extreme point\n\n\tSAMPLING can be\n\t\t(i) automatic sampling: SAMPLING is an integer, which is the number of sampling positions on a grid automatically defined on the DOMAIN\n\t\t(ii) point sampling: SAMPLING is a matrix px2 - each row is made of the X and Y positions of p sampling points\n\t\t(iii) grid sampling: SAMPLING is a matrix 2xp - row 1: is made of a p-sampling of the X line - row 2: is made of a p-sampling of the Y line\n\n\tMETHOD is an integer parameter which corresponds ot\n\t\t1 - hybrid method\n\t\t2 - simulated annealing method only\n\t\t3 - combinatorial method only\n\n\tN_SIMULATIONS is an integer parameter which corresponds to the number of simulations that will be used in the simulated annealing method\n\nOUTPUT PARAMETERS (their nature px1 or pxp depends on the SAMPLING parameter\n\n\tINF_PRED is the inferior kriging result\n\n\tSUP_PRED is the superior kriging result\n\n"

#define POINT 1
#define LINE 2
#define GRID 3

#define HYBRID 1
#define SIMULATED 2
#define COMBINATORIAL 3

#define Min(a,b) ( ( (a) < (b) ) ? (a) : (b) )
#define Max(a,b) ( ( (a) > (b) ) ? (a) : (b) )



/*
  #define NOYAU_TRIWEIGHT 6 
  #define NOYAU_COSINUS 7
  #define NOYAU_TRIANGLE_SUR_QUADRATIQUE 8
  #define NOYAU_GAUSSIEN 9
*/


void mexFunction(int NbOut, mxArray *PtOut[],int NbIn, const mxArray *PtIn[])
{

  CGeostat krig ;

  gsl_vector * current_position = gsl_vector_alloc(2);
  gsl_vector * current_data = gsl_vector_alloc(2);
  gsl_vector * current_params = gsl_vector_alloc(2);
  gsl_vector * current_prediction = gsl_vector_alloc(2);
  
  gsl_matrix * positions;
  gsl_matrix * inf_predictions;
  gsl_matrix * sup_predictions; 
  
  gsl_vector *lower = gsl_vector_alloc(2);
  gsl_vector *upper = gsl_vector_alloc(2);

  //gsl_matrix * data_positions;
  //gsl_vector * data;
  
  size_t estimation_nature, n_estimation_positions, n_data, k, i, j, method;
  bool automatic_domain_definition = false;
  bool automatic_grid_definition = false;
  //  size_t step = 0;

  double *pt,  *fin ;
  double Xm,Ym,Xstep,Ystep;
  double pred,var;

  /*
    mexPrintf("%d\n",step);
    step++;
    mexPrintf("(%d , %d)\n",mxGetM(PtIn[0]),mxGetN(PtIn[0]));
  */

 switch(NbOut)
   {
   case 2 : break ;
     // 2 is the final paramter
   default : mexErrMsgTxt(MESSAGE_D_ERREUR) ; return ; break ;
   }

  switch(NbIn)
   {
   case 6 : break ;
     // 3 is the final parameter
   default : mexErrMsgTxt(MESSAGE_D_ERREUR) ; return ; break ;
   }


 /*****************************************************
  *      input checking
  *****************************************************/     

  // data checking
  if(mxGetN(PtIn[0]) != 4)
    {
      mexPrintf("The data matrix should be nx4.\n") ; 
      mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
      return ;
    }

  // variogram parameter checking
  if(mxGetN(PtIn[1])*mxGetM(PtIn[1]) != 6 || Min(mxGetN(PtIn[1]),mxGetM(PtIn[1]))==1)
    {
      mexPrintf("The variogram parameter must be 3x2 or 2x3.\n") ; 
      mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
      return ;
    }
  // domain
  if(mxGetM(PtIn[2]) != 2 || mxGetN(PtIn[2]) != 2)
    {
      if(mxGetM(PtIn[2])* mxGetN(PtIn[2]) != 1)
	{
	  mexPrintf("the kriging domain must be of the form 2x2 where each line is the (X,Y) coordinates of two opposite extreme points of the rectangular domain\n") ; 
	  mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
	  return ;
	}
      automatic_domain_definition = true;
    }


  // sampling positions
  if(mxGetN(PtIn[3]) == 2)
    {
      /*
	vector with m=p lines and n=2 columns 
	p corresponds to the number of sampling positions on the line
	column 1 = the X positions
	column 2 = the Y positions
      */
      n_estimation_positions = mxGetM(PtIn[3]) ;
      if(n_estimation_positions == 1)
	estimation_nature = POINT;
      else
	estimation_nature = LINE;
      mexPrintf("%d estimation position(s)\n",n_estimation_positions);
    }
  else if(mxGetN(PtIn[3]) > 2 && mxGetM(PtIn[3]) == 2)
    {
      /*
	vector with m=2 lines and n=p columns 
	p corresponds to the number of sampling positions of the grid
	line 1 = the X positions
	line 2 = the Y positions
      */
      estimation_nature = GRID;
      n_estimation_positions = mxGetN(PtIn[3]);
      // mexPrintf("%d x %d estimation positions on a grid \n",n_estimation_positions,n_estimation_positions);  
    }
  else if (mxGetN(PtIn[3])* mxGetM(PtIn[3]) == 1)
    {
      
      estimation_nature = GRID;
      automatic_grid_definition = true;
      pt = mxGetPr(PtIn[3]);      
      n_estimation_positions = *pt;
      // mexPrintf("%d x %d estimation positions on a grid \n",n_estimation_positions,n_estimation_positions);  
    
    }
  else
    {
      mexPrintf("precise kriging on something else than a point\n") ; 
      mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
      return ;
    }


  if(mxGetN(PtIn[4])*mxGetM(PtIn[4]) != 1)
    {
      mexPrintf("The Method parameter must be an integer (1, 2 or 3).\n") ; 
      mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
      return ;
    }
  
  pt = mxGetPr(PtIn[4]);
  if(! ( (*pt == HYBRID) || (*pt == SIMULATED) || (*pt == COMBINATORIAL) ) )
    {
      mexPrintf("The Method parameter must be an integer (1, 2 or 3).\n") ; 
      mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
      return ;
    }

  if(mxGetN(PtIn[5])*mxGetM(PtIn[5]) != 1)
    {
      mexPrintf("The number of simulations must be an integer.\n") ; 
      mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
      return ;
    }  

  pt = mxGetPr(PtIn[5]);
  if( (*pt) > 1000000 )
    {
      mexPrintf("The number of simulations must be lower than 1000000.\n") ; 
      mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
      return ;
    }

  /*
   *  input retrieval
   */


  n_data = mxGetM(PtIn[0]) ;
  
  krig.setDataNature(2); // intervallist data (1: precise ; 3: fuzzy)
  //krig.setV;
  krig.allocate(n_data,2);
  krig.computeVariogramCloud();

  // data retrieval and transfer to krig
  pt = mxGetPr(PtIn[0]) ;
  fin = pt + n_data ;
  k = 0;
  while(pt<fin)
    {
      gsl_vector_set(current_position,0,*pt) ;
      gsl_vector_set(current_position,1,*(pt+n_data)) ;
      gsl_vector_set(current_data,0,*(pt+2*n_data)) ;
      gsl_vector_set(current_data,1,*(pt+3*n_data)) ;
      krig.setCoordinate(k,current_position);
      krig.setImpreciseData(k,current_data);
      krig.setData(k,(gsl_vector_get(current_data,1)+gsl_vector_get(current_data,0))/2);
      k++;
      pt++;
    }

  //  krig.display();
  

  // variogram parameters retrieval and transfer to krig
  pt = mxGetPr(PtIn[1]) ;
  if(mxGetM(PtIn[1]) == 2)
    {
      gsl_vector_set(current_params,0,*(pt)) ;
      gsl_vector_set(current_params,1,*(pt+1)) ;
      krig.setVariogramParameter(0,((*pt) + (*(pt+1)))/2);// mexPrintf("param 1 : %g\n",((*pt) + (*(pt+1)))/2);
      krig.setImpreciseVariogramParameter(0,current_params);

      gsl_vector_set(current_params,0,*(pt+2)) ;
      gsl_vector_set(current_params,1,*(pt+3)) ;
      krig.setVariogramParameter(1,((*(pt+2)) + (*(pt+3)))/2);
      krig.setImpreciseVariogramParameter(1,current_params);

      gsl_vector_set(current_params,0,*(pt+4)) ;
      gsl_vector_set(current_params,1,*(pt+5)) ;
      krig.setVariogramParameter(2,((*(pt+4)) + (*(pt+5)))/2);// mexPrintf("param 3 : %g\n",((*(pt+4)) + (*(pt+5)))/2);
      krig.setImpreciseVariogramParameter(2,current_params);
    }
else
    {
      
      gsl_vector_set(current_params,0,*(pt)) ;
      gsl_vector_set(current_params,1,*(pt+3)) ;
      krig.setVariogramParameter(0,((*pt) + (*(pt+3)))/2);// mexPrintf("param 1 : %g\n",((*pt) + (*(pt+1)))/2);
      krig.setImpreciseVariogramParameter(0,current_params);

      gsl_vector_set(current_params,0,*(pt+1)) ;
      gsl_vector_set(current_params,1,*(pt+4)) ;
      krig.setVariogramParameter(1,((*(pt+1)) + (*(pt+4)))/2);
      krig.setImpreciseVariogramParameter(1,current_params);

      gsl_vector_set(current_params,0,*(pt+2)) ;
      gsl_vector_set(current_params,1,*(pt+5)) ;
      krig.setVariogramParameter(2,((*(pt+2)) + (*(pt+5)))/2);// mexPrintf("param 3 : %g\n",((*(pt+4)) + (*(pt+5)))/2);
      krig.setImpreciseVariogramParameter(2,current_params);

    }

  // domain parameter retrieval
  if(automatic_domain_definition)
    {
      // mexPrintf("domain automatically defined\n");
      krig.setDomain();
    }
  else
    {
      // mexPrintf("domain user defined\n");
      pt = mxGetPr(PtIn[2]) ;
      gsl_vector_set(lower,0,Min(*pt,*(pt+1)));
      gsl_vector_set(upper,0,Max(*pt,*(pt+1)));
      gsl_vector_set(lower,1,Min(*(pt+2),*(pt+3)));
      gsl_vector_set(upper,1,Max(*(pt+2),*(pt+3)));
      
      krig.setDomain(lower, upper);
   
    }

  
  /*
  //display of the position vector
  krig.getDomain(lower,upper);
  mexPrintf("domain = lower : %g , %g ; upper : %g , %g\n",gsl_vector_get(lower,0),gsl_vector_get(lower,1),gsl_vector_get(upper,0),gsl_vector_get(upper,1));
  */
    
  
  // positions retrieval
  if(estimation_nature == POINT || estimation_nature == LINE )
    {
      positions = gsl_matrix_alloc(n_estimation_positions,2);
      pt = mxGetPr(PtIn[3]) ;
      fin = pt + n_estimation_positions ;
      k = 0;
      while(pt<fin)
	{
	  gsl_matrix_set(positions,k,0,(*pt++)) ;
	  k++;
	}
      k = 0;
      fin = pt + n_estimation_positions;
      while(pt<fin)
	{
	  gsl_matrix_set(positions,k,1,(*pt++)) ;
	  k++;
	}
    }
  else if(automatic_grid_definition)
    {
      positions = gsl_matrix_alloc(2,n_estimation_positions);  
      krig.getDomain(lower,upper);
      Xm = gsl_vector_get(lower,0);
      Xstep = (gsl_vector_get(upper,0)-Xm)/n_estimation_positions;
      Ym = gsl_vector_get(lower,1);
      Ystep = (gsl_vector_get(upper,1)-Ym)/n_estimation_positions;
      for(i=0;i<n_estimation_positions;i++)
	{
	  gsl_matrix_set(positions,0,i,Xm+i*Xstep) ;
	  gsl_matrix_set(positions,1,i,Ym+i*Ystep) ;
	}
    }
  else
    {
      positions = gsl_matrix_alloc(2,n_estimation_positions);  
      pt = mxGetPr(PtIn[3]) ;   
      k = 0;
      fin = pt + 2*n_estimation_positions;
      while(pt<fin)
	{
	  gsl_matrix_set(positions,0,k,(*pt++)) ;
	  gsl_matrix_set(positions,1,k,(*pt++)) ;
	  k++;
	}
    }

  /*
  //display of the position vector
  if(estimation_nature == GRID)
    {
      for(i=0;i<n_estimation_positions;i++)
	{
	  mexPrintf("positions(0,%d) = %g\n", i,gsl_matrix_get(positions,0,i));
	  mexPrintf("positions(1,%d) = %g\n", i,gsl_matrix_get(positions,1,i));
	}
    }
  */

  // method retrieval
  method = *(mxGetPr(PtIn[4])) ;


  // number of simulations retrieval
  //  pt = mxGetPr(PtIn[5]) ;
  krig.setSimulation(*(mxGetPr(PtIn[5])));


 /*****************************************************
  *      processing
  *****************************************************/

  if(estimation_nature == POINT || estimation_nature == LINE )
    {

      inf_predictions = gsl_matrix_alloc(n_estimation_positions,1);
      sup_predictions = gsl_matrix_alloc(n_estimation_positions,1);
      for (size_t i = 0;i < n_estimation_positions;i++)  
	{
	  // target parameter domX
	  gsl_vector_set(current_position, 0, gsl_matrix_get(positions,i,0));
	  gsl_vector_set(current_position, 1, gsl_matrix_get(positions,i,1));	
	  if(method == HYBRID)
	    {
	      krig.getOptimumBound_hybrid(current_prediction, current_position);	       
	    }
	  else if(method == SIMULATED)
	    {
	      krig.getOptimumBound_simulated_annealing(current_prediction, current_position);	       
	    }
	  else
	    {
	      krig.getOptimumBound_bounds_combinatorial(current_prediction, current_position);	       
	    }
	  gsl_matrix_set(inf_predictions,i,0,gsl_vector_get(current_prediction,0));
	  gsl_matrix_set(sup_predictions,i,0,gsl_vector_get(current_prediction,1));
	}
    }
  else {
    inf_predictions = gsl_matrix_alloc(n_estimation_positions,n_estimation_positions);
    sup_predictions = gsl_matrix_alloc(n_estimation_positions,n_estimation_positions);
    for(i=0;i<n_estimation_positions;i++)
      {
	gsl_vector_set(current_position, 0, gsl_matrix_get(positions,0,i));
	for(j=0;j<n_estimation_positions;j++)
	  {
	    gsl_vector_set(current_position, 1, gsl_matrix_get(positions,1,j));	
	    if(method == HYBRID)
	      {
		krig.getOptimumBound_hybrid(current_prediction, current_position);	       
	      }
	    else if(method == SIMULATED)
	      {
		krig.getOptimumBound_simulated_annealing(current_prediction, current_position);	       
	      }
	    else
	      {
		krig.getOptimumBound_bounds_combinatorial(current_prediction, current_position);	       
	      }
	    gsl_matrix_set(inf_predictions,i,j,gsl_vector_get(current_prediction,0));
	    gsl_matrix_set(sup_predictions,i,j,gsl_vector_get(current_prediction,1));
	    
	  }
      }
}


 /*****************************************************
  *     output creation
  *****************************************************/
   

if(estimation_nature == POINT || estimation_nature == LINE )
  {
    PtOut[0] = mxCreateDoubleMatrix(n_estimation_positions,1,mxREAL) ;
    pt = mxGetPr(PtOut[0]) ;
    fin = pt + n_estimation_positions ;

    k = 0;
    while(pt<fin)
      {
	(*pt++) = gsl_matrix_get(inf_predictions,k,0) ;
	k++;
      } 
   
    PtOut[1] = mxCreateDoubleMatrix(n_estimation_positions,1,mxREAL) ;
    pt = mxGetPr(PtOut[1]) ;
    fin = pt + n_estimation_positions ;
    k = 0;
    while(pt<fin)
      {
	(*pt++) = gsl_matrix_get(sup_predictions,k,0) ;
	k++;
      }
   

  }
 else if(estimation_nature == GRID )
   {
     PtOut[0] = mxCreateDoubleMatrix(n_estimation_positions,n_estimation_positions,mxREAL) ;
     pt = mxGetPr(PtOut[0]) ;
     //fin = pt ;
     for (i = 0; i < n_estimation_positions; i++)
     {
       k = 0;
       fin = pt + n_estimation_positions;
       while(pt<fin)
	 {
	   (*pt) = gsl_matrix_get(inf_predictions,i,k) ;
	   pt++;
	   k++;
	 }
     }
  
     PtOut[1] = mxCreateDoubleMatrix(n_estimation_positions,n_estimation_positions,mxREAL) ;
     pt = mxGetPr(PtOut[1]) ;
     fin = pt ;
     for ( i =0; i < n_estimation_positions; i++)
     {
       k = 0;
       fin = fin + n_estimation_positions;
       while(pt<fin)
	 {
	   (*pt++) = gsl_matrix_get(sup_predictions,i,k) ;
	   k++;
	 }
     }
   }
 else
   {
     mexErrMsgTxt(MESSAGE_D_ERREUR) ; 
     return ;
   }


 // Recopie des probabilites calculees
/*
 pt = mxGetPr(PtOut[0]) ;
 fin = pt + N ;
 ptt = DensiteDeProbabilite ;

 while(pt<fin)
 {
  (*pt++) = (*ptt++) ;
 }

 delta = ( borne_sup - borne_inf ) / (double)( N - 1 ) ;
 pt = mxGetPr(PtOut[1]) ;
 fin = pt + N ;
 (*pt++) = borne_inf ;
 while(pt<fin)
 {
  (*pt) = (*(pt-1)) + delta ; pt++ ;
 }
*/
 

 gsl_matrix_free(positions) ;
 gsl_matrix_free(inf_predictions) ;
 gsl_matrix_free(sup_predictions) ;
 gsl_vector_free(current_position);
 gsl_vector_free(current_data);
 gsl_vector_free(current_params);
 gsl_vector_free(current_prediction);
 gsl_vector_free(upper);
 gsl_vector_free(lower);
 
}

