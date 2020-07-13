#ifndef Solver_H_
#define Solver_H_


#include "commonData.h"
#include <fstream>
#include <string>
//#include <cusparse_v2.h>
class Solver {
	public :
    vector<double> solve3Diag (const vector <double> & lDiag, 
	                           const vector <double> & Diag, 
							   const vector <double> & uDiag,
							   const vector <double> & rHS) ;

    vector<double> SOR3DiagPeriodic (const vector  <bool>  & nodeIsActive,
									 const vector <double> & lDiag, 
	                                 const vector <double> & Diag, 
								     const vector <double> & uDiag,
								     const vector <double> & rHS,
								     const vector <int>    & prevIndex,
								     const vector <int>    & nextIndex,
									 vector <double> & firstGuess) ; 

}; 

#endif
