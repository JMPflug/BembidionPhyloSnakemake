#ifndef MATRICES_H
#define MATRICES_H

#include <vector>

using namespace std;

/***********************************************************************************
function init_blosum62
----------------------

The function fills blosum62 substitution values into a (CHAR_MAX x CHAR_MAX)-
matrix to allow a fast access to the substitution-values.

Parameters:

int (*scorematrix)[CHAR_MAX] : This is an 2D array with values for all possible 
			       substitutions.

Return values:

void

************************************************************************************/

void init_blosum62(int (*scorematrix)[CHAR_MAX]);

#endif
