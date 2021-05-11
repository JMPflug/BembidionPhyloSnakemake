#ifndef MATH_OP_H
#define MATH_OP_H

#include <climits>

using namespace std;

/***********************************************************************************			

function binom
--------------   

Calculates n choose k.The algorithm is taken from: Foundations of Algorithms in 
C++ Pseudocode (Neapolitan,Namipour) and uses dynamic programming.

Parameters:

int n
int k

Return values:

long int : n choose k

************************************************************************************/

long int binom(int n, int k);

/***********************************************************************************

function geometry_mapping
-------------------------

The function gets references to 4 sequences (string a-d), a reference to a 
substitution-matrix and calculates the barycentric coordinates by geometry mapping
(see Nieselt-Struwe, v. Haeseler, 2001). 

Parameters:

string& a, string& b, string& c, string&d 	: quartet sequences
int (*s)[CHAR_MAX]				: scorematrix
double *barycent				: barycent contains the calculated 
						  barycentric coordinates

Return values:

void

************************************************************************************/

void geometry_mapping(string& a, string& b, string& c, string& d, 
		      int (*s)[CHAR_MAX], double *barycent );


/***********************************************************************************

function euclid_dist
--------------------

The function calculates the Euclidean distance for 2 vectors.

Parameters:

float* veca : vector a
float* vecb : vector b
int vecsize : number of elements of vector

Return values:

void

************************************************************************************/
double euclid_dist(double* veca, double* vecb, int vecsize);

#endif
