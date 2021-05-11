#ifndef MARE_H
#define MARE_H

#include <vector>

#define NUM_OF_QUA   20000   //number of quartets drawn randomly
#define LOWER_LIM    3	     //connectivity between sequences and taxa


using namespace std;


//start_end stores the beginning and the end of a partition
  typedef struct{
    int start;
    int stop;
   }start_end; 


//quart stores 4 integers that represent 4 taxa
  typedef struct{
    int a;
    int b;
    int c;
    int d;
  }quart;


/***********************************************************************************

function get_taxa
-----------------   

The function extracts taxon names.

Parameters:

char* filename    		: filename of the MSA

vector<string>& taxon_name      : names of the taxa

Return values:

void
************************************************************************************/

void get_taxa(char* filename, vector<string>& taxon_name);

/***********************************************************************************

function pres_abs_matrix
------------------------   

The function extracts a partition from the msa and fills the matrix with presence-
absense values. If a taxon has a sequence, the corresponding entry is 1, 
otherwise -1. Only sequences with less than 1/3 gaps are excepted.  


Parameters:

int start	  : start of the partition 	
int stop 	  : end of the partition
char* filename    : filename of the MSA
int i		  : holds the number of the partition

vector< vector <double> >& matrix: 	For every new taxon in the msa the matrix is 
					extended by a vector of type pres_abs. It 
					gets a 1 as an entry if the sequence of the 
					taxon is available, otherwise -1.
 
vector<double> pres_abs	   	:   	pres_abs contains the values of one taxon. At 
					first 1 for an available sequence and -1 if it 
					is not	available.

vector<string>& partition_seq   :	If a sequence exists in the current partition 
					it is added to partition_seq 

vector<string>& taxon_name      :	names of the taxa

Return values:

void
************************************************************************************/


void pres_abs_matrix(int start,int stop,char* filename,
		     vector< vector <double> >& matrix, vector <double> pres_abs, 
		     int i, vector<string>& partition_seq, 
		     vector<string>& taxon_name);


/***********************************************************************************

function get_partitions
-----------------------   

get_partitions extracts start and end of a partition in the msa and the identifier 
of a partition from the charset file.

Parameters:

char* filename				: filename of charset file
vector<start_end>& partition_bound	: partition_bound is a vector of structs of 
					  type start_end 
vector<string>& part_name		: partition name

Return values:

void

************************************************************************************/

void get_partitions(char *filename,vector<start_end>& partition_bound,
		    vector<string>& part_name);

/***********************************************************************************

function get_barycent
---------------------   

function get_barycent is called if there are less than NUM_OF_QUA quartets to calcuate.
It calculates all possible quartets by calling function quartet_mapping. The vector 
barycent is filled with barycentric coordinates.  

Parameters:

int n				  : result of n choose 4, with n = number of taxa.
vector<string>& partition_seq     : partition_seq holds the sequences of a partition
int (*scorematrix)[CHAR_MAX]   	  : scorematrix can be BLOSUM62 or some other sub-
				    stitution matrix 
double (*barycent)[3]	          : barycent contains the 3 barycentric coordinates of
				    the evaluated quartets of one partition.
Return values:

void

************************************************************************************/

void get_barycent(int n,vector<string>& partition_seq, int (*scorematrix)[CHAR_MAX],
		  double (*barycent)[3]);

/***********************************************************************************

function get_barycent_rand
--------------------------   

function get_barycent_rand is called if there are more than NUM_OF_QUA quartets 
to calcuate. It calculates NUM_OF_QUA randomly drawn quartets by calling function 
quartet_mapping. The vector barycent is filled with barycentric coordinates.  

Parameters:

int n				: number of taxa.
vector<string>& partition_seq   : partition_seq holds the sequences of a partition
int (*scorematrix)[CHAR_MAX]   	: scorematrix can be BLOSUM62 or some other sub-
				  stitution matrix 
double (*barycent)[3]	        : barycent contains the 3 barycentric coordinates of
				  the evaluated quartets of one partition.

Return values:

void

************************************************************************************/

void get_barycent_rand(int n,vector<string>& partition_seq, 
		       int (*scorematrix)[CHAR_MAX], double (*barycent)[3]);

/***********************************************************************************


function evaluate_partition
---------------------------   

The function evaluates the partition by counting all vectors (barycentric coordinates) 
which are closest to the barycentric coordinates of Tree 1, 2, 3, 12, 13 and  23. The 
value is divided by the number of all vectors.

Parameters:

(double (*barycent)[3]		  : barycentric coordinates
long int num_of_vec		  : number of vectors
vector< vector <double> >& matrix : presence absence matrix
int y				  : number of partition

Return values:

double : Evaluation of partition

************************************************************************************/

double evaluate_partition(double (*barycent)[3], long int num_of_vec, 
			  vector< vector <double> >& matrix, int y);

/***********************************************************************************

function eval_taxa
------------------   

eval_taxa sums up the weighted entries of the matrix for a taxon. This value 
is divided by the number of partitions. The results is stored in the last column 
of the matrix. 

Parameters:

vector< vector <double> >& matrix : matrix of taxa and genes with weighted entries.

Return values:

void

************************************************************************************/

void eval_taxa(vector< vector <double> >& matrix);

/************************************************************************************

function eval_gene
------------------   

eval_gene sums up the weighted entries of the matrix for a partition. This 
value is divided by the number of partitions. The result is stored in the last line
of the matrix.

Parameters:

vector< vector <double> >& matrix : matrix of taxa and genes with weighted entries.

Return values:

void

************************************************************************************/

void eval_gene(vector< vector <double> >& matrix);

/************************************************************************************

function sort_matrix
--------------------   

The function sorts the matrix, so that the highest evaluated taxon is stored in the 
first column second highest in the second column etc. The highest evaluated gene is 
stored in the first line etc.    

Parameters:

vector< vector <double> >& matrix : matrix of taxa and genes with weighted entries.
vector<int>& fixed_taxa		  : fixed taxa

Return values:

void

************************************************************************************/

void sort_matrix(vector< vector <double> >& matrix, vector<int>& fixed_taxa);

/************************************************************************************

function cluster_taxa
---------------------   

If a taxon has 'lower_limit' sequenes in common with at least one other taxon of a 
Cluster, they are part of the same Cluster. The function returns the number of Clusters 
obtained from the matrix.

Parameters:

vector< vector <double> >& matrix : matrix of taxa and genes  with weighted entries
int lower_limit		     	  : see above 

Return values:

int : number of Clusters

************************************************************************************/

int cluster_taxa(vector< vector <double> >& matrix, int lower_limit);

/************************************************************************************

function cluster_partitions
---------------------------  

If a partition has 'lower_limit' sequenes in common with at least one other partition of a 
Cluster, they are part of the same Cluster. The function returns the number of Clusters 
obtained from the matrix.

Parameters:

vector< vector <double> >& matrix : matrix of taxa and genes  with weighted entries
int lower_limit		     	  : see above 

Return values:

int : number of Clusters

************************************************************************************/

int cluster_partitions(vector< vector <double> >& matrix, int lower_limit);

/************************************************************************************

function reduce_matrix
----------------------

reduce matrix finds the lowest evaluated partition or taxon and drops it from the 
matrix. In case of ties partitions are dropped. Every sequence which does not 
contain >4 taxa is dropped. The function reduces matrix as long as 3 taxa and 1 
partition is left.

Parameters:

vector< vector <float> >& matrix : matrix of taxa and genes with weighted entries
vector<int>& fixed_taxa		 : fixed taxa

Return values:

int 1: more than 3 taxa and 1 partition are left.
    0: 3 taxa and 1 partition are left.

************************************************************************************/

int reduce_matrix(vector< vector <double> >& matrix, vector<int>& fixed_taxa, vector<int>& fixed_genes, 
		  double tax_weight);

/************************************************************************************

function matrix_pot_info
------------------------

The function gets the average potential information content of the matrix

Parameters:

vector< vector <double> >& matrix : matrix of taxa and genes with weighted entries

Return values:

double : average potential information content of the matrix

************************************************************************************/

double matrix_pot_info(vector< vector <double> >& matrix);

/************************************************************************************

function optimality_f
---------------------

Optimality function f of MARE

Parameters:

int cms   	   : current matrix size
int oms   	   : original matrix size
double ap  	   : average potential information
double		   : weighting of information content

Return values:

double : result of optimality function f

************************************************************************************/

double optimality_f(double cms, double oms, double ap, double iw);

/************************************************************************************

function unw_optimality_f
-------------------------

unweighted Optimality function f of MARE

Parameters:

int cms   : current matrix size
int oms   : original matrix size
float ap  : average potential information

Return values:

double : result of optimality function f

************************************************************************************/

double unw_optimality_f(double cms, double oms, double ap);

/************************************************************************************

function get_fastaout
---------------------

get_output extracts the sequences from the original MSA and generates a new FASTA 
file 

Parameters:

char* filename		        : filename of original MSA
vector <double> partition_bound : partitionbounds
vector <double> partitions      : number of partitions
vector <double> taxa	        : number of taxa of the reduced MSA

Return values:

void

************************************************************************************/

void get_fastaout(char* filename, vector<start_end> partition_bound, 
		  vector <double> partitions, vector <double> taxa);

/************************************************************************************

function get_charsetout
-----------------------

get_charsetout extracts the partition bounds from the original charset file and
a new charset file

Parameters:

char* filename                  : filename of original MSA
vector <double> partition_bound : partitionbounds
vector <double> partitions      : number of partitions of the reduced MSA

Return values:

void

************************************************************************************/

void get_charsetout(char* filename, vector<start_end> partition_bound, 
		    vector <double> partitions);


/************************************************************************************

function opt_taxa_constraints
-----------------------------

opt_taxa_constraints is an option to fix taxa while reduction. Fixed taxa are not 
discarded while reduction.

Parameters:

char** argv		 : pointer to argument array
int argc		 : number of command line arguments
vector <char> taxon_name : taxon_name contains the names of the fixed taxa
vector <int> fixed_taxa  : number of fixed taxa

Return values:

void

************************************************************************************/

void opt_taxa_constraints (char** argv, int argc, vector<string> taxon_name, 
			   vector<int>& fixed_taxa);

/************************************************************************************

function opt_gene_constraints
-----------------------------

opt_gene_constraints is an option to fix genes while reduction. Fixed genes are not 
discarded while reduction.

Parameters:

char** argv		 : pointer to argument array
int argc		 : number of command line arguments
vector <char> part_name  : part_name contains the names of the fixed genes
vector <int> fixed_gene  : number of fixed genes

Return values:

void

************************************************************************************/

void opt_gene_constraints (char** argv, int argc, vector<string> part_name, 
			   vector<int>& fixed_gene);

/************************************************************************************

function opt_taxa_weighting
---------------------------

opt_taxa_weighting weights the taxa if the option is used. The user should enter a 
number that is larger than or equal to 0. 

Parameters:

char** argv		 : pointer to argument array
int argc		 : number of command line arguments

Return values:

double weighting of taxa

************************************************************************************/

void opt_reuse_matrix (char** argv, int argc, vector< vector <double> >& matrix, vector<double> pres_abs);

/************************************************************************************

function opt_reuse_matrix
-------------------------

opt_reuse_matrix reuses a precalculated matrix

Parameters:

char** argv		 		: pointer to argument array
int argc		 		: number of command line arguments
vector< vector <double> >& matrix	: For every new taxon in the msa the matrix is 
                                          extended by a vector of type pres_abs. It 
                                          gets a 1 as an entry if the sequence of the 
                                          taxon is available, otherwise -1.
 
vector<double> pres_abs         	: pres_abs contains the values of one taxon. At 
                                          first 1 for an available sequence and -1 if it 
                                          is not  available.
Return values:

void

************************************************************************************/

double opt_taxa_weighting(char** argv, int argc);

/************************************************************************************

function opt_info_weighting
---------------------------

opt_info_weighting weights the information content if the option is used. 
The user should enter a number that is larger than or equal to 0.

Parameters:

char** argv		 : pointer to argument array
int argc		 : number of command line arguments

Return values:

double weighting of information content

************************************************************************************/

double opt_info_weighting(char** argv, int argc);

/************************************************************************************

function opt_help
-----------------

The option outputs a small text of how to use MARE

Parameters:

char** argv		 : pointer to argument array
int argc		 : number of command line arguments

Return values:

void

************************************************************************************/

void opt_help (char** argv, int argc);

/************************************************************************************

function opt_pics
-----------------

The option outputs the simplex graphs as *.svg.
Parameters:

char** argv		 : pointer to argument array
int argc		 : number of command line arguments

Return values:

void

************************************************************************************/

double opt_pics (char** argv, int argc);

/************************************************************************************
function opt_matrix
-------------------

The option outputs the reduced and unreduced matrix as *.svg.
Parameters:

char** argv		 : pointer to argument array
int argc		 : number of command line arguments

Return values:

void

************************************************************************************/

double opt_matrix (char** argv, int argc);

/************************************************************************************

function draw_quartetmap
------------------------

draw_quartetmap generates an *svg file It contains the simplex graph with the 
barycentric coordinates of a partition

Parameters:

int num_of_vec	 	 : number of vectors
double (*barycent)[3]	 : barycentric coordinates
int y		 	 : number of partitions
string part_name	 : partition identifier
double rel_info		 : relative information of partition

Return values:

void

************************************************************************************/

void draw_quartetmap (int y, double (*barycent)[3], int num_of_vec, string part_name,
		      double rel_info);

/************************************************************************************

function draw_matrix
--------------------

draw_matrix generates an *svg file. It contains the matrix passed to the function. 

Parameters:

vector< vector <double> >& matrix : matrix of taxa and genes with weighted entries
string filename			  : filename
vector<string> taxon_name	  : names of taxa
vector<string> part_name	  : names of partitions

Return values:

void

************************************************************************************/

void draw_matrix(vector< vector <double> >& matrix, string filename, 
		 vector<string> taxon_name, vector<string> part_name);

#endif

