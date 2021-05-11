/**************************************************************************  
 * 									  * 
 * MARE (Matrix Reduction) 2009						  *
 * 									  *
 * Department of Entomology 						  *
 * Biozentrum Grindel & Zoologisches Museum Hamburg 2009		  *
 *									  * 
 *  									  *
 **************************************************************************/

/*
  BUG-LOG
	
	1. overwork of matrix.dat and matrix_red.dat
	2. matrix.dat renamed to matrix.txt
	3. taxon names read and constraints checked before evaluation of partitions	
	4. matrix_red.dat renamed to matrix_red.txt
	5. reduced and unreduced matrices as well as simplex graphs can optionally 
	   be created as *.svg
	6. taxon weight can be chosen by user
	7. info weight can be chosen by user
	8. in mare.cpp reduce_matrix(): new code added to put constraints on genes
	9. matrix read from file includes last two columns
	10. bugs fixed concerning linefeeds of the input files
	11. color code of the matrix output modified


  TOFIX BUG-LOG

	1. taxon weights < 1
*/


#include <iostream> 
#include <fstream>
#include <iomanip>
#include <vector>
#include <climits>
#include <float.h>
#include <cstdlib>
#include <limits>
#include <sys/stat.h>
#include "mare_a.h"
#include "math_op.h"
#include "matrices.h"
#include "error_handling.h"
#include <cerrno>
#include <cstring>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
  opt_help (argv, argc);

  if(argc<3)
    {
     cerr << "Not enough arguments! MARE requires two files as arguments.\n" 
	  << "For further help type ./MARE -h\n"; 
     exit(1);
    }

//This version of MARE only works with files of Unix format
  if(check_linefeed(argv[1])!=1)
    {
     cerr << "Incorrect format of input file!" << argv[1] << " Unix linefeeds required.\n";
     exit(1);
    }

  if(check_linefeed(argv[2])!=1)
    {
     cerr << "Incorrect format of input file!" << argv[2] << " Unix linefeeds required.\n";
     exit(1);
    }

//Checking format of FASTA file
  int part_end_fas = check_fas(argv[2]);

//Checking format of charset file
  int part_end_cha = check_part(argv[1]);

//Check if both files return the same alignment length
  if(part_end_fas!=part_end_cha)
    {
     cerr << "Incorrect format of input files! Files define different alignment lengths.\n";
     exit(1);
    }

//vector taxon_name contains the names of the taxa extracted from the
//concatenated multiple sequence alignment.
  
  vector<string> taxon_name;

  get_taxa(argv[2],taxon_name);

//vector fixed_taxa contains the numbers of the fixed taxa
  vector<int> fixed_taxa;

//opt_taxa_constraints reads the names of the taxa that have 
//to remain in the matrix 
  opt_taxa_constraints(argv, argc, taxon_name, fixed_taxa);

  double  iw = opt_info_weighting(argv,argc); 

//If option -t is used tax_weight is equal to 5/3
  double tax_weight  = opt_taxa_weighting(argv,argc);

  if (tax_weight < 1) 
     {
       perror("Please enter a value > 1 when choosing taxon weighting option -t");
       exit(errno);
     }

//taxon_name.clear();

//vector partition_bound stores start and end of each partition from the
//concatenated msa. start_end is a struct consisting of 2 integers 
//for start and end of a partition.
  vector<start_end> partition_bound;

//vector part_name contains the names of the partitions
  vector<string> part_name;

//get_partitions extracts the start and end of a partition
  get_partitions(argv[1], partition_bound, part_name); 

//vector fixed_genes contains the numbers of the fixed genes
  vector<int> fixed_genes;

//opt_gene_constraints reads the names of the genes that have 
//to remain in the matrix 
  opt_gene_constraints(argv, argc, part_name, fixed_genes);


//vector matrix represents the weighted matrix. It is composed of pres_abs 
//vectors. The pres_abs vector has as many entrys as partitions are present. 
//An entry represents the evaluation of a partition. The entries are initiated
//with -1 at first.
  vector<double> pres_abs (partition_bound.size()+2,-1); 
  vector< vector <double> > matrix;

//vector partition_seq contains the available sequences of a partition
  vector<string> partition_seq;

//barycent contains 20000 by weighted geometry mapping 
//calculated barycentric coordinates.
  double barycent[NUM_OF_QUA][3];

//scorematrix contains the substitution values after initiation. CHAR_MAX
//is the maximal value of char, usually 255.
  int scorematrix [CHAR_MAX][CHAR_MAX];

//initiation of the blosum matrix
  init_blosum62(scorematrix);

//vector all_quart holds NUM_OF_QUA quartets drawn randomly from partition_seq.
//quart is a struct that contains 4 integers representing 4 different taxa.
  vector <quart> all_quart;

//A directory is created to store the quartet maps

  if(mkdir("results",S_IRWXU)==-1)
    if(errno != EEXIST)
      {
       perror("Error! Cannot create new directory: ");
       exit(errno);
      }

 
opt_reuse_matrix(argv, argc, matrix, pres_abs);


if(matrix.size()==0)
	{
	cout << "\nCalculation of potential information content:\n\n" << flush;

//The following loop fills each column of the matrix with -1 (sequence not available) 
//and 1 (available) by calling  function pres_abs_matrix. If a partition is extracted, 
//20000 randomly drawn sequence-quartets are evaluated by weighted geometry-
//mapping(get_braycent(_rand)). Afterwards evaluate_partition determines the
//information content of a sequence.  

  	for(int y=0;y<partition_bound.size();y++)
	{
		cout << "working on       : " << part_name[y] << "\n" << flush;
		pres_abs_matrix(partition_bound[y].start,partition_bound[y].stop,argv[2],
		matrix, pres_abs,y,partition_seq, taxon_name); 
		int n = partition_seq.size();
		double rel_info = 0;	

		if(n>=4) 
	        { 
			long int bc = binom(n,4);

			if(bc<=NUM_OF_QUA)
			{
				get_barycent(n,partition_seq,scorematrix,barycent);
				rel_info = evaluate_partition(barycent,bc,matrix,y);

				if(opt_pics (argv, argc))
					draw_quartetmap (y, barycent, bc, part_name[y], rel_info);
			}
			else
           		{
				get_barycent_rand(n, partition_seq, scorematrix,barycent);
				rel_info = evaluate_partition(barycent,NUM_OF_QUA,matrix,y);

				if(opt_pics (argv, argc))
					draw_quartetmap (y, barycent, NUM_OF_QUA, part_name[y], rel_info);
			}
		}
		else
	        {
//sequences in partitions with less than 4 sequences are evaluated with 0
			for(int i=0;i<matrix.size();i++)
			{
				if(matrix[i][y]==1) matrix[i][y] = 0;
			}
		}

//print relative information: 
		cout << "rel. Information : " <<  rel_info << "\n\n" << flush;
		partition_seq.clear();
	}
	}

/*
//vector fixed_taxa contains the numbers of the fixed taxa
vector<int> fixed_taxa;

//opt_taxa_constraints reads the names of the taxa that have 
//to remain in the matrix 
opt_taxa_constraints(argv, argc, taxon_name, fixed_taxa);
*/


//Two lines are added to the matrix, so evaluations for every 
//partition and the number of the partition can be stored: 
  matrix.push_back(pres_abs);
  matrix.push_back(pres_abs);

  int p = matrix.size();
  int q = matrix[0].size();

  int j = 0;

  int n =  matrix.size()-2;
  int m =  matrix[0].size()-2;

//Numbering of the taxa:
  for(int i=0;i<n;i++)
     {
       matrix[i][m]=i;
     }

//Evaluation of the taxa:
  eval_taxa(matrix);

//Numbering of the partitions:
  for(int j=0;j<m;j++)
     {
       matrix[n][j]=j;
     }

//Evaluation of the partitions:
  eval_gene(matrix);


//write the matrix to stdout

  j = 0;

  string out1 = "results/";
  out1.append(argv[2]);
  out1.append("_matrix.txt");
  char * cstr_out1;
  cstr_out1 = new char [out1.size()+1];
  strcpy (cstr_out1, out1.c_str());


//write the matrix to matrix.txt
  ofstream fileout1 (cstr_out1,ios_base::out);
  if(!fileout1.good())
    {
     cerr << "Cannot open matrix.txt!\n";
     exit(1);
    }

  j = 0;

  fileout1 << setw(50) << "" << " \t" ;

  while(j < part_name.size())
       {
	fileout1 << setw(20) << left << part_name[j] << "\t";
        j++;
       }
  fileout1 << "\n";

  for(j=0;j<p-2;j++)
     {
      fileout1 << setw(50) << taxon_name[j] << " \t" ;
      for(int i=0;i<q-2;i++)
         {
          fileout1 << setw(20) << setfill(' ') << setprecision(5) << left  << matrix[j][i] << "\t";
         }
      fileout1 << "\n";
     }
  fileout1.close();


//Presence absence matrix, added for simulations
  vector< vector <double> > matrix_pres_abs (matrix);
  float entry=0.0;

  for(j=0;j<p-2;j++)
     {
      for(int i=0;i<q-2;i++)
         {
          if(matrix[j][i]>=0) 
		{
		 matrix_pres_abs[j][i]=1;
		 entry++;
		}
         }
     }


//Sorting of the matrix:
  sort_matrix(matrix_pres_abs,fixed_taxa);


/*
  string filename_pa = "results/pres_abs_matrix_sort.svg";
  draw_matrix(matrix_pres_abs, filename_pa, taxon_name, part_name);
*/

//Sorting of the matrix:
  sort_matrix(matrix,fixed_taxa);

//  Usually matrices are too large to output them as *.svg

if(opt_matrix (argv, argc))
  {
     string filename_dat = "results/";
	    filename_dat.append(argv[2]);
	    filename_dat.append("_matrix_unred_sort.svg");
   draw_matrix(matrix, filename_dat, taxon_name, part_name);
  }

//original matrix size is stored in oms 
//current matrix size is stored in cms
  double oms = (p-2)*(q-2);
  double cms = (p-2)*(q-2); 

//a stores the current result of the optimality function. ap is the 
//average potential information of the matrix. f stores the current 
//optimum
  double a, ap;
  double f = 0;

//matrix_pot_info returns the average information content
  ap =  matrix_pot_info(matrix);

//If all entries are evaluated with 1 the complete matrix is optimal
  if(ap==1)
    {
     cerr << "The input matrix is optimal. There is no reason to reduce it!\n";
     return 0;
    }

//reduction counts the reduction steps
  int reduction = 0;

//Before the reduction starts f is calculated 
//cluster_taxa checks the connectivity of the matrix


  if(cluster_partitions(matrix,LOWER_LIM))
    {
//    if(iw!=0) 	f = optimality_f(cms, oms, ap, iw);
	 	f = optimality_f(cms, oms, ap, iw);
//    else		f = unw_optimality_f(cms, oms, ap);
    }

//Now the reduction of the matrix begins. reduce_matrix returns 0 if 
//there are <4 taxa and <2 sequences left. The number of the taxa and 
//partitions, when an optimum is found are stored in partitions and
//taxa
  vector <double>  partitions;
  vector <double>  taxa;

  string out2 = "results/";
  out2.append(argv[2]);
  out2.append("_plot.txt");
  char * cstr_out2;
  cstr_out2 = new char [out2.size()+1];
  strcpy (cstr_out2, out2.c_str());


//write plot.txt
  ofstream fileout2 (cstr_out2,ios_base::out);
  if(!fileout2.good())
    {
     cerr << "Cannot open plot.txt!\n";
     exit(1);
    }

  fileout2 << "#p      |lam     |f       |parts   |taxa    |clus    |red     |\n"
           << "#-------+--------+--------+--------+--------|--------|--------|\n";

  p   = matrix.size();
  q   = matrix[0].size();
  cms = double(p-2)*double(q-2);
  ap  = matrix_pot_info(matrix);
  int noc = cluster_partitions(matrix,LOWER_LIM);

  fileout2 << " " << setw(8) << setprecision(3) << left << ap << setw(9) << setprecision(3) 
           << left << cms/oms << setw(9) << setprecision(3) << left << f << setw(9) << left 
	   << q-2 << setw(9) << left << p-2 << setw(9) << left << noc << setw(9) << left 
	   << reduction << "\n";

//The following variables will contain the results of the output matrix
  double ap_res = ap;
  double ma_rat = cms/oms;
  int no_tax   = p-2;
  int no_par   = q-2;


  string out3 = "results/";
  out3.append(argv[2]);
  out3.append("_info.txt");
  char * cstr_out3;
  cstr_out3 = new char [out3.size()+1];
  strcpy (cstr_out3, out3.c_str());

//write info.txt
  ofstream fileout3 (cstr_out3,ios_base::out);
  if(!fileout3.good())
    {
     cerr << "Cannot open info.txt!\n";
     exit(1);
    }
  
fileout3 << "MARE OUTPUT\n"
	 << "-----------\n\n";

vector<string> command;

fileout3 << "Command:\n";

for(int i=0;i < argc ;i++)
	{
	   fileout3 << argv[i] << " ";
	}
	   fileout3 << "\n\n";

fileout3   << "original multiple sequence alignment\n"
	   << "------------------------------------\n"
	   << "information content of B     : " << setprecision(3) << setw(4) << left 
	   << ap_res << "\n"
	   << "matrix saturation            : " << entry/(no_tax*no_par) << "\n"
	   << "nr of taxa                   : " << no_tax << "\n"
	   << "nr of partitions             : " << no_par << "\n\n";


//Reduction of the matrix. The matrix is reduced until there is only one partition and
//three taxa left. Then reduce_matrix returns 0

  cout <<"Reduction of the matrix..." << flush;

//Copy the current matrix
  vector < vector <double> > res_matrix(matrix);

  while(reduce_matrix(matrix, fixed_taxa,fixed_genes, tax_weight))
       {
        p   = matrix.size();
        q   = matrix[0].size();
        cms = double(p-2)*double(q-2);
        ap  = matrix_pot_info(matrix);

//	cout << f << "\n" << flush;
//	cout << a << "\n" << flush;


//        if(iw!=0)
            a = optimality_f(cms, oms, ap, iw);
//        else
//          a = unw_optimality_f(cms, oms, ap);

        noc = cluster_partitions(matrix,LOWER_LIM);
        reduction++;
	
	fileout2 << " " << setw(8) << setprecision(3) << left << ap << setw(9) 
		 << setprecision(3) << left << cms/oms << setw(9) << setprecision(3) 
		 << left << a << setw(9) << left << q-2 << setw(9) << left << p-2 
		 << setw(9) << left << noc << setw(9) << left << reduction << "\n";

//If the number of clusters is 1, the result of the optimality function is checked
        if(noc==1)
          {

//If the current result of the optimality function is better, get the current
//partitions and taxa
           if(a>f)
	     {
	      f = a;
	      ap_res = ap;
	      ma_rat = cms/oms;
	      no_tax = p-2;
	      no_par = q-2;

//Get the current partition:
	      partitions.resize(q-2);
	      for(int i=0; i<q-2; i++)
	         {
	          partitions[i] = matrix[p-2][i];
	         }

//Get the current taxa:
	      taxa.resize(p-2);
	      for(int i=0; i<p-2; i++)
	         {
	          taxa[i] = matrix[i][q-2]; 
	         }

//Copy of the current optimal matrix for results
	      res_matrix.resize(matrix.size());
	      entry = 0;
	      for(int j=0;j<p;j++)
     		 {
	          for(int i=0;i<q;i++)
         	     {
		      res_matrix[j].resize(matrix[j].size());
          	      res_matrix[j][i] = matrix[j][i];
			if ((matrix[j][i]>=0) && (j < p-2) && (i < q-2))entry++;
		     }
 		 }
	     }
          }
       }


  fileout2.close();

//The resulting sequences are extracted from the original alignment
//and stored in a new one
 get_fastaout(argv[2], partition_bound, partitions, taxa);

//A new charset is made
  get_charsetout(argv[1], partition_bound, partitions);

  fileout3 << "Parameters\n"
  	   << "----------\n";
//  if(iw!=0)
  fileout3 << "weighting of inform. content : " << iw << "\n";
//  else  
//  fileout3 << "weighting of inform. content : unweighted\n";
  fileout3 << "weighting of taxa            : " << setprecision(3)
	   << tax_weight << "\n"
  	   << "number of evaluated quartets : " << NUM_OF_QUA << "\n\n"

  	   << "reduced multiple sequence alignment\n"
  	   << "-----------------------------------\n"
  	   << "information content of B'            : " << setw(4) << setprecision(3)
	   <<  ap_res << "\n"
	   << "nr of taxa                           : " << no_tax << "\n"
	   << "nr of partitions                     : " << no_par << "\n"
	   << "ratio of reduced and complete matrix : " << setw(4) << setprecision(3)
	   << ma_rat << "\n"
	   << "matrix saturation                    : " << entry/(no_tax*no_par) << "\n";
  fileout3.close();

//Evaluation of the partitions:
  eval_gene(matrix);

//write the matrix to matrix_red.txt

  string out4 = "results/";
  out4.append(argv[2]);
  out4.append("_matrix_red.txt");
  char * cstr_out4;
  cstr_out4 = new char [out4.size()+1];
  strcpy (cstr_out4, out4.c_str());

//write matrix_red.txt
  ofstream fileout4 (cstr_out4,ios_base::out);
  if(!fileout4.good())
    {
     cerr << "Cannot open matrix_red.txt!\n";
     exit(1);
    }

  j = 0;

  fileout4 << setw(50) << "" << " \t" ;

  while(j < partitions.size())
       {
	fileout4 << setw(20) << left << part_name[int(partitions[j])] << "\t";
        j++;
       }
  fileout4 << "\n";

  for(j=0;j<taxa.size();j++)
     {
      fileout4 << setw(50) << taxon_name[int(taxa[j])] << " \t";
      for(int i=0;i<res_matrix[0].size()-2;i++)
         {
          fileout4 << setw(20) << setfill(' ') << setprecision(5) << left  << res_matrix[j][i] << "\t";
         }
      fileout4 << "\n";
     }
  fileout4.close();


  if(opt_matrix(argv, argc))
    {
     string filename = "results/";
	    filename.append(argv[2]);
	    filename.append("_matrix_red_sort.svg");
     draw_matrix(res_matrix, filename, taxon_name, part_name);
    }

  cout << "ok!\n\n";
  cout << "done!\n\n";

  return 0;
} 
