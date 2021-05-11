#include <iostream>
#include <fstream> 
#include <sstream>
#include <vector>
#include <string>
#include <climits>
#include <cfloat>
#include <ctime>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <locale>
#include <cstring>
#include <algorithm>
#include "mare_a.h"
#include "math_op.h"

using namespace std;

/*************************************************************************************************

fixed bugs
----------

11.07.09: line 297-304, 217-229; Sequences have to get shuffled, to avoid a pattern when mapping 
	  the quartet



**************************************************************************************************/


/*************************************************************************************************/

void get_partitions(char* filename,vector<start_end>& partition_bound, 
		    vector<string>& part_name)
{
  string buf;
  string name;
  char str[10];
  int k=0;
  start_end se;

  ifstream file(filename, ios::in);
  if(!file.good())
    {
     cerr << "Cannot open file " << filename << " ! Check if the filename "
          << "is written correctly.\n";
     exit(1);
    }

  while(!file.eof())
       {
        memset(str,'\0',10); 
	buf.empty();
        name.clear(); 
	int j=0;

//If a line is empty, it is ignored
        getline(file,buf);
	if(buf.empty()) continue;

//Translate buf to a string with only uppercase letters	
        locale loc;
        for (size_t l=0; l<buf.length(); ++l)
	    {
             name.push_back(toupper(buf[l],loc));
	    }

//Find the beginning of substring CHARSET	
	string substr ("CHARSET");
	size_t pos = name.find(substr);

//i marks the position after CHARSET in buf.
	int i = int(pos)+substr.size();
	part_name.resize(part_name.size()+1);

//The next loop gets the name of the partition without blanks
  	while(buf[i] != '=')
	     {
	      if(buf[i] != ' ')
		 {
	          part_name[k].push_back(buf[i]);
		 }
	      i++;
	     }

//Ignore everything until a '-' appears
  	while(buf[i]!='-')
	     {
	      i++;
	     } 

//The end of the partition is written between '-' and ';'
	while(buf[i]!=';')
	     {

//Extraction of numerals
	      if((buf[i])>='0'&&(buf[i]<='9'))
		{
		 str[j]=buf[i];
		 j++;
		 i++;
		}	
	      else i++;
	     }

         partition_bound.push_back(se);
         partition_bound[k].stop = atoi(str);
	 if(k>0) partition_bound[k].start = partition_bound[k-1].stop+1;
	 else partition_bound[k].start = 1;
	 k++;
       }
  file.close();
}
 
/************************************************************************************************/

void get_taxa(char* filename, vector<string>& taxon_name)

{
  ifstream file(filename, ios::in);
  if(!file.good())
    {
     cerr << "Cannot open file " << filename << " ! Check if the filename "
          << "is written correctly.\n";
     exit(1);
    }  

  string line, name;

  while(!file.eof())
       {
        int c = file.get();

        if(c=='>')
	  {
	   getline(file,name);
           taxon_name.push_back(name);
  	  }
       }
//vector fixed_taxa contains the numbers of the fixed taxa
  vector<int> fixed_taxa;

  file.close();
}
	
/*************************************************************************************************/
void pres_abs_matrix(int start,int stop,char* filename, vector< vector <double> >& matrix, 
		     vector<double> pres_abs, int i, vector<string>& partition_seq,
		     vector<string>& taxon_name)
{
  ifstream file(filename, ios::in);
  if(!file.good())
    {
     cerr << "Cannot open file " << filename << " ! Check if the filename "
          << "is written correctly.\n";
     exit(1);
    }  

  int j=0;
  int k=0;
  int m=0;

  string line, name;

  while(!file.eof())
       {
        int c = file.get();

//Empty lines are ignored:
	if(c=='\n') continue;

//taxon_name only gets the name of the taxon, when the loop iterates 
//over the MSA for the first time (i=0). For each taxon another line 
//in the presence-absence matrix is generated
        if(c=='>' && i==0)
	  {
	   getline(file,name);
//	   taxon_name.push_back(name);
	   matrix.push_back(pres_abs);
  	  }

//After the first iteration over the MSA, lines with taxon names in it 
//are ignored
	else if(c=='>' && i!=0)
	file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

//All other lines only contain sequences. See also function check_fas 	
	else if(!file.eof())
	  {
	   file.putback(c);
	   getline(file,line);
	   partition_seq.push_back(line.substr(start-1,stop-(start-1)));

//Calculation of the gappiness
	   double gapcount = 0;
	   for(int x=0;x<(stop-(start-1));x++) 	    
	      {
	       if((partition_seq[k][x]=='-')||(partition_seq[k][x]=='X')
		  ||(partition_seq[k][x]=='?')) gapcount += 1;
	      }
	   double gappiness = gapcount/(double)(stop-(start-1)); 

	   if(gappiness<(1))
	     {
	      matrix[j][i]=1;
	      k++;
	     }

//Otherwise the sequence is dropped
	   else partition_seq.pop_back();
	
	   j++;
	   m++;
        }
     }
  file.close();
}
	
/*************************************************************************************************/

void get_barycent(int n,vector<string>& partition_seq, int (*scorematrix)[CHAR_MAX],
		    double(*barycent)[3])
{
//The following 4 loops determine all possible quartets. geometry_mapping calculates the value of 
//the corresponding sequence-quartet. The result is stored in barycent.  
 
  int x=0;

  for(int i=0; i<n; i++)
     {
      for(int j=1; j<n; j++)
         {
	  if(j<=i) continue;
	  for(int k=2; k<n; k++)
	     {
	      if(k<=j) continue;
	      for(int l=3; l<n; l++)
		 {
		  if(l<=k) continue;
		  else
		    {

//To avoid a pattern of Tree 1 and Tree 3 in the mapped quartet. The sequence 
//numbers have to get shuffled, otherwise they would be in ascending order
	             vector <int> shuf_vec(4);
 
	             shuf_vec[0] = i;
                     shuf_vec[1] = j;
                     shuf_vec[2] = k;
                     shuf_vec[3] = l;

                     random_shuffle(shuf_vec.begin(), shuf_vec.end());

                     geometry_mapping(partition_seq[shuf_vec[0]],partition_seq[shuf_vec[1]],
                                      partition_seq[shuf_vec[2]],partition_seq[shuf_vec[3]],
                                      scorematrix, barycent[x]);
		     x++;
		    }
		 }
	     }
         }
     }
 }

/************************************************************************************************/

void get_barycent_rand(int n,vector<string>& partition_seq, int (*scorematrix)[CHAR_MAX],
		       double (*barycent)[3])

{
//A vector all_quart which contains the 4 integer values of each quartet is generated.  

  vector <quart> all_quart;
  int x = 0;

//all_quart has to store (n choose 4) times 4 integer values for all possible quartets. The 
//following 4 loops fill all_quart with all possible quartets.  

  all_quart.resize(binom(n,4));

  for(int i=0; i<n; i++)
     {
      for(int j=1; j<n; j++)
 	 {
 	  if(j<=i) continue;
	  for(int k=2; k<n; k++)
	     {
	      if(k<=j) continue;
	      for(int l=3; l<n; l++)
	 	 {
		  if(l<=k) continue;
		  else
		    {
		     all_quart[x].a = i;
		     all_quart[x].b = j;
		     all_quart[x].c = k;
		     all_quart[x].d = l;
 		     x++;
		    }
		 }
             }
         }
     }

  x=0;

  long int qua_count = all_quart.size();	

//The loop draws 20000 random quartets from all_quart

  long int j = 0;
  for(int i = NUM_OF_QUA;i>0;i--)
     {
//time(NULL) will nearly always be the same because it is the current
//time in seconds. That is why j is added
      srand(time(NULL)+j);

//r is a value between 0 and 1.
      double r = (double)rand()/RAND_MAX;

//p is a value between 0 and (n choose 4)
      double p = r*qua_count;

//It is converted to a long int. j lies between 0 and (n choose 4). 
      j = static_cast<long int>(p);

      if(j != qua_count) j++;

//j now lies between 1 and (n choose 4). All (n choose 4) possible quartets
//can be drawn with the same probability. all_quart[j-1] is calculated now: 


//To avoid a pattern of Tree 1 and Tree 3 in the mapped quartet. The sequence 
//numbers have to get shuffled, otherwise they would be in ascending order
 
      vector <int> shuf_vec(4);
 
      shuf_vec[0] = all_quart[j-1].a;
      shuf_vec[1] = all_quart[j-1].b;
      shuf_vec[2] = all_quart[j-1].c;
      shuf_vec[3] = all_quart[j-1].d;

      random_shuffle(shuf_vec.begin(), shuf_vec.end());

      geometry_mapping(partition_seq[shuf_vec[0]],partition_seq[shuf_vec[1]],
                       partition_seq[shuf_vec[2]],partition_seq[shuf_vec[3]],
                       scorematrix, barycent[x]);
      x++;

//It must be verified that a quartet cannot be drawn twice. So the last 
//quartet from all_quart gets the place of the drawn one.

      all_quart[j-1].a = all_quart[qua_count-1].a;
      all_quart[j-1].b = all_quart[qua_count-1].b;
      all_quart[j-1].c = all_quart[qua_count-1].c;
      all_quart[j-1].d = all_quart[qua_count-1].d;
      all_quart.pop_back();
      qua_count--;
     }
}
/************************************************************************************************/

double evaluate_partition(double (*barycent)[3], long int num_of_vec, 
			 vector< vector <double> >& matrix, int y)
{

//evaluate_partition calculates the tree support as described in "Nieselt-Struwe, 
//v. Haeseler 2001"

  double tree_sup = 0;

  for(int i=0;i<num_of_vec;i++)
     {
      double min = 0;	
      if((barycent[i][0] <= barycent[i][1]) && (barycent[i][0] <= barycent[i][2]))
  	min = barycent[i][0];
      else if ((barycent[i][1] <= barycent[i][0]) && (barycent[i][1] <= barycent[i][2])) 
	min = barycent[i][1];
      else
	min = barycent[i][2];

      if (min <= (1.0/6.0)) tree_sup += 1;
     }

  double p = tree_sup / double(num_of_vec);
  int n = matrix.size();
  for(int i=0;i<n;i++)
    {
     if(matrix[i][y]==1)
     matrix[i][y] = (matrix[i][y])*p;
    }
  return p;
}

/************************************************************************************************/

void eval_taxa(vector< vector <double> >& matrix)
{ 
  int n = matrix.size()-2;
  int m = matrix[0].size()-2;
  double min = 1;
  int num;

//The 2 loops iterate over all entries of the taxon matrix. The entries for one taxon are summed up 
//and divided by the number of partitions. The last column of taxon is filled with these values. 
//The second last contains the taxon number. 
  for(int i=0;i<n;i++)
     {
      matrix[i][m+1]=0;
      for(int j=0;j<m;j++)
	 {
	  if(matrix[i][j]!=(-1)) matrix[i][m+1] += matrix[i][j];
	 }
      matrix[i][m+1] = (matrix[i][m+1])/m;
     }
}

/************************************************************************************************/

void eval_gene(vector< vector <double> >& matrix)
{ 

//for some detailed info see eval_taxa. This function works in the same way but a partition is 
//evaluated.

  int n = matrix.size()-2;
  int m = matrix[0].size()-2;
  double min = 1;
  int num;

  for(int j=0;j<m;j++)
     {
      matrix[n+1][j]=0;
      for(int i=0;i<n;i++)
 	 {
	  if(matrix[i][j]!=(-1)) matrix[n+1][j] += matrix[i][j];
	 }
      matrix[n+1][j] = (matrix[n+1][j])/n;
     }
}

/************************************************************************************************/

void sort_matrix(vector< vector <double> >& matrix, vector<int>& fixed_taxa)
{
  int n = matrix.size();
  int m = matrix[0].size();

  int num;
  int i,j,x=0;
  double max = 0;

  vector < vector <double> > matrix2(matrix);

//The sorted lines (taxa) are stored in a temporary matrix2.

//Sorting of lines:

//The while-loop simply finds with each iteration the line that was evaluated with the highest 
//value by function eval_taxa. The whole line is than stored in the temporary matrix2.The 
//entry in matrix is set to -1. In the next iteration the second highest is copied to matrix2 etc.   

  while(x<(n-2))
       {
	for(i=0;i<(n-2);i++)
	   {
	    if(matrix[i][m-1]>=max)
              {
	       max=matrix[i][m-1];
	       num=i;
              }
	   }
        for(j=0;j<m;j++)
	   {
	    matrix2[x][j] = matrix[num][j];
	   }
        matrix[num][m-1]=(-1);
        x++;
        max=0;
       } 

//The last line holds the evaluations for each taxon. The second last contains the number of 
//a taxon. These lines are also copied to matrix2
  for(j=0;j<m;j++)
     {
      matrix2[n-2][j] = matrix[n-2][j];
      matrix2[n-1][j] = matrix[n-1][j];
     }

//Sorting of columns:

//Sorting of columns is performed in the same way as sorting of lines. But now  the columns are 
//copied from matrix2 back to matrix . After that procedure the matrix is completely sorted. 

  x=0;

  while(x<m-2)
       {
        for(j=0;j<m-2;j++)
           {
            if(matrix2[n-1][j]>=max)
              {
               max=matrix2[n-1][j];
               num=j;
              }
           }

        for(i=0;i<n;i++)
           {
            matrix[i][x] = matrix2[i][num];
           }
        matrix2[n-1][num]=(-1);
        x++;
        max=0;
       }

  for(j=0;j<n;j++)
     {
      matrix[j][m-2] = matrix2[j][m-2];
      matrix[j][m-1] = matrix2[j][m-1];
     }
}

/************************************************************************************************/

int cluster_taxa(vector< vector <double> >& matrix, int lower_limit)
{
  int count = 0;
  int n, m;
  int add_tax = 0;
  vector <int> tax_of_cluster(1);
  vector < vector<int> > cluster(1,tax_of_cluster);

//matrix[0] is the first taxon of cluster[0] 
  cluster[0][0]=0;

//clustering of taxa:

//if there is only one partition, there is just one cluster
  if(matrix[0].size()-2 == 1) return 1;

//Loop over all taxa
  for(int i=1;i<matrix.size()-2;i++)
     {

//Loop over all cluster
      for(int j=0;j<cluster.size();j++)
	 {

//Loop over all taxa of a cluster
	  for(int k=0;k<cluster[j].size();k++)
	     {

//Loop over all partitions  
	      for(int l=0;l<matrix[0].size()-2;l++)
	         {

//count increases if a taxon has one sequence in common with another taxon
//of a cluster 
		  if((matrix[i][l] >= 0) && (matrix[cluster[j][k]][l] >= 0))
		    {
		     count++;
		     if(count == lower_limit)
		       {

//If count is equal to lower_limit, cluster[j] gets another taxon
		        cluster[j].push_back(i);
			count = 0;
			add_tax=1;

//By setting l to matrix[0].size()-2 and k to cluster[j].size(), the last two 
//loops finish, and it is checked if the taxon is also part of the next Cluster
			l=matrix[0].size()-2;
			k=cluster[j].size();
		       }
		    }
		 }
	      count = 0;
	     }
 	 }

//If the taxon does not have two sequences in common, it becomes a new cluster
      if(!add_tax)
	 {
	  cluster.push_back(tax_of_cluster);
	  cluster[cluster.size()-1][0]=i;
	 }
       else add_tax=0;
     }

   if(cluster.size() == 1) return 1;

//clustering of clusters

//the clustering of taxa can lead to a result with taxa in two or more 
//different clusters. Therefore the clusters have to be clustered
  int add_clus = 0;

  for(int i=0;i<cluster.size();i++)
     {

//Comparison of the taxa in different clusters
      n = cluster[i].size();
      for(int j=1;j<cluster.size() && j>i;j++)
	 {
	  m = cluster[j].size();
	  for(int k=0; k<cluster[i].size(); k++)
	     {
	      for(int l=0; l<cluster[j].size(); l++)
		 {
		  
//If two clusters have the same taxa, it is erased in one of them, 
//afterwards the clusters merge and the number of clusters decreases 
		   if(cluster[j][l] == cluster[i][k])
		     {
		      cluster[i].erase(cluster[i].begin()+k);
		      add_clus=1;
		     }
		  }
	     }

//Merge cluster and erase one of them
	   if(add_clus)
	     {
	      add_clus=0;
	      for(int l=0;l<m;l++)
		 {
		  cluster[i].push_back(cluster[j][l]);
		 }
	      cluster.erase(cluster.begin()+j);
	     }
	 }

//As long as the size of a cluster changes there might be new taxa in it 
//that lead to an inclusion of other clusters. Therefore the process is 
//repeated.
      if(cluster[i].size()!=n) i--;
      if(cluster.size()==1) return 1;
     }
  return cluster.size();
}

/************************************************************************************************/

int cluster_partitions(vector< vector <double> >& matrix, int lower_limit)
{
  int count = 0;
  int n, m;
  int add_par = 0;
  vector <int> par_of_cluster(1);
  vector < vector<int> > cluster(1,par_of_cluster);

//matrix[][0] is the first partition of cluster[0] 
  cluster[0][0]=0;

//Loop over all partitions
  for(int i=1;i<matrix[0].size()-2;i++)
     {

//Loop over all cluster
      for(int j=0;j<cluster.size();j++)
	 {

//Loop over all partitions of a cluster
	  for(int k=0;k<cluster[j].size();k++)
	     {

//Loop over all taxa  
	      for(int l=0;l<matrix.size()-2;l++)
	         {

//count increases if a partition has one sequence in common with another partition
//of a cluster 
		  if((matrix[l][i] >= 0) && (matrix[l][cluster[j][k]] >= 0))
		    {
		     count++;
		     if(count == lower_limit)
		       {

//If count is equal to lower_limit, cluster[j] gets another partition
		        cluster[j].push_back(i);
			count = 0;
			add_par=1;

//By setting l to matrix.size()-2 and k to cluster[j].size(), the last two 
//loops finish, and it is checked if the partition is also part of the next Cluster
			l=matrix.size()-2;
			k=cluster[j].size();
		       }
		    }
		 }
	      count = 0;
	     }
 	 }

//If the partition does not have lower_limit sequences in common, it becomes a new cluster
      if(!add_par)
	 {
	  cluster.push_back(par_of_cluster);
	  cluster[cluster.size()-1][0]=i;
	 }
       else add_par=0;
     }

  if(cluster.size() == 1) return 1;

//clustering of clusters

//the clustering of partitions can lead to a result with partitions in two or more 
//different clusters. Therefore the clusters have to be clustered
  int add_clus = 0;

  for(int i=0;i<cluster.size();i++)
     {

//Comparison of the partitions in different clusters
      n = cluster[i].size();
      for(int j=1;j<cluster.size() && j>i;j++)
	 {
	  m = cluster[j].size();
	  for(int k=0; k<cluster[i].size(); k++)
	     {
	      for(int l=0; l<cluster[j].size(); l++)
		 {
		  
//If two clusters have the same partition, it is erased in one of them, 
//afterwards the clusters merge and the number of clusters decreases 
		   if(cluster[j][l] == cluster[i][k])
		     {
		      cluster[i].erase(cluster[i].begin()+k);
		      add_clus=1;
		     }
		  }
	     }

//Merge cluster and erase one of them
	   if(add_clus)
	     {
	      add_clus=0;
	      for(int l=0;l<m;l++)
		 {
		  cluster[i].push_back(cluster[j][l]);
		 }
	      cluster.erase(cluster.begin()+j);
	     }
	 }

//As long as the size of a cluster changes there might be new partitions in it 
//that lead to an inclusion of other clusters. Therefore the process is 
//repeated.
      if(cluster[i].size()!=n) i--;
      if(cluster.size()==1) return 1;
     }
  return cluster.size();
}

/************************************************************************************************/


int reduce_matrix(vector< vector <double> >& matrix, vector<int>& fixed_taxa, vector<int>& fixed_genes,double tax_weight)
{
  unsigned int boolean = 1;
  int n = matrix.size();
  int m = matrix[0].size();
  int num_tax; 
  int num_seq;
  double min_tax;
  double min_seq;

//find the lowest evaluated taxon. There should be at least 3 taxa left. The last two lines 
//are filled with taxon numbers and the evaluation. Therefore matrix.size() has to be 
//larger than 5. If this is not the case the  minimal evaluation is set to FLT_MAX (maximal float 
//value), so the taxon cannot be thrown out. 

  if(n>5 && n>fixed_taxa.size()+2)
    {
     int x=0;

//look if the first taxa of the matrix are fixed 

     for(int i=0; i<fixed_taxa.size();i++)
	{
	 if(matrix[x][m-2] == fixed_taxa[i])
	   {
	    x++;

//start the loop again
	    i = -1;
	   }
	}

     min_tax = matrix[x][m-1];

     for(int i=x+1; i<n-2; i++)
        {
	 boolean = 1;
         if(matrix[i][m-1]<=min_tax)
           {

//Look if the current minimal taxon is a fixed taxon
 	    for(int y=0; y<fixed_taxa.size(); y++)
	       {
		if(matrix[i][m-2] == fixed_taxa[y]) boolean = 0;
	       }

//boolean was initiated with 1. If it is still 1 the taxon is accepted
	    if(boolean)
	      {
               min_tax = matrix[i][m-1];
               num_tax=i;
	      }
           }
        }

//If the first chosen taxon which is not fixed is the minimal one, than get
//its position in the vector
      if(min_tax == matrix[x][m-1]) num_tax = x;  
    }

//By setting min_tax to DBL_MAX (maximal double value) it is avoided that the last taxon is thrown out. 
  else min_tax = DBL_MAX;

//find the lowest evaluated sequence. Here at least 1 sequence has to be left.
  int q;

  if(m>3)
    {
     min_seq = matrix[n-1][0];
     for(q=1;q<m-2;q++)
	{


/***************************************************
added code: check if gene is included in fixed_genes.
	    If yes, it cannot be excluded by setting
	    evaluation to DBL_MAX. Thus, different 
	    handling of fixed_genes and fixed_taxa. 
****************************************************/

	 for(int ft=0;ft<fixed_genes.size();ft++)
		{
		if(fixed_genes[ft]==matrix[n-2][q]) 
			{
			matrix[n-1][q]=DBL_MAX;
			break;
			}
		}

/***************************************************/

	 if(matrix[n-1][q]<=min_seq)
	   {
	    min_seq  = matrix[n-1][q];
	    num_seq = q;
	   }
        }
     if(min_seq == matrix[n-1][0]) num_seq = 0;  
    }

//By setting min_seq to DBL_MAX it is avoided that the last sequence is thrown out. 
  else min_seq = DBL_MAX;

//If there are more than 3 taxa and 1 sequence left, the lowest evaluated taxon or sequence 
//is thrown out. tax_weight weights the minimal taxon.

  double tw = tax_weight*min_tax;

//Erasure of the minimal sequence:
  if((min_seq <= tw) && !((min_tax == DBL_MAX)&&(min_seq == DBL_MAX)))
    {
     for(int i=0;i<n;i++)
	{
	 matrix[i].erase(matrix[i].begin()+num_seq);
	}
     eval_taxa(matrix);
    }

//Erasure of the minimal taxon
  else if((min_seq > tw) && !((min_tax == DBL_MAX)&&(min_seq == DBL_MAX)))  
    {

//If min_tax is DBL_MAX and the taxa are weighted with a low value or zero,
//tax_weight*DBL_MAX gets below min_seq. The lowest sequence has to be erased
     if(min_tax==DBL_MAX) 
       {
        for(int i=0;i<n;i++)
           {
            matrix[i].erase(matrix[i].begin()+num_seq);
           }
       }
     else
       {
	matrix.erase(matrix.begin()+num_tax);
	eval_gene(matrix);
       }
    }  

n = matrix.size();
m = matrix[0].size();

//If a sequence is available for less than 4 taxa, the partition is dropped.
//If it is the only partition left, it will be retained.

/****************************************************
added code: if n>5, otherwise -> Segmentation fault
	    !!!!not checked, if it runs consistently
	    !!!!Be aware!!!!
****************************************************/

if(n>5)
{
  if(matrix[0].size()>3)
    {
     int count=0;
     for(int j=0;j<m-2;j++)
        {
         for(int i=0;i<n-2;i++)
            {
             if(matrix[i][j]>=0)
               count++;
             if(count>3) break;
            }
         if(count<4)
	   {
	    if(matrix[0].size()>3)
	      {
               for(int i=0;i<n;i++)
	          {
	           matrix[i].erase(matrix[i].begin()+j);
	          }    
               eval_taxa(matrix);
	      }
	   }
	 count = 0;
	}
    }
}


if ((min_tax == DBL_MAX)&&(min_seq == DBL_MAX)) return 0;
else return 1;
}

/************************************************************************************************/

double matrix_pot_info(vector< vector <double> >& matrix)
{
  double a=0;
  int m = matrix[0].size();
  int n = matrix.size()-2;		
  for(int i=0;i<n;i++)
     {
      a += matrix[i][m-1];
     }
  double x = a/n;
  return x;
}

/************************************************************************************************/

double optimality_f(double cms, double oms, double ap, double iw) 
{
  double f;
  double lam = cms/oms;
  f = 1-abs(lam - pow(ap,iw*(1-ap)));
  return f;
}

/************************************************************************************************/

double unw_optimality_f(double cms, double oms, double ap)
{
  double f;
  double lam = cms/oms;
  f = 1-abs(lam - ap);
  return f;
}

/************************************************************************************************/

void get_fastaout(char* filename, vector<start_end> partition_bound, 
		  vector<double> partitions, vector<double> taxa)
{
 ifstream file(filename, ios::in);
 if(!file.good())
   {
    cout << "Cannot open " << filename << "!\n";
    exit(1);
   }
  string out = "results/";
  out.append(filename);
  out.append("_reduced");
  char * cstr_out;
  cstr_out = new char [out.size()+1];
  strcpy (cstr_out, out.c_str());

  ofstream fileout (cstr_out,ios_base::out);
  if(!fileout.good())
    {
     cout << "Cannot open out.fas!\n";
     exit(1);
    }

  int cur_par;
  int c;
  int start, stop;
  int startout, stopout;

  string line, name, seq;

  for(int j=0;j<taxa.size();j++)
     {
      int cur_tax = int(taxa[j]);
      int i = 0;
	
//Ignore taxa until a taxon of the resulting vector 
//taxa appears
      while(i<=cur_tax)
           {
	    c = file.get();

	    if(c=='>') 
	      {
	       i++;
	      }
	   }
      file.putback(c);

//Parsing of the header and the whole sequence. It was checked before
//that a header is followed by a sequence
      getline(file,name);
      fileout << name << "\n";
      getline(file,line); 

//Extracting the partitions of the result:
      for(int k=0; k<partitions.size(); k++)
	 {
	  cur_par  = int(partitions[k]);
	  startout = partition_bound[cur_par].start;
	  stopout  = partition_bound[cur_par].stop;
	  seq = line.substr(startout-1,stopout-(startout-1));

	  fileout << seq; 
	 }

      fileout << "\n";
      file.seekg(0,ios::beg);

	c = file.get();
      file.putback(c);

 }
  file.close();
fileout.close();
}

/************************************************************************************************/

void get_charsetout(char* filename, vector<start_end> partition_bound, vector <double> partitions)
{
 ifstream file(filename, ios::in);
 if(!file.good())
   {
    cout << "Cannot open " << filename << "!\n";
    exit(1);
   }

  string out = "results/";
  out.append(filename);
  out.append("_reduced");
  char * cstr_out;
  cstr_out = new char [out.size()+1];
  strcpy (cstr_out, out.c_str());

  ofstream fileout (cstr_out,ios_base::out);
  if(!fileout.good())
    {
     cout << "Cannot open out.charset!\n";
     exit(1);
    }

  int end_of_part = 0;
  int start_of_part = 0;
  int part, part_count, c, part_len, start, stop, startout, stopout;

  for(int j=0;j<partitions.size();j++)
     {
      part = (int)partitions[j];
      startout = partition_bound[part].start;
      stopout  = partition_bound[part].stop;
      part_len = stopout - (startout-1);

      part_count = 0;
      c = 0;

//If it is not the first partition that has to be read:
      if(part>0)
	{
	 while(1)
	      {
	       c = file.get();

//Ignore everything until the first ';' and after that a '\n' appears
	       if(c==';')
		 {
		  while(c!='\n')
		       {
			c = file.get();
		       }
		  part_count++;
		 }

//If part_count is equal to part, the partition bounds are read
	       if(part_count == part)break;
	      }
	}	 

//Read everything until a '=' appears and write it to fileout
      while(c != '=')
	   {
	    c = file.get();
	    fileout << (char)c;
	   }

//Get to the beginning of the file. Remember that partitions
//are not sorted in vector partitions
      file.seekg(0L,ios::beg); 
      end_of_part += part_len;

//Get the partition bounds
      if(j==0)
        fileout << " 1 - " << end_of_part << " ;\n";
     
      if(j!=0)
	fileout << " " << start_of_part+1 << " - " << end_of_part << ";\n";

      start_of_part += part_len;
     }
  file.close();
  fileout.close();
}


/************************************************************************************************/

void opt_taxa_constraints (char** argv, int argc, vector<string> taxon_name, 
     vector<int>& fixed_taxa)
{
  string tax;

  for(int i=3; i<argc; i++)
     {
      if(argv[i][0] == '-')
	{
	 if(argv[i][1] == 'c')
	   {

//If there is no argumet after -c an error message appears
	    if(i+1 == argc)
	      {
	       cerr << "Missing Identification of fixed taxon! Enter the FASTA header "
	       	       "of the fixed taxon without '>' at the beginning when using "
		       "option -c.\n";
	       exit(1); 
	      }

//If there is something except another option starting with '-', it is put into a string
	    else
	      {
	       i++;
	       while((i+1 <= argc)&&(argv[i][0] != '-'))
		    {
		     tax.append(argv[i]);

//The taxons name may consist of more than one word. Words are seperated with a whitespace 
		     tax.append(" ");
		     i++;
		    }
	       i--;

//If there is no name of a fixed taxon an error message appears
	       if(tax.empty())
		 {
		  cerr << "Missing Identification of fixed taxon! Enter the FASTA header "
                  	  "of the fixed taxon without '>' at the beginning when using "
			  "option -c.\n";
		  exit(1);
		 }
	       else
		 {

//The last whitespace has to to be deleted
		  tax.erase(tax.end()-1);
		  int x=0;      

//Search for the name in vector of extracted names
		  for(int j=0;j<taxon_name.size();j++)
		     {
		      if(tax.compare(taxon_name[j])==0)
			{
			 x++;
			 fixed_taxa.push_back(j);
		        }
		     }

//If there is no taxon that matches to the argument, an error message appears
		  if(x==0)
		    {
		     cerr << "Taxon " << tax << " not found in multiple sequence "
		          << "alignment!\n";
                     exit(1);
		    }
		  tax.clear();
		 }
	      }
	   }
	}
     }
}  

/************************************************************************************************/

void opt_gene_constraints (char** argv, int argc, vector<string> part_name, 
     vector<int>& fixed_gene)
{
  string gen;

  for(int i=3; i<argc; i++)
     {
      if(argv[i][0] == '-')
	{
	 if(argv[i][1] == 'g')
	   {

//If there is no argumet after -g an error message appears
	    if(i+1 == argc)
	      {
	       cerr << "Missing Identification of fixed gene when using option -g.\n";
	       exit(1); 
	      }

//If there is something except another option starting with '-', it is put into a string
	    else
	      {
	       i++;
	       while((i+1 <= argc)&&(argv[i][0] != '-'))
		    {
		     gen.append(argv[i]);

//The taxons name may consist of more than one word. Words are seperated with a whitespace 
		     gen.append(" ");
		     i++;
		    }
	       i--;

//If there is no name of a fixed taxon an error message appears
	       if(gen.empty())
		 {
		  cerr << "Missing Identification of fixed gene when using option -g.\n";
		  exit(1);
		 }
	       else
		 {

//The last whitespace has to to be deleted
		  gen.erase(gen.end()-1);
		  int x=0;      

//Search for the name in vector of extracted names
		  for(int j=0;j<part_name.size();j++)
		     {
		      if(gen.compare(part_name[j])==0)
			{
			 x++;
			 fixed_gene.push_back(j);
		        }
		     }

//If there is no taxon that matches to the argument, an error message appears
		  if(x==0)
		    {
		     cerr << "Gene " << gen << " not found in multiple sequence "
		          << "charset file!\n";
                     exit(1);
		    }
		  gen.clear();
		 }
	      }
	   }
	}
     }
}

/************************************************************************************************/
void opt_reuse_matrix (char** argv, int argc, vector< vector <double> >& matrix, vector<double> pres_abs)
{

  string mat;
  double d;
  char buffer[10];

  for(int i=3; i<argc; i++)
     {
      if(argv[i][0] == '-')
	{
	 if(argv[i][1] == 'r')
	   {

//If there is no argumet after -r an error message appears
	    if(i+1 == argc)
	      {
	       cerr << "Missing Identification of matrix when using option -r.\n";
	       exit(1); 
	      }

//If there is something except another option starting with '-', it is put into a string
	    else
	      {
	       i++;
	       while((i+1 <= argc)&&(argv[i][0] != '-'))
		    {
		     mat.append(argv[i]);
		     i++;
		    }
	       i--;

//If there is no filename
	       if(mat.empty())
		 	{
		  	cerr << "Missing Identification of matrix when using option -r.\n";
		  	exit(1);
		 	}
	       else
		 	{

			char* filename = new char [mat.size()+1];
			strcpy (filename, mat.c_str());

		  	ifstream file(filename, ios::in);
		  	if(!file.good())
    		 	{
		  	cerr << "Cannot open file " << filename << " ! Check if the filename " << "is written correctly.\n";
    		  	exit(1);
    			}

        		char c				;
			int lc 		= 0		;
			int entry 	= 0		;
			int count	= 0		;
			
			while(c != '\n')
				{
				c=file.get();
				count++;
				}
		
			char* cstr  = new char [count];
			char* cstr2 = new char [50];
			
			stringstream sstr;
 			string str;
 			double f;

  				while(!file.eof())
       				{

				pres_abs.clear();	

				file.getline(cstr,count);

				if(file.eof()) continue;
	
				for(int i=0;i<count;i++)
					{

					if(cstr[i] != ' ' && cstr[i] != '\0' && cstr[i] != '\t')
						{
						str.append(1,cstr[i]);
						}
					else 
						{
						if(str.size()>0)
							{
// 							sstr<<str;
//							sstr>>f
							strcpy (cstr2, str.c_str());
							f = atof(cstr2);
							str.clear();
							pres_abs.push_back(f);
							}
						else continue;
						}
					}
				pres_abs.erase(pres_abs.begin());
				pres_abs.push_back(-1);
				pres_abs.push_back(-1);
				matrix.push_back(pres_abs);
				}
				
				file.close();
			}
		}
	}
}
}
}


/************************************************************************************************/

double opt_taxa_weighting (char** argv, int argc)
{
  string str;

  double weight = 1.0;	

  for(int i=3; i<argc; i++)
     {
      if(argv[i][0] == '-')
	{
	 if(argv[i][1] == 't')
	   {

//If there is no argumet after -t an error message appears
	    if(i+1 == argc)
	      {
	       cerr << "Missing weighting of taxa after option -t.\n";
	       exit(1); 
	      }

//If there is something except another option starting with '-', it is put into a string
	    else
	      {
	       i++;
	       while((i+1 <= argc)&&(argv[i][0] != '-'))
		    {
		     str.append(argv[i]);

//Appending " " is a relict of the taxon constraint function which was template of this one 
		     str.append(" ");
		     i++;
		    }
	       i--;

//If there is no value an error message appears
	       if(str.empty())
		 {
		  cerr << "Missing weighting of taxa after option -t.\n";
		  exit(1);
		 }
	       else
		 {

//The last whitespace has to to be deleted
		  str.erase(str.end()-1);
		  int x=0;      

		  char * cstr;

  		  cstr = new char [str.size()+1];
  		  strcpy (cstr, str.c_str());

		  weight = strtod (cstr,NULL);
		 }
	      }
	   }
	}
     }
 return weight;
}  
/************************************************************************************************/

double opt_pics (char** argv, int argc)
{
  char a;
  for(int i=3; i<argc; i++)
     {
      if(argv[i][0] == '-')
        {
         if(argv[i][1] == 's') return 1;
	}
     }
  return 0;
}

/************************************************************************************************/

double opt_matrix (char** argv, int argc)
{
  char a;
  for(int i=3; i<argc; i++)
     {
      if(argv[i][0] == '-')
        {
         if(argv[i][1] == 'm') return 1;
	}
     }
  return 0;
}

/************************************************************************************************

double opt_info_weighting (char** argv, int argc)
{
  char a;
  for(int i=3; i<argc; i++)
     {
      if(argv[i][0] == '-')
        {
         if(argv[i][1] == 'd') return 0;
        }
     }
  return 3;
}

************************************************************************************************/

double opt_info_weighting (char** argv, int argc)
{
  string str;

//default setting = 3

  double weight = 3.0;
	
  for(int i=3; i<argc; i++)
     {
      if(argv[i][0] == '-')
	{
	 if(argv[i][1] == 'd')
	   {

//If there is no argumet after -d an error message appears
	    if(i+1 == argc)
	      {
	       cerr << "Missing weighting of information content after option -d.\n";
	       exit(1); 
	      }

//If there is something except another option starting with '-', it is put into a string
	    else
	      {
	       i++;
	       while((i+1 <= argc)&&(argv[i][0] != '-'))
		    {
		     str.append(argv[i]);

//Appending " " is a relict of the taxon constraint function which was template of this one 
		     str.append(" ");
		     i++;
		    }
	       i--;

//If there is no value an error message appears
	       if(str.empty())
		 {
		  cerr << "Missing weighting of information content after option -d.\n";
		  exit(1);
		 }
	       else
		 {

//The last whitespace has to to be deleted
		  str.erase(str.end()-1);
		  int x=0;      

		  char * cstr;

  		  cstr = new char [str.size()+1];
  		  strcpy (cstr, str.c_str());
		  weight = strtod (cstr,NULL);

		 }
	      }
	   }
	}
     }
  return weight;
}  
/************************************************************************************************/

void opt_help (char** argv, int argc)
{
  char a;
  for(int i=0; i<argc; i++)
     {
      if(argv[i][0] == '-')
        {
         if(argv[i][1] == 'h')
           {
cout << "\n\n"
     << "         HOW TO USE MARE (MATRIX REDUCTION)\n"
     << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-\n\n"

     << "RUNNING MARE\n"
     << "------------\n\n"

     << "MARE requires two files as command line arguments. The\n" 
     << "first one determines start and end of each partition\n" 
     << "(gene or protein) in the concatenated multiple sequence\n" 
     << "alignment. It must have the following format:\n\n"

     << "charset YourFirstPartition = 1 - 100 ;\n"
     << "charset YourSecondPartition = 101 - 255 ;\n"
     << "etc\n\n"

     << "The second file is the concatenated multiple sequence\n" 
     << "alignment which has to be a multiple sequence FASTA\n" 
     << "file.\n\n"

     << "To run MARE with default settings enter\n\n"
     << "./MARE file1 file2\n\n"
     << "OPTIONS\n"
     << "-------\n\n" 
     << "-c [taxon_1] -c [taxon_2] -c [taxon_3] etc\n\n"
     << "It is possible to put constraints on one or more taxa.\n" 
     << "These taxa cannot be thrown out while reduction. To use\n" 
     << "this option enter the FASTA header of the taxon without\n"
     << "'>' at the beginning.\n\n"
     << "-g [gene_1] -g [gene_2] -g [gene_3] etc\n\n"
     << "It is possible to put constraints on one or more genes.\n" 
     << "These genes cannot be thrown out while reduction. To use\n" 
     << "this option enter the gene name of the charset file.\n\n"
     << "-d [alpha]\n\n"
     << "By choosing this option the user can change the default\n" 
     << "weighting of information content in MARE (default=3.0). Generally higher\n" 
     << "weightings lead to smaller subsets of taxa and partitions with higher infor-\n"
     << "mation content.\n\n"
     << "-m\n\n"
     << "This option outputs the reduced and unreduced matrix as *.svg.\n\n"
     << "-s\n\n"
     << "This option outputs the simplex graphs as *.svg.\n"
     << "WARNING: A *.svg of a simplex graph requires more than 3 MB\n"
     << "         storage space. So 1000 partitions will take more\n"
     << "         than 3 GB.\n\n"
     << "-t [taxon weight] \n\n"
     << "Higher taxon weights cause stronger weighting of taxa while re-\n"
     << "duction. Thus more taxa remain in the resulting alignment (default=1).\n"
     << "Only values > 1 are allowed.\n\n"
     << "-r [filename] \n\n"
     << "By choosing this option a precalculated data matrix as output by MARE (matrix.txt)\n"
     << "can be used instead of calculating information contents again. It is required that\n" 
     << "taxa and genes in the matrix have the same order as in the msa and charset file.\n\n"
     << "Here is an example of how to use the options in MARE:\n\n"
     << "./MARE file1 file2 -c Homo_sapiens -c Aplysia_californica -t\n\n";
     
     exit(1);
	  }
	}
     }
}

/************************************************************************************************/

void draw_quartetmap (int y, double (*barycent)[3], int num_of_vec, string part_name, 
		      double rel_info)
{
  string out = part_name;
  size_t pos = part_name.size();
  out = out.insert(pos,"_map.svg");  
  out = out.insert(0,"results/");
  
  char * cstr;
  cstr = new char [out.size()+1];
  strcpy (cstr, out.c_str());

  ofstream fileout (cstr,ios_base::out);
  if(!fileout.good())
    {
     cerr << "Cannot open " << out << "!\n";
     exit(1);
    }

  string init_line = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" ;
  string gen_line  = "<!-- created by matrix_reduction.pl -->";

  double width     = 100.0;
  double height    = 50.0*sqrt(3.0);	
  double x_P2      = width/2.0;

  fileout << init_line << "\n"
          << gen_line << "\n\n";

  fileout << "<svg\n"
	  << "xmlns:svg=\"http://www.w3.org/2000/svg\"\n"
          << "xmlns=\"http://www.w3.org/2000/svg\"\n"
          << "version=\"1.0\"\n"
          << "width=" << "\"" << width << "\"\n"
          << "height=" << "\"" << height << "\"\n"
          << "id=\"svg2\">\n\n"

  	  << "<defs\n"
	  << "id=\"defs4\" />\n\n"

  	  << "<text x=\"5\" y=\"10\" font-size=\"3px\">     gene: " << part_name << "</text>\n"	
  	  << "<text x=\"5\" y=\"13\" font-size=\"3px\">rel. info: " << rel_info    << "</text>\n\n"
  
  	  << "<g\n"
	  << "transform=\"matrix(0.8,0,0,0.8,13,11)\"\n"
	  << "id=\"g0000\">\n\n";

  double x_center = 50.0;
  double y_center = (100.0/3.0)*(sqrt(3));

//unit vectors
 
  double x_e1 = -50.0;
  double x_e2 =  50.0;
  double x_e3 =   0.0; 
  double y_e1 = (50.0 * sqrt(3))/3.0;
  double y_e2 = (50.0 * sqrt(3))/3.0;
  double y_e3 = -2.0 * ((50 * sqrt(3))/3.0);

//outer triangle

  double T1[2];
  double T2[2];
  double T3[2];

  T1[0] = x_center + ((1) * x_e1) + ((0/6) * x_e2) + ((0/6) * x_e3);
  T1[1] = y_center + ((3/3) * y_e1) + ((0/6) * y_e2) + ((0/6) * y_e3);

  T2[0] = x_center + ((0/6) * x_e1) + ((3/3) * x_e2) + ((0/6) * x_e3);
  T2[1] = y_center + ((0/6) * y_e1) + ((3/3) * y_e2) + ((0/6) * y_e3);

  T3[0] = x_center + ((0/6) * x_e1) + ((0/6) * x_e2) + ((3/3) * x_e3);
  T3[1] = y_center + ((0/6) * y_e1) + ((0/6) * y_e2) + ((3/3) * y_e3);


fileout << "<path\n"
        << "d=\"M " << T1[0] << "," << T1[1] << " L " << T2[0]<<","<<T2[1]<<" L "<<T3[0]  
        << ","<<T3[1]<<" L "<<T1[0]<<","<<T1[1]<< " z\"\n"
        << "style=\"fill:none;fill-opacity:0;stroke:black;stroke-width:0.5;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
        << "id=\"path18\" />\n\n";
  
//For a given vector the following loop starts

  int id = 0;

  for(int i=0;i<num_of_vec;i++)
     {
      string colour = "#8080FF";	
      srand(time(NULL)+id);
//r is a value between 0 and 1.
      double r = (double)rand()/RAND_MAX;
//p is a value between 0 and (100000)
      double p = r*100000;
//It is converted to a int. j lies between 0 and (100000). 
      id = static_cast<int>(p);
	
      double x_c = x_center + (barycent[i][0] * x_e1) + (barycent[i][1] * x_e2) + 
	           (barycent[i][2] * x_e3) ;
      double y_c = y_center + (barycent[i][0] * y_e1) + (barycent[i][1] * y_e2) + 
 	           (barycent[i][2] * y_e3) ;

      double max;
      double min;

//Topology 1

      if(barycent[i][1]>=barycent[i][2]) 
	{ 
	 max = barycent[i][1];
	 min = barycent[i][2];
	}
      else 
	{
	 max = barycent[i][2];
	 min = barycent[i][1];
	}
      if(((barycent[i][0]-max) > 1.0/2.0) && (3.0/4.0 < (barycent[i][0]+min)))
	{
	 colour = "#000080"; 
        }
//Topology 2

      if(barycent[i][0]>=barycent[i][2]) 
	{ 
	 max = barycent[i][0];
	 min = barycent[i][2];
	}
      else 
	{
	 max = barycent[i][2];
	 min = barycent[i][0];
	}
      if(((barycent[i][1]-max)>1.0/2.0) && (3.0/4.0 < (barycent[i][1]+min)))
	{
	 colour = "#000080"; 
	}
//Topology 3

      if(barycent[i][0]>=barycent[i][1]) 
	{ 
	 max = barycent[i][0];
	 min = barycent[i][1];
	}
      else 
	{
	 max = barycent[i][1];
	 min = barycent[i][0];
	}

      if(((barycent[i][2]-max)>1.0/2.0) && (3.0/4.0 < (barycent[i][2]+min)))
	{
	 colour = "#000080"; 
	}

      if((barycent[i][0]<=barycent[i][1]) && (barycent[i][0]<=barycent[i][2]))
	 min = barycent[i][0];
      else if((barycent[i][1]<=barycent[i][0]) && (barycent[i][1]<=barycent[i][2]))
	     min = barycent[i][1];
      else min = barycent[i][2];

      if(min>1.0/6.0) colour = "#AAAAFF";
  
      fileout << "<circle\n"
	      << "cx=\"" << x_c <<"\" cy=\"" << y_c << "\" r=\"0.5\"\n"
     	      << "style=\"fill:" << colour << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:0\"\n"
              << "id=\"circle" << id << "\" />\n";
     }


//inner triangle points: 	t1_(2/3,1/6,1/6)
//			t1_(1/6,2/3,1/6)
//			t1_(1/6,1/6,2/3)
//			a(3/4,1/4,0)
//			b(1/4,3/4,0)
//			c(0,3/4,1/4)
//			d(0,1/4,3/4)
//			e(1/4,0,3/4)
//			f(3/4,0,1/4)

//inner triangle 

  double t1[2];
  double t2[2];
  double t3[2];

  t1[0] = x_center + ((2.0/3.0) * x_e1) + ((1.0/6.0) * x_e2) + ((1.0/6.0) * x_e3);
  t1[1] = y_center + ((2.0/3.0) * y_e1) + ((1.0/6.0) * y_e2) + ((1.0/6.0) * y_e3);

  t2[0] = x_center + ((1.0/6.0) * x_e1) + ((2.0/3.0) * x_e2) + ((1.0/6.0) * x_e3);
  t2[1] = y_center + ((1.0/6.0) * y_e1) + ((2.0/3.0) * y_e2) + ((1.0/6.0) * y_e3);

  t3[0] = x_center + ((1.0/6.0) * x_e1) + ((1.0/6.0) * x_e2) + ((2.0/3.0) * x_e3);
  t3[1] = y_center + ((1.0/6.0) * y_e1) + ((1.0/6.0) * y_e2) + ((2.0/3.0) * y_e3);

//Voronoi cell borders of tree attractors

  double a[2];
  double b[2];
  double c[2];
  double d[2];
  double e[2];
  double f[2];

  a[0] = x_center + ((3.0/4.0) * x_e1) + ((1.0/4.0) * x_e2) + ((0.0/4.0) * x_e3);
  a[1] = y_center + ((3.0/4.0) * y_e1) + ((1.0/4.0) * y_e2) + ((0.0/4.0) * y_e3);

  b[0] = x_center + ((1.0/4.0) * x_e1) + ((3.0/4.0) * x_e2) + ((0.0/4.0) * x_e3);
  b[1] = y_center + ((1.0/4.0) * y_e1) + ((3.0/4.0) * y_e2) + ((0.0/4.0) * y_e3);

  c[0] = x_center + ((0.0/4.0) * x_e1) + ((3.0/4.0) * x_e2) + ((1.0/4.0) * x_e3);
  c[1] = y_center + ((0.0/4.0) * y_e1) + ((3.0/4.0) * y_e2) + ((1.0/4.0) * y_e3);

  d[0] = x_center + ((0.0/4.0) * x_e1) + ((1.0/4.0) * x_e2) + ((3.0/4.0) * x_e3);
  d[1] = y_center + ((0.0/4.0) * y_e1) + ((1.0/4.0) * y_e2) + ((3.0/4.0) * y_e3);

  e[0] = x_center + ((1.0/4.0) * x_e1) + ((0.0/4.0) * x_e2) + ((3.0/4.0) * x_e3);
  e[1] = y_center + ((1.0/4.0) * y_e1) + ((0.0/4.0) * y_e2) + ((3.0/4.0) * y_e3);

  f[0] = x_center + ((3.0/4.0) * x_e1) + ((0.0/4.0) * x_e2) + ((1.0/4.0) * x_e3);
  f[1] = y_center + ((3.0/4.0) * y_e1) + ((0.0/4.0) * y_e2) + ((1.0/4.0) * y_e3);

  fileout << "<path\n"
	  << "d=\"M " << t1[0] << "," << t1[1] << " L " << t2[0] << "," << t2[1] << " L "
	  << t3[0] << "," << t3[1] << " L " << t1[0] << "," << t1[1] << " z\"\n"
          << "style=\"fill:none;fill-opacity:0;stroke:black;stroke-width:0.2;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"path18a\" />\n";

  fileout << "<path\n"
	  << "d=\"M " << a[0] << "," << a[1] << " L " << t1[0] << "," << t1[1] << " L "
	  << f[0] << "," << f[1] << "\"\n"
          << "style=\"fill:none;fill-opacity:0;stroke:black;stroke-width:0.2;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"path18b\" />\n";

  fileout << "<path\n"
	  << "d=\"M " << b[0] << "," << b[1] << " L " << t2[0] << "," << t2[1] << " L "
	  << c[0] << "," << c[1] << "\"\n"
          << "style=\"fill:none;fill-opacity:0;stroke:black;stroke-width:0.2;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"path18c\" />\n";

  fileout << "<path\n"
	  << "d=\"M " << d[0] << "," << d[1] << " L " << t3[0] << "," << t3[1] << " L "
	  << e[0] << "," << e[1] << "\"\n"
          << "style=\"fill:none;fill-opacity:0;stroke:black;stroke-width:0.2;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"path18d\" />\n\n";

  fileout << "</g>\n"
	  << "</svg>\n";
}

/************************************************************************************************/

void draw_matrix(vector< vector <double> >& matrix, string filename, 
		 vector<string> taxon_name, vector<string> part_name)
{
  int nrows = matrix.size()-2;
  nrows *= 10;
 
  int ncolumns = matrix[0].size()-2;
  ncolumns *= 10;

  string init_line = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>";
  string gen_line  = "<!-- created by matrix_reduction.pl -->";

  char * cstr;
  cstr = new char [filename.size()+1];
  strcpy (cstr, filename.c_str());

  ofstream fileout (cstr,ios_base::out);
  if(!fileout.good())
    {
     cout << "Cannot open " << filename << "!\n";
     exit(1);
    }

  int width = ncolumns + 2;
  int height= nrows    + 2;

  fileout << init_line << "\n"
          << gen_line << "\n";

  fileout << "<svg\n"
	  << "xmlns:svg=\"http://www.w3.org/2000/svg\"\n"
          << "xmlns=\"http://www.w3.org/2000/svg\"\n"
          << "version=\"1.0\"\n"
          << "width=" << "\"" << width << "\"\n"
          << "height=" << "\"" << height << "\"\n"
          << "id=\"svg2\">\n\n"

  	  << "<defs\n"
	  << "id=\"defs4\" />\n\n"

	  << "<rect\n"
     	  << "width=\"" << width << "\"\n"
     	  << "height=\"" << height << "\"\n"
     	  << "x=\"0\"\n"
     	  << "y=\"0\"\n"
     	  << "style=\"opacity:1;fill:white;fill-opacity:0;stroke:black;stroke-width:1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
     	  << "id=\"rect001\" />\n\n";

  int y=1;

  vector <string> c(12); //km: changed from c(11) to c(12), because a class more included
		
		c[0]  = "white"; //km: gene not present 
		c[1]  = "#00112b"; //km: tree-likeness: >0.9-1
		c[2]  = "#025"; //km: tree-likeness: >0.8-0.9
		c[3]  = "#003380"; //km: tree-likeness: >0.7-0.8    
		c[4]  = "#04a"; //km: tree-likeness: >0.6-0.7
		c[5]  = "#0055d4"; //km: tree-likeness: >0.5-0.6
		c[6]  = "#06f"; //km: tree-likeness: >0.4-0.5
		c[7]  = "#2a7fff"; //km: tree-likeness: >0.3-0.4
		c[8]  = "#59f"; //km: tree-likeness: >0.2-0.3
		c[9]  = "#80b3ff"; //km: tree-likeness: >0.1-0.2
		c[10] = "#acf"; //km: tree-likeness: >0-0.1 
		c[11] = "#d40000"; //km: tree-likeness: = 0: color: red; c[11] new added!

  int id = 0;   

  int par = matrix[0].size()-2;
  int tax = matrix.size()-2;

  for(int i=0; i<tax; i++)
     {	  
      int x=1;
      for(int j=0; j<par; j++)
	 {
	  double rect = matrix[i][j];
	  	  
	  srand(time(NULL)+id);
//r is a value between 0 and 1.
          double r = (double)rand()/RAND_MAX;
//p is a value between 0 and (100000)
          double p = r*100000;
//It is converted to a int. j lies between 0 and (100000). 
          id = static_cast<int>(p);

          if(rect == -1) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[0] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	    //km: one class added / defined: if tree-likeness is 0 -> color it red!!!
	    	  else if(rect == 0) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[11] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	    //km: end addition of new class
	    	  else if( rect > 0 && rect <= 0.1 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[10] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.1 && rect <= 0.2 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[9] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.2 && rect <= 0.3 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[8] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.3 && rect <= 0.4 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[7] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.4 && rect <= 0.5 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[6] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.5 && rect <= 0.6 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[5] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.6 && rect <= 0.7 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[4] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.7 && rect <= 0.8 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[3] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.8 && rect <= 0.9 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[2] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else if( rect > 0.9 && rect <= 1 ) 
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[1] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  else
	    {
  fileout << "<rect\n"
          << "width=\"10\"\n"
	  << "height=\"10\"\n"
	  << "x=\"" << x << "\"\n"
	  << "y=\"" << y << "\"\n"
	  << "style=\"fill:" << c[0] << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\"\n"
	  << "id=\"rect" << id << "\" />\n";
	    }
	  x += 10;
	 }
      y += 10;
     }   

  y = 10;

//write taxa to matrix

  fileout << "<text\n"
	  << "x=\"-5\"\n"
	  << "y=\"10\"\n"
	  << "style=\"font-size:8px;font-style:normal;font-weight:normal;text-align:end;text-anchor:end;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans\"\n"
	  << "id=\"text1892\"\n"
	  << "xml:space=\"preserve\">\n";


  for(int i=0; i<tax; i++)
     {
      srand(time(NULL)+id);
//r is a value between 0 and 1.
      double r = (double)rand()/RAND_MAX;
//p is a value between 0 and (100000)
      double p = r*100000;
//It is converted to a int. j lies between 0 and (100000). 
      int id = static_cast<int>(p);
	
  fileout << "<tspan\n"  
	  << "x=\"-5\"\n"       
	  << "y=\"" << y << "\"\n"
	  << "style=\"text-align:end;text-anchor:end\"\n"
	  << "id=\"tspan" << id << "\">" << taxon_name[int(matrix[i][par])] << "</tspan>";
   y += 10;
  }

  fileout << "</text>\n\n";

//printing gene names, with shift removes first entry of genes which says just taxa

  int x = 9;

  for(int i=0; i<par; i++)
     {
      srand(time(NULL)+id);
//r is a value between 0 and 1.
      double r = (double)rand()/RAND_MAX;
//p is a value between 0 and (100000)
      double p = r*100000;
//It is converted to a int. j lies between 0 and (100000). 
      int id = static_cast<int>(p);
	
  fileout << "<text\n"
	  << "x=\"" << x << "\"\n"       
	  << "y=\"-5\"\n"
     	  << "transform=\"rotate(-90," << x << ",-5)\"\n"     	
          << "style=\"font-size:8px;font-style:normal;font-weight:normal;fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:Bitstream Vera Sans\"\n"
	  << "id=\"text8747\"\n"
	  << "xml:space=\"preserve\"><tspan\n"
          << "x=\"" << x << "\"\n"
	  << "y=\"-5\"\n"
	  << "id=\"tspan" << id << "\">" << part_name[int(matrix[tax][i])] << "</tspan></text>\n";
	 
      x += 10;
     }

  fileout << "</svg>\n";
}

