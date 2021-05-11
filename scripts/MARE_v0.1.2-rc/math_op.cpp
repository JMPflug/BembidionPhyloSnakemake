#include <iostream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "math_op.h"

using namespace std;

/*************************************************************************************************/

long int binom(int n,int k)
{
  if(n<k)
    {
     cerr << "Error: n must be bigger or equal k, when calculating 'n choose k'!";
     exit(1);
    }  
  int i,j;
  long int B[n+1][k+1];
  int min;
  for(i=0; i<=n;i++)
     {
      if (i<=k) 
      min = i;
      else min = k;
      for(j=0;j<=min;j++)
         {
	   if (j==0 || j==i) B[i][j]=1;
          else
	     B[i][j] = B[i-1][j-1]+B[i-1][j];
         }
     }
  return B[n][k];
}
/*************************************************************************************************/

/*************************************************************************************************/

void geometry_mapping(string& a, string& b, string& c, string&d, int (*s)[CHAR_MAX],
		      double *barycent)
 {

  int len = a.size();
  int i = 0;
  int count = 0;

  double sum_gab_cd = 0;
  double sum_gac_bd = 0;
  double sum_gad_bc = 0;

  for(int i=0;i<len;i++)
     {

//Gaps are not scored
      if(a[i] == '-' || b[i] == '-' || c[i] == '-' || d[i] == '-')continue;
	
      if(s[a[i]][b[i]] == INT_MAX ||
	 s[a[i]][c[i]] == INT_MAX ||
	 s[a[i]][d[i]] == INT_MAX ||
	 s[c[i]][d[i]] == INT_MAX ||
	 s[b[i]][d[i]] == INT_MAX ||
	 s[b[i]][c[i]] == INT_MAX)continue;
	
	 int ab = 11-s[a[i]][b[i]];
	 int ac = 11-s[a[i]][c[i]];
	 int ad = 11-s[a[i]][d[i]];
	 int cd = 11-s[c[i]][d[i]];
	 int bd = 11-s[b[i]][d[i]];
	 int bc = 11-s[b[i]][c[i]];

	 int x = ab + cd;
	 int y = ac + bd;
	 int z = ad + bc;

	 int max;

	 if ((x>=y)&&(x>=z)) max = x;
	 else if ((y>=x)&&(y>=z)) max = y;
	 else max = z;

	 double gab_cd = 0.5*(max-ab-cd);
	 double gac_bd = 0.5*(max-ac-bd);
	 double gad_bc = 0.5*(max-ad-bc);

	 sum_gab_cd += gab_cd;
	 sum_gac_bd += gac_bd;
	 sum_gad_bc += gad_bc;

	count++;
     }
   if((sum_gab_cd == 0)&&(sum_gac_bd == 0)&&(sum_gad_bc==0)) 
     {
      barycent[0] = 1.0/3.0;
      barycent[1] = 1.0/3.0;
      barycent[2] = 1.0/3.0;
     }
   else if (count<50)
     {
      barycent[0] = 1.0/3.0;
      barycent[1] = 1.0/3.0;
      barycent[2] = 1.0/3.0;
     }
   else
     {
      barycent[0] = sum_gab_cd / (sum_gab_cd + sum_gac_bd + sum_gad_bc); 
      barycent[1] = sum_gac_bd / (sum_gab_cd + sum_gac_bd + sum_gad_bc);
      barycent[2] = sum_gad_bc / (sum_gab_cd + sum_gac_bd + sum_gad_bc);
     } 

}

/*************************************************************************************************/

/*************************************************************************************************/

double euclid_dist(double* veca, double* vecb, int vecsize)
{
  double dist = 0;
  for(int i=0; i<vecsize; i++)
     {
      dist += pow((veca[i] - vecb[i]),2);
     }    
  return sqrt(dist);
}

