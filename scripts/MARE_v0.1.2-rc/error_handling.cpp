#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <locale>
#include <cstring>

using namespace std;

/*************************************************************************************************/

int check_linefeed(char* filename)
{
  ifstream file(filename, ios::in);
  if(!file.good())
    {
     cerr << "Cannot open file " << filename << " ! Check if the filename " 
          << "is written correctly.\n";
     exit(1);
    }

  int a;
  int c = file.get();

  while(!file.eof())
       {
        if(c=='\r') 
	  {
	   a = file.get();

//Check if the file has Windows format, otherwise the file has 
//MAC format
	   if(a == '\n') return 3;
	   else return 2;
	  }

//If there is no '\r' in the file, but '\n' it should be of Unix 
//format
	else if((c=='\n') || (file.eof())) return 1;
	c = file.get();
       }
  file.close();
  return 1;
}

/*************************************************************************************************/


/*************************************************************************************************/

int check_fas(char* filename)
{
  ifstream file(filename, ios::in);
  if(!file.good())
    {
     cerr << "Cannot open file " << filename << " ! Check if the\n"
          << "filename is written correctly.\n";
     exit(1);
    }

  int a;
  int c = 0;
  int i = 0;
  int linecount = 0;
  int seqcount = 0;
  int alicheck;
  
  string sequence;

  while(!file.eof())
       {
	sequence.empty();
	a = file.get();

//Empty lines are ignored
	if(a=='\n' || file.eof()) linecount++;

//If a header appears there has to be a sequence in the following line.
	else if(a=='>')
	  {
	   seqcount++;
	   linecount++;
	   file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	   linecount++;
	   getline(file,sequence);

	   if (file.eof())
		{
              	cerr << "Error in line " << linecount << " of " << filename
		     << ": Linefeed is missing at end of sequence!\n";
		exit(1);
		};

//All sequence lengths are compared to the length of the first sequence 
	   if(i==0) 
	     {
	      i++;
	      alicheck = sequence.size();
	     }
	   else if(alicheck != sequence.size())
             {
              cerr << "Error in line " << linecount << " of " << filename
		   << ": Sequences are not aligned.\n";
              exit(1);
	     }
	  }
	else

//Everything besides header and the sequence in the following line is not part of a FASTA file
          {
           cerr << "Error in line " << linecount << " of " << filename
		   << ": Check format of your FASTA file.\n";
           exit(1);
	  }
       }

//It takes at least 5 sequences to run MARE.
  if(seqcount<=4)
    {
     cerr << "Error in " << filename << ": Alignment cannot be reduced."
	  << " MARE requires more than 4 sequences.\n";
     exit(1);
    }
  file.close();
  return alicheck;
}

/*************************************************************************************************/

int check_part(char* filename)
{
  ifstream file(filename, ios::in);
  if(!file.good())
    {
     cerr << "Cannot open file " << filename << " ! Check if the filename "
          << "is written correctly.\n";
     exit(1);
    }

  string buf;
  string copy_buf;
  char str[10];
  int linecount = 0;
  int name_flag = 0;
  int cont_flag = 0;
  int part_end;

  while(!file.eof())
       {
        buf.empty();
        memset(str,'\0',10);
        copy_buf.clear();

//If a line is empty, it is ignored
        getline(file,buf);
	linecount++;
        if(buf.empty()) continue;

//If there is any content in the file, a content flag is set
        else cont_flag = 1;

//Translate buf to a string with only uppercase letters
        locale loc;
        for (size_t l=0; l<buf.length(); ++l)
            {
             copy_buf.push_back(toupper(buf[l],loc));
            }

//Find the beginning of substring CHARSET
        string substr ("CHARSET");
        size_t pos = copy_buf.find(substr);

//Check if the word charset appears in the file
	if(pos==string::npos) 
	  {
           cerr << "Incorrect format of " << filename << "! Word 'charset' is missing "
                << "in line " << linecount <<".\nFor further help type ./MARE -h.\n";
	   exit(1);
	  }

//i marks the position after CHARSET in buf.
        int i = int(pos)+substr.size();

//Check if '=' and a partition identifier appears in the line
        while(buf[i] != '=' )
             {
              if(buf[i] != ' ') name_flag = 1;
	      i++;
	      if(i==buf.size())
		{
		 cerr << "Incorrect format of " << filename << "! '=' is missing in line " 
		      <<  linecount << ".\nFor further help type ./MARE -h.\n";
		 exit(1);
	        }
	     }
	if(!name_flag)
	  {
	   cerr << "Incorrect format of " << filename << "! Partition identifier is missing " 
		<< "in line " << linecount << ".\nFor further help type ./MARE -h.\n";
	   exit(1);
	  }
	else name_flag = 0;

//Ignore everything until a '-' appears including the beginning of a partition
        while(buf[i]!='-')
             {
              i++;
              if(i==buf.size())
                {
                 cerr << "Incorrect format of " << filename << "! '-' is missing in line "
                      <<  linecount << ".\nFor further help type ./MARE -h.\n";
                 exit(1);
		}
             }

	int j=0;

//Check if a numeral for the end of a partition and ';' appears
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
              if(i==buf.size())
                {
                 cerr << "Incorrect format of " << filename << "! ';' is missing in line "
                      <<  linecount << ".\nFor further help type ./MARE -h.\n";
                 exit(1);
	        }
             }
        if(strlen(str)==0)
           {
            cerr << "Incorrect format of " << filename << "! partition bound is missing in line "
                 <<  linecount << ".\nFor further help type ./MARE -h.\n";
            exit(1);
	   }
	else part_end = atoi(str);
       }
  if(!cont_flag)
    {
     cerr << "Incorrect format of " << filename << "! No content in charset file.\n"
	  << "For further help type ./MARE -h.\n";
     exit(1);
    }
  file.close();
  return part_end;
}

