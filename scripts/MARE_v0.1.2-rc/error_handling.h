#ifndef ERROR_HAND_H
#define ERROR_HAND_H

/***********************************************************************************

function check_linefeed
-----------------------

The function checks if the file uses Windows, Mac or Unix linefeeds. 

Parameters:

char* filename    : name of the file to be checked

Return values:

int 1 : Unix
int 2 : Mac
int 3 : Windows

************************************************************************************/

int check_linefeed(char* filename);

/***********************************************************************************

function check_fas
------------------

The function checks if the input file has correct Fasta format 

Parameters:

char* filename    : name of the file to be checked

Return values:

int : end of last partition

************************************************************************************/

int check_fas(char* filename);

/***********************************************************************************

/***********************************************************************************

function check_part
-------------------

The function checks if the input file has correct charset format

Parameters:

char* filename    : name of the file to be checked

Return values:

int : end of last partition

************************************************************************************/

int check_part(char* filename);

#endif
