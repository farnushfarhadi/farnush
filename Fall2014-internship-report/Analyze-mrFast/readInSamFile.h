#ifndef RISF_H
#define RISF_H
#include <string>
#include <iostream>
using namespace std;
class readInSamFile
{
  public:
  	//kmer name and seq
    string name , seq;
    //kmer count for mapping
    int count ;

    //constructors
  	readInSamFile ();
    readInSamFile (string , string , int  );


    //defining operator
    bool operator < (const readInSamFile &) const;

};


#endif