#ifndef MHS_H
#define MHS_H
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <stdlib.h>
using namespace std;
class mersHash
{
  public:
    //the 'Count' attribute of a 15-mer
    int countInGenome;
    // the hashed value of a 15-mer
    bool* hashSeq;
    int k;

    //constructors
  	mersHash ();
    mersHash (int , int , bool*);


    //defining operator
    bool operator < (const mersHash &) const;
    //assign a array to hashSeq
    void assignArray(const bool* );

};


#endif