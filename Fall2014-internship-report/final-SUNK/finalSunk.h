#ifndef FS_H
#define FS_H
#include <string>
#include <iostream>
#include "mersHash.h"
using namespace std;
class finalSunk
{

public:


    //constructors
  	finalSunk ();


    void generateFinalSunk(const string &fin_15_mers, const string &fin_semi_sunks, const int &k, const int &th);
private:
    // This function hashes a k-length string made of A,C,G and T to a 2*k-length bool array
    bool* hashKmer( const string kmer , const int &k );
    //This function generate a 'mersHash' set from all 15 mers. 
    //it hashes the sequence of 15 mers to use lower amount of memory
    //The set is stored in myset address
    void readFileAndGenerateSet (set<mersHash> &myset ,  const string &fileName , const int &k);
    //This function split an id with identifier 'c' and returns the part specified by 3rd arg 
    string splitKmerId (const char &identifier , const string &id ,const int &part);
    //This function reads all SUNKs, sums its frequencies and then decide for the final SUNK based on th
    void SUNKto15MersAndCalculateCount (const set<mersHash> &myset , const string &fileName , const int &k , const int &th);
    //This function devides a k-length string to k'-length strings starting from from first and sliding by one char
    // and hashes k'-length strings, then finds its count variable in myset which contains all k'-length in their hashed format
    // at the end, it sums all count related to k'-length strings
    int devideMerToSmallerMers (const set<mersHash> &myset , const int &smallMerK , const string &mer);



};


#endif