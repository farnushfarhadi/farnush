

#ifndef MF_H
#define MF_H
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#include <stdlib.h>


#include "Kmer15.h"

using namespace std;
class mersFinder
{
  public:

    //constructor
    mersFinder();
    //distructor
    ~mersFinder();
    //This function gets all information of the genome including chromosome name,
    // chromosome size and chromosome sequences.
    void readFile (const string fName);
    //This function maily does what is described in 'Generate-15mers' of part 2.3
    void generateKmersAllGenome(const int &k);

  private:

    //storing chromosome infs
    vector<string> * allNames;
    vector<int> * allSizes;
    vector<string> * allSeq;
    //storing mapped name of chromosomes 
    //because the merged name gets so big in size through merg function, 
    //we map the name for avoding to get output file's name problems
    vector<string> mapChrNameToString ;


    //This function splits the string with '_' returns the last part
    string splitLastPart(const string &s);
    //This function splits the string with '_' returns the first part
    string splitFirstPart(const string &s);
    //This function splits an id with identifier 'c' and returns the part specified by 3rd arg
    string splitKmerId (const char &identifier , const string &id ,const int &part);
    //This function gets the genome information including chromosome names and size 
    void findFileInformation (const string fName);
    //This function convert a string to char* type
    void stringToCharArray (char* &s , const string &seq , const int &start , const int &l);
    //This function remove all files specified in 2nd arg in relatve address 'dest'
    void removeFiles(const string  &dest, const vector<string> &names);
    //This function generate the set of all k-mers within chromosome and call 
    // writeKmersTofile function
    void generateChrKmersAndWriteTofile (const int &k , const string &name , const string &seq);
    //This function writes the privously generated set to a file following the format
    // described  in the documentation
    void writeKmersTofile (const set<Kmer15> * myset , const string &name);
    //This function merges diffenerent kmers in the way described in documentation
    void mergeKmers(  vector<string> &middleFiles , const vector<string> &names);

   


};
#endif