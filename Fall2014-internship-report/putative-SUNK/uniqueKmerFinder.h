#ifndef UKF_H
#define UKF_H
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "Kmer.h"
using namespace std;
class uniqueKmerFinder
{
  public:
    //constructor
    uniqueKmerFinder();
    //distructor
    ~uniqueKmerFinder();
    //This function gets all information of the genome. categorizes names based on original or hap
    //stores all sequences needed for further processing
    void readFile(const string fName);
    //This function mainly does what is described in part 2.1
    void generateUniqueKmersAllGenome(const int &k);
  private:
    //storing ooriginal chromosome infs
    vector<string> * originalNames;
    vector<string> * originalSeq;
    vector<int> * originalSize;

    
    vector<string> * allNames;
    vector<string> * allSeq;
    vector<int> * allSizes;


    //storing hap chromosome infs
    vector<string> * hapNames;
    vector<int> * hapSizes;
    vector<string> * hapSeq;
    
    //storing mapped name of chromosomes 
    //because the merged name gets so big in size through merg function, 
    //we map the name for avoding to get output file's name problems
    vector<string> mapChrNameToString , originalMap , hapMap;

    //This function splits the string with '_' returns the last part
    string splitLastPart(const string &s);
    //This function splits the string with '_' returns the first part
    string splitFirstPart(const string &s);
    //This function splits an id with identifier 'c' and returns the part specified by 3rd arg
    string splitKmerId (const char &identifier , const string &id ,const int &part);
    //This function decide whether a chromosome names is considered as original or not
    bool nameIsOriginal(const string &name);
    //This function gets the genome information including chromosome names and size 
    void findFileInformation (const string fName);
    //This function convert a string to char* type
    void stringToCharArray (char* &s , const string &seq , const int &start , const int &l);
    //This function search for a element in vector
    bool findInVector (vector<string> &v , string &x);
    //This function remove all files specified in 2nd arg in relatve address 'dest'
    void removeFiles(const string  &dest, const vector<string> &names);
    //This function generate the set of locally unique kmers within chromosome and call 
    // writeKmersTofile function
    void generateChrKmersAndWriteTofile (const int &k , const string &name , const string &seq );
    //This function writes the privously generated set to a file following the format
    // described  in the documentation
    void writeKmersTofile (const set<Kmer> * myset , const string &name );
    //This function merges diffenerent kmers in 4 different methods/ described in documentation
    void mergeKmers(  vector<string> &middleFiles , const vector<string> &names , const int &original);
};
#endif