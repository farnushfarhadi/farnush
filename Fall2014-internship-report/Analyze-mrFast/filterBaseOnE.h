#ifndef FBOE_H
#define FBOE_H
#include <vector>
#include <fstream>
#include <sstream>
#include <set>
#include <stdlib.h>
#include "readInSamFile.h"
using namespace std;
class filterBaseOnE
{
  public:


    //constructors
  	filterBaseOnE ();

  	//This function filters the kmers based on mrfast results
    void filter(const bool &e, const int &mrfast_res_num, const string &add_reference, const string &add_mrfast_res, const int &limit);


   private:
   	//This function splits an id with identifier 'c' and returns splited parts
	void splitWithChar (const char &identifier , const string &id ,vector<string> &res);
	//This function get the reference names and mappes them to sum numbers
	void getRefNames(vector<string> &refNames , vector<string> &mapNames ,  const string &fName);
	//This function generates a set for each chromosome of mrfast_results
	void generateSetForEachChr(ofstream &chrFile, const bool forE0 , const int mrfast_res_number , const string &chrName , const string &add );
	//This function merge the kmers mapped in different chromosomes together while updating the 'Count' variable
	void merge(vector<string> &middleFiles , const vector<string> &names , const bool &forE0 , const int &mrfast_res_number ,const string &add);
	//This function buid the reverse complemented of a string, used for backward mapping
	void reverseString (string &seqR , const string &seq , const int &k );
	//This function remove all files specified in 2nd arg in relatve address 'dest'
	void removeFiles(const string  &dest, const vector<string> &names);
	//This functions works for fintering based on e=1 or2, it checks the count for > limit
	// and for e=0, it matches the 'Count' to number of kmer's names and decides it should be filtered or not
	bool checkCountLimit (const char &c , const string &name , const bool &forE0 , const int &limit);
	//This function calls checkCountLimit and categorizes kmers and write to files
	void checkKmersAndWriteToFile (const bool &forE0 , const int &num , const string &add, const int &limit);

};


#endif