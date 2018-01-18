#ifndef GC_H
#define GC_H
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <stdlib.h>
using namespace std;
class genomicCoordinates
{
public:
    //adrees of .tab file
    string genomic_SD_add;
    //address of align_both directory
    string align_both_directory_add;


    genomicCoordinates (string, string);
    void reportMismatchIndelCoordinatesOfSDalignments();
private:
    //This function finds when alignment sequence starts
    void findSeq ( const string &id ,string &res);
    //This function splits an id with identifier 'c' and returns splited parts
    void splitWithChar (const char &identifier , const string &id ,vector<string> &res);
    //This function read the alignment file and returns seq1, seq2 and alignment pattern
    void getSDalignInfFromFile (const string &file_name , string &allSeqAlignPatt, string &allFirstSeq, string &allSecSeq);
    //This function reports mismatch based on alignment pattern 
    void reportMismatch (string &align_pattern,string &first_seq,string &sec_seq,const string &chr1,
    const string &chr2,const string &strand,const int &chr1_s,const int &chr2_s,const int &chr1_e,const int &chr2_e);
    //This function reports indels based on alignment pattern while considering seq1 and seq2
    void reportIndel (string &align_pattern,string &first_seq, string &sec_seq,const string &chr1,
    const string &chr2,const string &strand,const int &chr1_s,const int &chr2_s,const int &chr1_e , const int &chr2_e);
    //This function gets one row of 17 columns .tab file and extract basic informations
    void getAlignAttributes (string &align_file, string &chr1, string &chr2 , int &chr1_s, 
    int &chr1_e, int &chr2_s, int  &chr2_e , string &strand , const string &line_tab);



};


#endif