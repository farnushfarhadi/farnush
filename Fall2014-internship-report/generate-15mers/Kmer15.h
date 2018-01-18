#ifndef K15_H
#define K15_H
#include <string>
using namespace std;
class Kmer15
{
  public:
    //constructors
  	Kmer15 ();
    Kmer15 (string , unsigned long long);

    //getting class variables
  	string getSeq() const;
    unsigned long long getCount() const;

    //setting class variables
  	void setSeq(string );
    void setCount(unsigned long long);

    //defining operator
    bool operator < (const Kmer15 &) const;
    bool operator ==(const Kmer15 &) const;
  private:
    //kmer sequence
  	string seq;
    //kmer count
    long long count;


};


#endif