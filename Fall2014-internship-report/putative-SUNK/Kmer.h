#ifndef KS_H
#define KS_H
#include <string>
using namespace std;
class Kmer
{
  public:
    //constructors
  	Kmer ();
    Kmer (string , bool , string  , string , string);

    //getting class variables
  	string getSeq() const;
    string getName() const;
    bool getUnique() const;
    string getStart() const;
    string getEnd() const;

    //setting class variables
  	void setSeq(string );
    void setName(string );
    void setUnique(bool );
    void setStart(string);
    void setEnd (string);

    //defining operator
    bool operator < (const Kmer &) const;
    bool operator ==(const Kmer &) const;
  private:
    //k length string as a k-mer
  	string seq;
    //chromosome name from which kmer is coming
    string chr_name;
    //unique flag shows the kmer is locally unique within chromosme
  	bool unique;
    //0-based starting position of kmer
    string start;
    //0-based ending position of kmer
    string end;


};


#endif