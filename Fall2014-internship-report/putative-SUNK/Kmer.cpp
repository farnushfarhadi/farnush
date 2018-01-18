
#include "Kmer.h"

Kmer::Kmer()
{

}



Kmer::Kmer(string seq_ , bool unique_ , string name , string s, string e)
{
  unique = unique_;
  seq = seq_;
  chr_name=name;
  start = s;
  end = e;
}


string Kmer::getSeq() const
{
  return seq;
}

string Kmer::getName() const
{
  return chr_name;
}

bool Kmer::getUnique() const
{
  return unique;
}

string Kmer::getStart() const
{
  return start;
}

string Kmer::getEnd() const
{
  return end;
}

void Kmer::setSeq(string s)
{
  seq=s ;
}

void Kmer::setName(string n)
{
  chr_name=n ;
}

void Kmer::setUnique(bool u)
{
  unique = u;
}

void Kmer::setStart(string s)
{
  start = s;
}

void Kmer::setEnd(string e)
{
  end = e ;
}


bool Kmer :: operator < (const Kmer & object) const
{
  return (seq < object.seq);
}

bool Kmer ::  operator ==(const Kmer & object) const
{
  return (seq == object.seq);  
}
