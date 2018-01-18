
#include "Kmer15.h"

Kmer15::Kmer15()
{

}



Kmer15::Kmer15(string seq_ , unsigned long long count_)
{
  seq = seq_;
  count = count_;
}


string Kmer15::getSeq() const
{
  return seq;
}


unsigned long long Kmer15::getCount() const
{
  return count;
}
void Kmer15::setSeq(string s)
{
  seq=s ;
}


void Kmer15::setCount(unsigned long long c)
{
  count = c;
}

bool Kmer15 :: operator < (const Kmer15 & object) const
{
  return (seq < object.seq);
}

bool Kmer15 ::  operator ==(const Kmer15 & object) const
{
  return (seq == object.seq);  
}