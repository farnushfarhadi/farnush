
#include "mersHash.h"

mersHash::mersHash()
{

}



mersHash::mersHash(int c , int k_ , bool* h)
{
  countInGenome = c;
  k = k_;
  hashSeq = h;
}



bool mersHash :: operator < (const mersHash & object) const
{
  for (int i=0 ; i<2*k ; i++)
  {
    if (object.hashSeq[i]==hashSeq[i])
      continue;
    else if (hashSeq[i] < object.hashSeq[i])
      return true;
    else
      return false;
  }
  return false;
}

