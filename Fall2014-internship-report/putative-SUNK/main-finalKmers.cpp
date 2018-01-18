#include "uniqueKmerFinder.h"
using namespace std;

int main(int argc,char *argv[])
{

  int k;
  string genome;

  k = atoi(argv[1]);
  genome = argv[2];
 
  uniqueKmerFinder * humanGenome = new uniqueKmerFinder;
  humanGenome -> readFile(genome);
  humanGenome -> generateUniqueKmersAllGenome( k );
  delete humanGenome;
  return 0;
}