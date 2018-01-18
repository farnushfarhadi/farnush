#include "mersFinder.h"
using namespace std;
int main(int argc,char *argv[])
{

  int k;
  string add_genome; 
  add_genome = argv[1];
  k = atoi(argv[2]);

  mersFinder * humanGenome = new mersFinder;

  humanGenome -> readFile(add_genome);
  humanGenome -> generateKmersAllGenome(k);

  delete humanGenome;
}