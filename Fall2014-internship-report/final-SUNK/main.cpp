#include "finalSunk.h"
using namespace std;
int main (int argc , char* argv[])
{

	string fin_15_mers = argv[1];
  string fin_semi_sunks = argv[2];
  int k = atoi (argv[3]);
  int th = atoi (argv[4]);

  finalSunk object;
  object.generateFinalSunk(fin_15_mers, fin_semi_sunks, k, th);
	return 0;
}


