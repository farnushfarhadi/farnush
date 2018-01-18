#include "filterBaseOnE.h"
using namespace std;
int main(int argc,char *argv[])
{

	string add_reference  , add_mrfast_res;
    int mrfast_res_num , limit;
    bool e;

	e = atoi(argv[1]);
	add_reference = argv[2];
	add_mrfast_res = argv[3];
	mrfast_res_num = atoi(argv[4]);
	limit = atoi (argv[5]);
	
	filterBaseOnE mrfast_res;
	mrfast_res.filter( e, mrfast_res_num , add_reference , add_mrfast_res , limit);

	return 0;
}
