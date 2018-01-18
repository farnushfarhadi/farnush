#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
using namespace std;
int main(int argc , char* argv[])
{
	ofstream sunk;
	ifstream fin;
	ostringstream convert;
	string line;
	string reads_count_add = argv[1];
	string fout_sunk = argv[2];
	int num = atoi (argv[3]);

	sunk.open ((fout_sunk+"semi_final_sunk.fasta").c_str());
	for (int i=1  ; i<=num ; i++ )
	{
		convert << i;
		fin.open ( (reads_count_add+ "reads-count-"+ convert.str() + "OK.fasta").c_str() );

		while (! fin.eof())
		{
			getline (fin , line);
			if (line =="")
				continue;
			sunk << line << endl; 
		}
		convert.str("");
		fin.close ();
	}
	sunk.close();
	return 0;
}