#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>

using namespace std;
int main (int argc,char *argv[])
{


	int file_cap , file_ctr , line_ptr , parts , file_num;
	string add_input , line;
	ifstream fin;
	ofstream fout;
	ostringstream convert;
	add_input = argv[1];
	parts = atoi (argv[2]);

	fin.open (add_input.c_str());

	
	getline (fin , line );
	line_ptr = 1; 
	while (!fin.eof())
	{
		getline ( fin , line);
		if (line != "")
			line_ptr++;
	}
	fin.close();

	file_cap = (line_ptr/2)/(parts);
	fin.open (add_input.c_str());
	file_ctr = 0;
	file_num = 1;
	convert << file_num;
	fout.open ((add_input+"-"+convert.str()+".fasta").c_str());
	convert.str("");
	while (!fin.eof())
	{
		getline ( fin , line);
		fout << line << endl;
		getline (fin , line);
		fout << line << endl;
		file_ctr ++;
		if (file_num != parts)
		{
			if (file_ctr == file_cap)
			{
				file_ctr = 0 ;
				file_num++;
				fout.close();
				convert << file_num;
				fout.open ((add_input+"-"+convert.str()+".fasta").c_str());
				convert.str("");
			}
		}
	}
	fout.close();

	if (line_ptr % (2*parts) == 0)
	{
		convert << file_num;
		remove ((add_input+"-"+convert.str()+".fasta").c_str());
		convert.str("");
	}
	fin.close();


	return 0;
}