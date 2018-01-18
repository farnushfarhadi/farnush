#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include "genomicCoordinates.h"
using namespace std;

int main(int argc , char *argv[])
{
	string genomic_SD_add = argv[1];
	string align_both_directory_add = argv[2];

	genomicCoordinates human( genomic_SD_add , align_both_directory_add);
	human.reportMismatchIndelCoordinatesOfSDalignments();



	return 0;
}

