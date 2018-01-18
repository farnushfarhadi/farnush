#include "readInSamFile.h"

readInSamFile::readInSamFile()
{

}

readInSamFile::readInSamFile(string n , string s, int c)
{
	name = n ;
	count = c ;
	seq = s;
}

bool readInSamFile :: operator < (const readInSamFile &object) const
{
	return (seq < object.seq );
}
