
#include "finalSunk.h"

finalSunk::finalSunk()
{

}



void finalSunk::readFileAndGenerateSet ( set<mersHash> &myset , const string &fileName , const int &k)
{

	string line , name , tmp;
	ifstream file;
	
	
	 
	file.open (fileName.c_str());
	while (!file.eof())
	{
		getline (file , line);
		if (line == "")
			continue;
		if (line[0]== '>')
		{
			name = line; 
			getline (file , line);
		}
		tmp = splitKmerId('|' , name , 1);
		int merCount = atoi( (tmp.substr(1, tmp.size()-1)).c_str() );
		bool * hash = hashKmer (line , k );  
    mersHash mer (merCount , k , hash );
		myset.insert (mer); 
	}


}


bool* finalSunk::hashKmer( const string kmer , const int &K )
{
  //A 00
  //C 01
  //G 10
  //T 11
  bool * hash_res = new bool[2*K];
  for (int i=0 ; i<kmer.size() ; i++)
  {
    if (kmer[i]=='A')
    {
      hash_res[2*i]=0;
      hash_res[2*i+1]=0;
    }
    else if(kmer[i]=='C')
    {
      hash_res[2*i]=0;
      hash_res[2*i+1]=1;      
    }
    else if(kmer[i]=='G')
    {
      hash_res[2*i]=1;
      hash_res[2*i+1]=0;      
    }
    else
    {
      hash_res[2*i]=1;
      hash_res[2*i+1]=1;      
    }
  }

  return hash_res;
}


string finalSunk::splitKmerId (const char &identifier , const string &id ,const int &part)
{
  int pre = 0;
  vector<string> res;
  for (int i=0 ; i<id.size() ; i++)
  {
    if (id[i] == identifier)
    {
      res.push_back(id.substr (pre , i-pre));
      pre = i+1;
    }
  }
  res.push_back(id.substr (pre , id.size()-pre));

  for (int i=0 ; i<res.size() ; i++)
    if (i+1 == part)
      return res[i];
}


void finalSunk::SUNKto15MersAndCalculateCount (const set<mersHash> &myset , const string &fileName , const int &k , const int &th)
{
  ifstream file; 
  ofstream file_out_ok , file_out_nok ;
  string line , name;
  file.open (fileName.c_str());
  system("mkdir \"res\"");
  file_out_ok.open ("res/final_sunk.fasta");
  file_out_nok.open ("res/sunkNot.fasta");
  while (!file.eof())
  {
    getline (file , line);
    if (line[0]=='>')
    {
      name = line; 
      getline (file , line );
      int sumCount = devideMerToSmallerMers (myset , k , line);
      if (sumCount < th)
      {
        file_out_ok << name << "|"<< sumCount << endl << line << endl; 
      }
      else
      {
        file_out_nok << name << "|"<< sumCount << endl << line << endl;
      }
    }
  }

  file.close();
}


int finalSunk::devideMerToSmallerMers (const set<mersHash> &myset , const int &smallMerK , const string &mer)
{
  int cur_count , sumCount;
  
  sumCount = 0 ; 
  for (int i=0 ; i<= mer.size()-smallMerK ; i++)
  {
    string merSub = mer.substr(i , smallMerK);
    bool * hash = hashKmer(merSub , smallMerK);
    mersHash obj (0 , smallMerK , hash);
    if (myset.count(obj)==0)
    { // since the k'-length and k-length strings are made of the same genome and k > k',
      // 15 mers are completely cover the k-length strings
      cout << "ERROR, 15 mers not complemte. " << merSub << " is not in the set!" << endl;
      break;
    }
    cur_count = myset.find(obj)-> countInGenome;
    sumCount += cur_count;
  }

  return sumCount;

}

void finalSunk::generateFinalSunk(const string &fin_15_mers, const string &fin_semi_sunks, const int &k, const int &th)
{
	set<mersHash> myset;
	readFileAndGenerateSet ( myset , fin_15_mers , k );
  	SUNKto15MersAndCalculateCount (myset , fin_semi_sunks , k, th );
}