
#include "mersFinder.h"
//tested 

mersFinder::mersFinder()
{

}


mersFinder::~mersFinder()
{
  delete allNames;
  delete allSizes;
  delete allSeq;
  allNames = NULL;
  allSizes = NULL;
  allSeq = NULL;
}
//tested
void mersFinder::findFileInformation(const string fName)
{
	allNames  = new vector<string>;
	allSizes = new vector<int>;
	ifstream genome;
	string line;
	int size = 0 ; 
	genome.open (fName.c_str());
	while (!genome.eof())
  {
    getline (genome , line);
    while (line[0]=='>')
    {
      allNames->push_back(line);
      getline (genome , line);
      while (line[0]!='>' )
      {
        if (genome.eof())
        {
          size += line.size();
          break;
        }
        else
        {
          size += line.size();
          getline (genome , line);
        }
      }

      allSizes->push_back(size);
      size = 0;
    }
  }

  genome.close();
}

//tested 
void mersFinder::readFile(const string fName)
{
	findFileInformation (fName);

	allSeq = new vector<string>;

	ifstream genome;
	string line , newName, finalSeq ;//, reversed;
	int sToW , size ,index ,  l , file_idx , count; 
  ostringstream convert;

	genome.open (fName.c_str());
  	while (!genome.eof())
  	{
    	sToW = 0 ;
    	file_idx = 0 ; 
    	getline (genome , line);
    	//choosing _hap chromosome + original chrs which their hap are not availabele + random
    	while (line[0]=='>')
    	{ 
      		size = (*allSizes)[file_idx] ; 
      		char* seq = new char[size];
      		getline(genome,line);
      		while (line[0]!='>' )
      		{
        		l = line.size();
        		if (genome.eof())
        		{
          			stringToCharArray (seq , line , sToW , l);
          			sToW += line.size();
          			break;
       			}
        		else
        		{
          			stringToCharArray (seq , line , sToW , l);
          			sToW +=  line.size();
          			getline (genome , line);
        		}
      		}

      		finalSeq = string(seq).substr(0, sToW);
      		allSeq->push_back(finalSeq);
      		delete[] seq;
      		seq = NULL;
      		sToW = 0;
      		file_idx++;
    	}
  	}

  	genome.close();

  count = 1;
  for (int i=0 ; i< allNames -> size() ; i++)
  {
    convert <<  count;
    mapChrNameToString.push_back(">" + convert.str());
    count++;
    convert.str("");
  }


}


void mersFinder::generateKmersAllGenome(const int &k)
{
  vector<string> middleFiles, originalNamesR , chromosome , chrHaps , finalFiles, finalChrHap;
  string toFind , chrNameForHap;
  ostringstream convert;

  system("mkdir \"res\"");
  //generate kmers for all chromosomes
  cout << "starts to generate different chromosome kmers... " << endl;
  for (int i=0 ; i< allNames-> size() ; i++)
    generateChrKmersAndWriteTofile(k , (*allNames)[i] , (*allSeq)[i] );

  //merging chromosomes together
  for (int i=0 ; i< allNames -> size() ; i++)
    rename ( ("res/"+(*allNames)[i] + ".txt").c_str() , ("res/"+mapChrNameToString[i] + ".txt").c_str()); 


  cout << "merging chromosomes ..."<< endl;
  mergeKmers ( middleFiles , mapChrNameToString);

  cout << "removing middle files ..." << endl;
  removeFiles ("" , middleFiles);


  for (int i=0 ; i<mapChrNameToString.size() ; i++)
    mapChrNameToString[i] = mapChrNameToString[i] + ".txt";
  cout << "removing chromosomes 15mers ..." << endl; 
  removeFiles ("res/" , mapChrNameToString); 



}

void mersFinder::generateChrKmersAndWriteTofile (const int &k , const string &name , const string &seq )
{
  unsigned long long cur_count;
  int cur_size ;
  string kMerSeq , new_Kmer_name , cur_seq , cur_name;
  Kmer15 new_Kmer;
  set<Kmer15> * my_Kmer = new set<Kmer15>;

  cur_name = name;
  cur_size = seq.size();
  cur_seq = seq;

  cout<<"starting to generate kmers of "<<cur_name<<" ..."<<endl;

  bool* k_start = new bool [cur_size];
  for (int i=0 ; i<cur_size ; i++)
   k_start[i]=1;


  //scanning sequence for N and finding candidate of k-Mers starting position
  for (int i=0; i<cur_size ; i++)
    if ( cur_seq[i]=='N')
        for (int j=max(i-k+1,0) ; j<=i ;j++)
          k_start[j]=0;


  new_Kmer.setCount(1);
  // generating unique k-mers
  for (int i=0 ; i<cur_size ; i++)
  { 
    if ((k_start[i]==1) && ((cur_size-i)>= (k)) )
    {//start position
      kMerSeq = cur_seq.substr (i,k); 
      new_Kmer.setSeq(kMerSeq);

      if (my_Kmer->count(new_Kmer)!=0)
      { // if k-mer is found in set, increase count by one 
        cur_count = (my_Kmer->find(new_Kmer))-> getCount();
        my_Kmer->erase(my_Kmer->find(new_Kmer));
        cur_count++;
        Kmer15 obj (kMerSeq , cur_count);
        my_Kmer -> insert (obj);
      }
      else
        my_Kmer->insert(new_Kmer);

      kMerSeq = "";
    }

  }

  cout<<"kmers of "<<cur_name<<" finished."<<endl;
  cout << "size of set (ALL) is "<< my_Kmer->size()<< "."<<endl;      

  //making the memory free
  delete[] k_start;
  k_start = NULL;

  writeKmersTofile (my_Kmer , cur_name );


  cout << "Deleting set of " << cur_name << " starts..." << endl; 
  delete my_Kmer;
  my_Kmer = NULL;

  cout << "Deleting set of " << cur_name << " finished." << endl ; 
  
}





void mersFinder::writeKmersTofile (const set<Kmer15> * myset , const string &name )
{

  ofstream f;
  set<Kmer15>::iterator it ;
  ostringstream convert ;
  unsigned long long ctr = 1;

  f.open ( ("res/" + name + ".txt").c_str());  
  cout << "writing " << myset->size() << " kmers to file..."<< endl ;

  for (it = myset->begin() ; it != myset->end() ; it++)
  {

    convert << ctr;
    ctr++;

    f << ">" << it -> getCount() << "|" << convert.str() << endl << it->getSeq() << endl;
    //f << "|";
    //f << convert.str() << endl << it->getSeq() << endl;
    convert.str("");
  }

  cout << "Writing to file finished" << endl;
  f.close();
}

void mersFinder::mergeKmers(  vector<string> &middleFiles , const vector<string> &names)
{

  vector<string> fun_names;
  ifstream f1 , f2;
  ofstream fout;
  string line1 , line2 , name1 , name2 , kmerName1 , kmerName2 , kmerName , c1_str , c2_str;
  ostringstream convert;
  unsigned long long ctr , c1 , c2  , c;


  if (names.size()==1)
  {
    return;
  }
  else
  {

    for (int i=0 ; i<names.size() ; i= i+2)
    {
      ctr = 1;
      if (i== names.size() -1 )
      {
        fun_names.push_back(names[i]);
        cout << names[i] << " not merged yet." << endl;
      }
      else
      {
        f1.open (("res/"+names[i]+".txt").c_str());
        f2.open (("res/"+names[i+1]+".txt").c_str());
        if (names.size()==2)
          fout.open("res/final15Mers.fasta");
        else
        {
          middleFiles.push_back("res/"+names[i]+names[i+1]+".txt");
          fout.open(("res/"+names[i]+names[i+1]+".txt").c_str());
        }

        //reading file starts
        getline (f1 , line1 );
        getline (f2 , line2 );
        while (!f1.eof() && !f2.eof())
        {
          if (line1=="" && line2=="")
            continue;

          if (line1[0]=='>' &&  line2[0]=='>')
          {
            name1 = line1;
            name2 = line2;
            getline (f1 , line1);
            getline (f2 , line2);
            continue;
          }

          else if (line2[0]=='>')
          {
            name2 = line2;
            getline (f2 , line2);
            continue;
          }

          else if (line1[0]=='>')
          {
            name1 = line1 ; 
            getline (f1 , line1);
            continue;
          }
          else // here we write to output files
          {
            convert << ctr;
            ctr++;
            if (line1 < line2)
            {
              kmerName1 = splitKmerId('|' , name1 , 1) + "|" + convert.str() ;
              fout<< kmerName1 << endl << line1 << endl;
              getline (f1 , line1);

            }
            else if (line1 > line2)
            {
              kmerName2 = splitKmerId('|' , name2 , 1) + "|" + convert.str() ;
              fout << kmerName2 << endl << line2 << endl;
              getline (f2 , line2);
            }
            else 
            {// kmers are same
              c1_str =  splitKmerId('|' , name1 , 1) ;
              c1_str = c1_str.substr(1 , c1_str.size()-1 );
              c2_str =  splitKmerId('|' , name2 , 1) ; 
              c2_str = c2_str.substr(1 , c2_str.size()-1 );
              c1 = atoll(c1_str.c_str());
              c2 = atoll( c2_str.c_str());
              c = c1 + c2;
              fout << ">" << c << "|" << convert.str() << endl << line1 << endl; 
              getline (f1 , line1);
              getline (f2 , line2);

            }
            convert.str("");
          }
        }

        while (!f2.eof())
        {
          if (line2 == "")
          continue;

          if (line2[0]=='>')
          {// now line2 is a name
            name2 = line2;
            getline (f2 , line2);
          }
          convert << ctr;
          ctr++;
          
          kmerName2 = splitKmerId('|' , name2 , 1) + "|" + convert.str() ;
          fout << kmerName2 << endl << line2 <<endl; 
          convert.str("");
          getline (f2 , line2);
          
        }

        while (!f1.eof())
        {
          if (line1 == "")
          continue;
          // now line1 is a name, it happens when end file2 is similiar to another one of file1 kmer. getline in file1 gives a name
          // but when end of file2 is < one of file1 kmer, line1 is that kmer now and the name is saved in name1
          if (line1[0]=='>')
          {// now line1 is a name
            name1 = line1;
            getline (f1 , line1);
          }
          convert << ctr;
          ctr++;
          kmerName1 = splitKmerId('|' , name1 , 1) + "|" + convert.str() ;
          convert.str("");
          fout << kmerName1 << endl << line1 << endl;
          getline (f1 , line1);
        }

        fun_names.push_back(names[i] + names[i+1]);
        cout << "merging " << names[i] << " and " << names[i+1] << " finished" << endl;
      }
      fout.close() ;
      f1.close();
      f2.close();
    }
    mergeKmers (middleFiles , fun_names );
  }
}


void mersFinder::stringToCharArray (char* &s , const string &seq ,const  int &start , const int &l)
{
  for (int i=start ; i<(start+l) ; i++)
    s[i] = seq[i-start];
}


string mersFinder::splitLastPart(const string &s)
{
  for (int i=s.size()-1 ; i>=0 ; i--)
    if (s[i]=='_')
      return s.substr(i+1,s.size()-i-1);
  return s; 
}


string mersFinder::splitFirstPart(const string &s)
{
  for (int i=0 ; i<s.size() ; i++)
    if (s[i]=='_')
      return s.substr(0,i);
  return s; 

}


string mersFinder::splitKmerId (const char &identifier , const string &id ,const int &part)
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

void mersFinder::removeFiles(const string  &dest, const vector<string> &names)
{
  for (int i=0 ; i<names.size() ; i++ )
  {
    remove ( (dest + names[i]).c_str());
    cout << names[i] << " deleted." << endl;
  }
}


