
#include "uniqueKmerFinder.h"

//tested 
uniqueKmerFinder::uniqueKmerFinder()
{

}
uniqueKmerFinder::~uniqueKmerFinder()
{
  delete originalNames;
  delete hapNames;
  delete originalSize;
  delete hapSizes;
  delete originalSeq;
  delete hapSeq;
  delete allNames;
  originalNames = NULL;
  hapNames = NULL;
  originalSize= NULL;
  hapSizes= NULL;
  originalSeq= NULL;
  hapSeq= NULL;
  allNames = NULL;

}

void uniqueKmerFinder::findFileInformation(const string fName)
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
void uniqueKmerFinder::readFile(const string fName)
{
  cout << "reading file with address " << fName << endl ; 
	findFileInformation (fName);

	allSeq = new vector<string>;

	originalSize = new vector<int>;
	originalNames = new vector<string>;
	originalSeq = new vector<string>;

	hapSizes = new vector<int>;
	hapNames = new vector<string>;
	hapSeq = new vector<string>;


	ifstream genome;
	string line , newName, finalSeq ;//, reversed;
	int sToW , size ,index ,  l , file_idx; 


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

//categorizig information 
  	for (int i=0 ; i<allNames -> size() ; i++)
  	{
  		if (nameIsOriginal( (*allNames)[i]) )
		  {// original ones
  			originalNames -> push_back ((*allNames)[i]);
  			originalSize -> push_back ((*allSizes)[i]);
  			originalSeq -> push_back ((*allSeq)[i]);
		  }
		  else 
		  {
  			hapNames -> push_back ((*allNames)[i]);
  			hapSizes -> push_back ((*allSizes)[i]);
  			hapSeq -> push_back ((*allSeq)[i]);		
		  }
  	}

    //delete allNames;
    delete allSeq;
    delete allSizes;

    //allNames = NULL;
    allSeq = NULL; 
    allSizes = NULL;

  /*  cout << "haps:"<< endl;
    for (int i=0 ; i<hapNames->size() ; i++)
      cout << (*hapNames)[i] << " " << (*hapSeq)[i] << " " << (*hapSizes)[i] << endl; 

    cout << "originals:" << endl;
    for (int i=0 ; i< originalNames->size() ; i++)
      cout << (*originalNames)[i] << " " << (*originalSeq)[i] << " " << (*originalSize)[i] << endl; 
      */

}


void uniqueKmerFinder::generateUniqueKmersAllGenome(const int &k)
{
  vector<string> middleFiles , chromosome , chrHaps , finalFiles, finalChrHap;
  string toFind , chrNameForHap;
  ostringstream convert;
  system("mkdir \"res\"");
  //generate kmers for original chromosomes
  cout << "starts to generate different original chromosome kmers... " << endl;
  for (int i=0 ; i< originalNames -> size() ; i++)
    generateChrKmersAndWriteTofile( k , (*originalNames)[i] , (*originalSeq)[i]);

  //generate kmers for hap chromosome
  cout << "starts to generate different haplotypes kmers... " << endl;
  for (int i=0 ; i< hapNames-> size() ; i++)
    generateChrKmersAndWriteTofile(k , (*hapNames)[i] , (*hapSeq)[i]);

  //merging original chromosomes together
  int count = 1;
  for (int i=0 ; i< allNames -> size() ; i++)
  {
    convert <<  count;
    mapChrNameToString.push_back(">" + convert.str());
    rename ( ("res/"+(*allNames)[i] + ".txt").c_str() , ("res/"+mapChrNameToString[i] + ".txt").c_str()); 
    if (nameIsOriginal ( (*allNames)[i]))
      originalMap.push_back (mapChrNameToString[i]);
    else
      hapMap.push_back( mapChrNameToString[i] );
    count++;
    convert.str("");
  }
  
  cout << "merging original chromosomes STEP 1..."<< endl;
  mergeKmers ( middleFiles , originalMap , 1 );
  cout << " merging original chromosome finished" << endl ; 
  cout << "removing middle files STEP 1..." << endl;
  removeFiles ("" , middleFiles);
  middleFiles.clear();


  cout << "merging different haps of chromosome together STEP 2" << endl;
  //finding different chromosomes which have haplotypes sequence too
  for (int i=0 ; i<hapNames -> size() ; i++)
  {
    toFind = splitFirstPart((*hapNames)[i]);
    if (! findInVector(chromosome , toFind))
      chromosome.push_back(splitFirstPart((*hapNames)[i]));
  }


// finding each chromosome haplotypes and merge all of them together
  for (int i=0 ; i<chromosome.size() ; i++)
  {
    for (int j=0 ; j<hapNames->size() ; j++)
      if (splitFirstPart((*hapNames)[j]) == chromosome[i])
        chrHaps.push_back( hapMap[j] );
  
    cout << "merging " << splitFirstPart(chrHaps[0]) << " haps" << endl;

    // making final_chr_hap names map to their string final_>2_hap for example
    for (int j=0 ; j<chrHaps.size() ; j++)
      chrNameForHap = chrNameForHap + chrHaps[j] ;
    finalChrHap.push_back( chrNameForHap);

    mergeKmers (middleFiles , chrHaps , 2 );

    // removing middle files for merging haps
    cout << "removing middle files STEP 2..." << endl ; 
    removeFiles ("" , middleFiles);

    chrHaps.clear();
    middleFiles.clear();
    chrNameForHap = "" ; 
  }

  for (int i=0 ; i<finalChrHap.size() ; i++)
    finalChrHap[i] = "final_"+ finalChrHap[i] +  "_hap" ;
  //removing kmers of each chromosome from original and haps
  //HAPS
  cout << "removing hap chromosome kmer files..." << endl;
  for (int i=0 ; i<hapNames -> size() ; i++ )
    middleFiles.push_back( hapMap[i] + ".txt");
  removeFiles ("res/" , middleFiles);
  middleFiles.clear();

  //ORIGINAL
  for (int i=0 ; i<originalNames -> size() ; i++ )
    middleFiles.push_back(originalMap[i] + ".txt");
  removeFiles ("res/" , middleFiles);
  middleFiles.clear();


  cout << "merging chrs-with-haps together STEP 3..." << endl;
  mergeKmers (middleFiles , finalChrHap , 3 );
  cout << "removing middle files STEP 3..." << endl;
  removeFiles ("" , middleFiles);
  middleFiles.clear();

  for (int i=0 ; i<finalChrHap.size() ; i++)
    finalChrHap[i] = finalChrHap[i]+ ".txt" ;
  removeFiles ("res/" , finalChrHap);
  finalChrHap.clear();


  finalFiles.push_back("finalOriginal");
  finalFiles.push_back("finalHap");

  cout << "merging finalOriginal with finalHaps STEP 4 ..." << endl;
  mergeKmers (middleFiles , finalFiles , 4);
  for (int i=0 ; i<finalFiles.size() ; i++)
    finalFiles[i] = finalFiles[i]+ ".txt" ;

  //cout << "removing final files! :D " << endl;
  removeFiles("res/" , finalFiles);
  finalFiles.clear(); 
}

void uniqueKmerFinder::generateChrKmersAndWriteTofile (const int &k , const string &name , const string &seq)
{
  int cur_size ;
  string kMerSeq , new_Kmer_name , cur_seq , cur_name;
  Kmer new_Kmer;
  ostringstream convert_s , convert_e;
  set<Kmer> * my_Kmer = new set<Kmer>;

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


  new_Kmer.setUnique(1);
  new_Kmer.setName(cur_name);
  // generating unique k-mers
  for (int i=0 ; i<cur_size ; i++)
  { 
    if ((k_start[i]==1) && ((cur_size-i)>= (k)) )
    {//start position
      kMerSeq = cur_seq.substr (i,k); 
      new_Kmer.setSeq(kMerSeq);
      convert_s << i;
      convert_e << (i+k-1);

      new_Kmer.setStart(convert_s.str());
      new_Kmer.setEnd(convert_e.str());

      convert_s.str("");
      convert_e.str("") ;

      if (my_Kmer->count(new_Kmer)!=0)
      { // if k-mer is found in set, change the unique flag
          if ((my_Kmer->find(new_Kmer))-> getUnique())
          {// if it is the 1st time of repeatition of kmer change the flag
            my_Kmer->erase(my_Kmer->find(new_Kmer));
            Kmer obj (kMerSeq, 0 , cur_name , "0" , "0");
            my_Kmer->insert(obj);
          }
      }
      else
      { 
        //cout << "YES"<<endl; 
        my_Kmer->insert(new_Kmer);
      }

      kMerSeq = "";
    }

  }
  cout<<"kmers of "<<cur_name<<" finished. size of set (ALL) is "<< my_Kmer->size()<< "."<<endl; 

  //making the memory free
  delete[] k_start;
  k_start = NULL;

  writeKmersTofile (my_Kmer , cur_name);


  cout << "Deleting set of " << cur_name << " starts..." << endl; 
  delete my_Kmer;
  my_Kmer = NULL;

  cout << "Deleting set of " << cur_name << " finished." << endl ; 
  
}





void uniqueKmerFinder::writeKmersTofile (const set<Kmer> * myset , const string &name)
{

  ofstream f;
  set<Kmer>::iterator it ;
  string newName ;
  ostringstream convert , u ;
  long long ctr = 1;


  f.open ( ("res/" + name + ".txt").c_str());  
  cout << "writing all forward kmers to file..."<<endl;     

  for (it = myset->begin() ; it != myset->end() ; it++)
  {
    u << it->getUnique();
    convert << ctr;
    //if (R)
      //newName = it -> getName() + "|" + it->getStart() + "|" + it->getEnd() + "|R|" + u.str() + "|" + convert.str();
    //else
      newName = it -> getName()+ "|" + it->getStart() + "|" + it->getEnd() + "|" + u.str() + "|" + convert.str();

    convert.str("");
    u.str("");
    ctr++;
    f << newName << endl << it->getSeq() << endl;
  }

  cout << "Writing to file finished" << endl;
  f.close();
}

//original =1 -> merging original chrs together
//original =2 -> merging a chromosome different haplotypes together
//original =3 -> merging different chromosome-with-haps together leads to finalHap
//original =4 -> merging finalOriginal with finalHap leads to finalKmers
void uniqueKmerFinder::mergeKmers(  vector<string> &middleFiles , const vector<string> &names , const int &mergeType)
{

  vector<string> fun_names;
  ifstream f1 , f2;
  ofstream fout;
  string line1 , line2 , name1 , name2 , kmerName1 , kmerName2 , kmerName;
  ostringstream convert;
  long long ctr ;


  if (names.size()==1)
  {
    if (mergeType == 1)
      return;
    else if (mergeType ==2 ) // merging a chromosome different haps
    {
      rename ( ("res/"+names[0]+".txt").c_str() , ("res/final_"+splitFirstPart(names[0]) + "_hap.txt").c_str() );
      return;
    }
    else if (mergeType ==3)
    {
      rename (("res/"+names[0]+".txt").c_str() , "res/finalHap.txt" );
    }
    else
    {
      rename (("res/"+names[0]+".txt").c_str() , "res/putativeSunk" );
    }
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
        { // defining outpur name
          if (mergeType==1)
            fout.open("res/finalOriginal.txt");
          else if (mergeType ==2 )
            fout.open(("res/final_"+names[i] + names[i+1] + "_hap.txt").c_str());    
          else if (mergeType == 3)
            fout.open("res/finalHap.txt");    
          else
            fout.open("res/putativeSunk");
        } // defining output name finished

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
            if (line1 < line2)
            {
              if (mergeType == 4)
              {
                if (splitKmerId ('|' , name1 , 4) == "0")
                {
                  getline (f1 , line1);
                  convert.str("");
                  continue;
                }

              }
              ctr++;
              kmerName1 = splitKmerId('|' , name1 , 1) + "|" + splitKmerId('|' , name1 , 2) + "|" + splitKmerId('|' , name1 , 3) + "|" ;
              if (mergeType == 4)
                kmerName1 = kmerName1  + convert.str();
              else
                kmerName1 = kmerName1 + splitKmerId('|' , name1 , 4) + "|" + convert.str();
              fout<< kmerName1 << endl << line1 << endl;
              getline (f1 , line1);

            }
            else if (line1 > line2)
            {
              if (mergeType == 4)
              {
                if (splitKmerId ('|' , name2 , 4) == "0")
                {
                  getline (f2 , line2);
                  convert.str("");
                  continue;
                }

              }
              ctr++;
              kmerName2 = splitKmerId('|' , name2 , 1) + "|" + splitKmerId('|' , name2 , 2) + "|" + splitKmerId('|' , name2 , 3) + "|";
              if (mergeType == 4)
                kmerName2 = kmerName2  + convert.str();
              else
                kmerName2 = kmerName2 + splitKmerId('|' , name2 , 4) + "|" + convert.str();
              fout << kmerName2 << endl << line2 << endl;
              getline (f2 , line2);
            }
            else 
            {// kmers are same
              // check for unique and update start and end
              if (mergeType ==1 || mergeType==3) // merging original chromosomes and merging chrs-with-haps together
              {
                ctr++;
                kmerName = ">chr23|0|0|0|" + convert.str();
              }
              else if (mergeType ==2) // merging a chromosome different haps together
              {
                ctr++;
                if (splitKmerId('|' ,name1 , 4)=="1" &&  splitKmerId('|' ,name2 , 4) =="1")
                { // both kmer are unique in each haplotype
                  kmerName = splitKmerId('|',name1,1) + "," + splitKmerId('|',name2,1) + "|" + splitKmerId('|',name1,2) + "," + splitKmerId('|',name2,2) + "|";
                  kmerName = kmerName + splitKmerId('|',name1,3) + "," + splitKmerId('|',name2,3) + "|1|" +convert.str();
                }
                else
                  kmerName = ">chr23|0|0|0|" + convert.str();                
              }
              else // merging finalHap with finalOriginal to finalKmers
              {
                if (splitKmerId('|' ,name1 , 4)=="1" &&  splitKmerId('|' ,name2 , 4) =="1")
                { // both kmer are unique in haps and original
                  if (splitFirstPart(splitKmerId('|' , name1, 1 )) == splitFirstPart(splitKmerId('|' , name2, 1 )) )
                  {
                    ctr++;
                    kmerName = splitKmerId('|',name1,1) + "," + splitKmerId('|',name2,1) + "|" + splitKmerId('|',name1,2) + "," + splitKmerId('|',name2,2) + "|";
                    kmerName = kmerName + splitKmerId('|',name1,3) + "," + splitKmerId('|',name2,3) + "|" +convert.str();

                  }
                  else // not unique
                  {
                    getline (f1 , line1);
                    getline (f2 , line2);  
                    convert.str("");
                    continue;
                  }
                }   
                else // not unique
                {
                  getline (f1 , line1);
                  getline (f2 , line2);
                  convert.str("");
                  continue;
                }
              }

              fout << kmerName << endl << line1 << endl; 
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
          if (splitKmerId ('|' , name2 , 4)=="1")
          {
            kmerName2 = splitKmerId('|' , name2 , 1) + "|" + splitKmerId('|' , name2 , 2) + "|" + splitKmerId('|' , name2 , 3) + "|" ;
            if (mergeType == 4)
              kmerName2 = kmerName2 + convert.str();
            else
              kmerName2 = kmerName2 + splitKmerId('|' , name2 , 4) + "|" + convert.str();
          }
          else
          {
            if (mergeType == 4)
            {
              convert.str(""); 
              getline (f2 , line2);
              continue;
            }
            kmerName2 = ">chr23|0|0|0|" + convert.str();
          }
          convert.str("");
          fout << kmerName2 << endl << line2 <<endl;
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
          if (splitKmerId ('|' , name1 , 4)=="1")
          {
            kmerName1 = splitKmerId('|' , name1 , 1) + "|" + splitKmerId('|' , name1 , 2) + "|" + splitKmerId('|' , name1 , 3) + "|" ;
            if (mergeType == 4)
              kmerName1 = kmerName1 + convert.str();
            else
              kmerName1 = kmerName1 + splitKmerId('|' , name1 , 4) + "|" + convert.str();
          }
          else
          {
            if (mergeType == 4)
            {
              convert.str(""); 
              getline (f1 , line1);
              continue;
            }            
            kmerName1 = ">chr23|0|0|0|" + convert.str();
          }
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
    mergeKmers (middleFiles , fun_names , mergeType);
  }
}


bool uniqueKmerFinder::findInVector (vector<string> &v , string &x)
{
  for (int i=0 ; i< v.size() ; i++)
    if (v[i] == x)
      return 1;
  return 0;
}

void uniqueKmerFinder::stringToCharArray (char* &s , const string &seq ,const  int &start , const int &l)
{
  for (int i=start ; i<(start+l) ; i++)
    s[i] = seq[i-start];
}


string uniqueKmerFinder::splitLastPart(const string &s)
{
  for (int i=s.size()-1 ; i>=0 ; i--)
    if (s[i]=='_')
      return s.substr(i+1,s.size()-i-1);
  return s; 
}


string uniqueKmerFinder::splitFirstPart(const string &s)
{
  for (int i=0 ; i<s.size() ; i++)
    if (s[i]=='_')
      return s.substr(0,i);
  return s; 

}

bool uniqueKmerFinder::nameIsOriginal(const string &name)
{
  if (splitLastPart(name)=="random" || name.size()==5 || name.size()==6  || splitFirstPart(name) == ">chrUn")
    return 1;
  return 0;
}


string uniqueKmerFinder::splitKmerId (const char &identifier , const string &id ,const int &part)
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

void uniqueKmerFinder::removeFiles(const string  &dest, const vector<string> &names)
{
  for (int i=0 ; i<names.size() ; i++ )
  {
    remove ( (dest + names[i]).c_str());
    cout << names[i] << " deleted." << endl;
  }
}
