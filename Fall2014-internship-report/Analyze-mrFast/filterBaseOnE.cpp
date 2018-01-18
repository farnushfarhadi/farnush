#include "filterBaseOnE.h"

filterBaseOnE::filterBaseOnE()
{

}


void filterBaseOnE::filter(const bool &e, const int &mrfast_res_num, const string &add_reference, const string &add_mrfast_res , const int &limit)
{
	vector<string> refNames , mapNames , mapNames_temp, middleFiles , names;
	ostringstream convert; 
	ofstream chrFile;

	system( ("mkdir \""+add_mrfast_res+"e0\"").c_str() );
	system( ("mkdir \""+add_mrfast_res+"e1-2\"").c_str() );

	//get the information of the reference (which the reads were mapped in)
	getRefNames (refNames , mapNames, add_reference);

	//For each of 8 output files of mrFast result
	for (int j = 1 ; j<=mrfast_res_num ; j++)
	{
		convert << j ;
		//write reads mapped to refName[i] in seperate files.
		for (int i=0 ; i< refNames.size() ; i++)
		{
			 
			if (!e)
				chrFile.open ( (add_mrfast_res+"e0/"+mapNames[i] +"-" + convert.str()+ ".fasta").c_str()  );
			else
				chrFile.open ( (add_mrfast_res+"e1-2/"+mapNames[i] +"-" +convert.str() + ".fasta").c_str()  );

			//name of the files to be merged
			names.push_back(mapNames[i] +"-" +convert.str());
			
			cout << "generating set for " << mapNames[i] << " : " << refNames[i] <<" of mrFast res " << j << endl; 
			generateSetForEachChr( chrFile, (!e) , j , refNames[i] , add_mrfast_res);
			chrFile.close();
		}

		//merge reads mapped to different chromosomes (different files) together
		cout << "merging starts..." << endl; 
		merge (middleFiles , names , (!e) , j , add_mrfast_res);

		//removing unnecessary files
		cout << "removing middle files ..." << endl;
		 if (!e)
		  	removeFiles ((add_mrfast_res+"e0/") , middleFiles);
		else
			removeFiles ((add_mrfast_res+"e1-2/") , middleFiles);

		for (int i=0 ; i<mapNames.size() ; i++)
		    mapNames_temp.push_back(mapNames[i] +"-" +convert.str()+ ".fasta");
		cout << "removing chromosomes reads ..." << endl; 
		if (!e)
		  	removeFiles ((add_mrfast_res+"e0/") , mapNames_temp);
		else
			removeFiles ((add_mrfast_res+"e1-2/") , mapNames_temp);
		mapNames_temp.clear();
		convert.str("");
		names.clear();

		//finding those kmers which cant pass this level of filtering
		checkKmersAndWriteToFile ((!e) , j , add_mrfast_res, limit);
		for (int j = 1 ; j<=mrfast_res_num ; j++)
		{
			convert << j;
			if (e)
			{
				remove ( (add_mrfast_res+"e1-2/reads-count-"+convert.str()+".fasta").c_str() );
			}
			else
			{
				remove ( (add_mrfast_res+"e0/reads-count-"+convert.str() + ".fasta").c_str() );
			}
			convert.str("");
		}
	}

	/*
//For checking maximum number of mapping
	vector<string> res;
	ifstream file_in;
	//ostringstream convert;
	string line ;
	int maxC = 1;
	for (int j = 1 ; j<=res_num ; j++)
	{
		int maxC_file=1; 
		convert << j;
		file_in.open ( (add+"e1-2/reads-" +convert.str() + ".txt").c_str()  );

		while (!file_in.eof())
		{
			getline (file_in , line);
			if (line[0]== '>')
			{
				splitWithChar ('|' , line , res);
				int c = atoi (res[3].c_str());
				if (c > maxC_file)
				{
					maxC_file = c;
				}
				res.clear();
			}

		}

		cout << "maximum count after reads-" << convert.str() << ".txt was " << maxC_file << endl ; 
		convert.str("");
		file_in.close();
	}

*/

}


void filterBaseOnE::splitWithChar (const char &identifier , const string &id ,vector<string> &res)
{
  int pre = 0;
  for (int i=0 ; i<id.size() ; i++)
  {
    if (id[i] == identifier)
    {
      res.push_back(id.substr (pre , i-pre));
      pre = i+1;
    }
  }
  res.push_back(id.substr (pre , id.size()-pre));

}

//This function extract sequence names from a fasta file and also maps the names to numbers fo further analysis.
void filterBaseOnE::getRefNames(vector<string> &refNames , vector<string> &mapNames, const string &fName)
{
	ifstream genome;
	string line;
	int size = 0 ; 
	int count = 1; 
	ostringstream convert; 
	genome.open (fName.c_str());
	while (!genome.eof())
	{
	  getline (genome , line);
	  if (line[0]=='>')
	  {
	  	convert << count;
	    refNames.push_back(line.substr(1,line.size()-1));
	    mapNames.push_back (">"+ convert.str());
	    count++;
	    convert.str("");
	  }
	}
  genome.close(); 
}

//This function generate a set from reads of out-mrfast_res_number which are mapped to chrName and writes it to a file (chrFile) with specific name.
void filterBaseOnE::generateSetForEachChr(ofstream &chrFile, const bool forE0 , const int mrfast_res_number , const string &chrName , const string &add )
{
	ifstream sam;
	ostringstream convert;
	string line , readN , readS , readChr ;
	int cur_count , map_FLAG;
	char nm;
	bool chrPass;
	vector<string> line_res , name_res;
	readInSamFile cur_read;
	set<readInSamFile>::iterator it;


	chrPass = 0 ; 
	set<readInSamFile> myset;
	convert << mrfast_res_number;
	sam.open ( (add+"out"+convert.str()).c_str() );
	while (! sam.eof())
	{
		getline (sam , line);
		if (line[0]=='@' || line=="")
			continue; 
		splitWithChar('\t' , line , line_res);

		readN = line_res[0];
		readChr = line_res[2];
		map_FLAG = atoi (line_res[1].c_str());


		if (readChr!= chrName)
		{
			if (chrPass)
			{
				line_res.clear();
				break;
			}
			else
			{
				line_res.clear(); 
				continue;
			}
		}
		nm = (line_res[11])[5];
		//checking fwrd mapping or reversed
		if (map_FLAG==16)
		{
			reverseString (readS , line_res[9] , int (line_res[9].size()));
		}
		else
		{
			readS = line_res[9];
		}


		if (forE0)
		{
			if (nm=='1' || nm=='2')
			{ 
				line_res.clear();
				continue;
			}

		}
		else
		{
			if (nm=='0')
			{
				line_res.clear();
				continue;
			}
			else
				nm = '1';
		}


		chrPass = 1; 
		cur_read.name = readN ;
		cur_read.seq = readS;
		cur_read.count=1; 

		// generating the set
		if (myset.count (cur_read)!= 0 )
		{
			cur_count = (myset.find(cur_read))->count;
			cur_count++; 
	        myset.erase(myset.find(cur_read));
	        readInSamFile obj (readN , readS , cur_count   );
	        myset.insert (obj);
		}
		else
		{
			myset.insert (cur_read);
		}
			
		line_res.clear();
		convert.str("");
	}


	cout << "writing for " << chrName << " of file out" << mrfast_res_number <<  " set size: " << myset.size() << endl; 
	for (it=myset.begin() ; it!= myset.end() ; it++)
	{
		splitWithChar ('|' , it->name ,name_res);
		chrFile << ">"<<name_res[0] << "|" << name_res[1]<<"|" << name_res[2] << "|"<<it->count << endl << it->seq << endl;
		name_res.clear();
	}
	myset.clear();
	sam.close();


}


void filterBaseOnE::merge(  vector<string> &middleFiles , const vector<string> &names , const bool &forE0 , const int &mrfast_res_number , const string &add)
{

  vector<string> fun_names, id_res1 , id_res2;
  ifstream f1 , f2;
  ofstream fout;
  string line1 , line2 , name1 , name2 , kmerName1 , kmerName2 , kmerName ;
  ostringstream convert;
  int c1 , c2  , c;


  if (names.size()==1)
  {
    return;
  }
  else
  {

    for (int i=0 ; i<names.size() ; i= i+2)
    {
      if (i== names.size() -1 )
      {
        fun_names.push_back(names[i]);
        cout << names[i] << " not merged yet." << endl;
      }
      else
      {
      	if (forE0)
      	{
	        f1.open ((add+"e0/"+names[i]+".fasta").c_str());
	        f2.open ((add+"e0/"+names[i+1]+".fasta").c_str());
    	}
    	else
    	{ 
	        f1.open ((add+"e1-2/"+names[i]+".fasta").c_str());
	        f2.open ((add+"e1-2/"+names[i+1]+".fasta").c_str());    		
    	}
        if (names.size()==2)
        {
        	if (forE0)
        	{
        		convert << mrfast_res_number;
          		fout.open( (add+"e0/reads-count-"+convert.str()+".fasta").c_str() );
          		convert.str("");
        	}
          	else
          	{
          		convert << mrfast_res_number;
          		fout.open ( (add+"e1-2/reads-count-"+convert.str()+".fasta").c_str() );
          		convert.str("");
          	}
        }
        else
        {
          middleFiles.push_back(names[i]+names[i+1]+".fasta");
          if (forE0)
          	fout.open((add+"e0/"+names[i]+names[i+1]+".fasta").c_str());
          else
          	fout.open((add+"e1-2/"+names[i]+names[i+1]+".fasta").c_str());
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
            if (line1 < line2)
            {
                fout<< name1 << endl << line1 << endl;
                getline (f1 , line1);

            }
            else if (line1 > line2)
            {
                fout << name2 << endl << line2 << endl;
                getline (f2 , line2);
            }
            else 
            {// kmers are same
            	splitWithChar ('|' , name1 , id_res1);
            	splitWithChar ('|' , name2 , id_res2);
            	//if ( (id_res1[0]+id_res1[1]+id_res1[2]) == (id_res2[0]+ id_res2[1]+ id_res2[2]) )
            	//{
	              c1 = atoi(id_res1[3].c_str());
	              c2 = atoi( id_res2[3].c_str());
	              c = c1 + c2;
	              cout << c1 << " " << c2 << " " << c << endl;
	              fout << id_res1[0] << "|" << id_res1[1] << "|" << id_res2[2] << "|" << c << endl << line1 << endl; 
	              id_res2.clear();
	              id_res1.clear();
	              getline (f1 , line1);
	              getline (f2 , line2);
            	//}
            	//else
            		//cout << "same seq BUT different read IDs, ERROR" << endl ; 


            }
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
          
          fout << name2 << endl << line2 <<endl; 
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
          fout << name1 << endl << line1 << endl;
          getline (f1 , line1);
        }

        fun_names.push_back(names[i] + names[i+1]);
        cout << "merging " << names[i] << " and " << names[i+1] << " finished" << endl;
      }
      fout.close() ;
      f1.close();
      f2.close();
    }
    merge (middleFiles , fun_names  , forE0 , mrfast_res_number , add);
  }
}


void filterBaseOnE::removeFiles(const string  &dest, const vector<string> &names)
{
  for (int i=0 ; i<names.size() ; i++ )
  {
    remove ( (dest + names[i]).c_str());
    cout << names[i] << " deleted." << endl;
  }
}

void filterBaseOnE::reverseString (string &seqR , const string &seq , const int &k )
{

  //reversing 
  char* res = new char [k];

  for (int i=k-1 ; i>=0 ; i--)
  {
    res[k-i-1] = seq[i];
  }

  seqR = string (res).substr(0, k);

  delete[] res;
  res = NULL;

  for (int i=0 ; i<seqR.size() ; i++)
  {
    if (seqR[i]=='A')
      seqR[i]='T';
    else if (seqR[i]=='T')
      seqR[i]='A';
    else if (seqR[i] == 'C')
      seqR[i] ='G';
    else if (seqR[i] == 'G')
      seqR[i] ='C';
    else 
      seqR[i] = 'N';
  }
}

void filterBaseOnE::checkKmersAndWriteToFile (const bool &forE0 , const int &num , const string &add, const int &limit)
{

	ifstream file_in;
	ofstream file_out_okay , file_out_nokay ;
	ostringstream convert;
	string line , name ;
	convert << num;
	if (forE0)
	{

		file_in.open ( (add+"e0/reads-count-" + convert.str()+ ".fasta").c_str()  );
		file_out_okay.open ((add+"e0/reads-count-" + convert.str()+ "OK.fasta").c_str() );
		file_out_nokay.open ( (add+"e0/reads-count-" + convert.str()+ "NOK.fasta").c_str() );
	}
	else
	{
		file_in.open ( (add+"e1-2/reads-count-" +convert.str() + ".fasta").c_str()  );
		file_out_okay.open ((add+"e1-2/reads-count-" +convert.str() + "OK.fasta").c_str() );
		file_out_nokay.open ( (add+"e1-2/reads-count-" +convert.str() + "NOK.fasta").c_str()  );
	}


	while (!file_in.eof())
	{
		getline (file_in , line);
		if (line[0]== '>')
		{
			name = line; 
			getline (file_in , line );
			if (checkCountLimit(',' , name , forE0, limit))
			{
				file_out_okay << name << endl << line << endl; 
			}
			else
			{
				file_out_nokay <<  name << endl << line <<  endl ;
			}
		}

	}

	convert.str("");
	file_in.close();
	file_out_okay.close();
	file_out_nokay.close();
}

bool filterBaseOnE::checkCountLimit (const char &c , const string &name , const bool &forE0 , const int &limit)
{
	vector<string> nameSplit;
	int nOfComma , count;
	nOfComma = 1;
	splitWithChar ('|' , name , nameSplit);
	count = atoi (nameSplit[3].c_str());
	if (forE0)
	{// not more than once
		for (int i=0 ; i<nameSplit[0].size() ; i++)
			if ((nameSplit[0])[i] == c )
				nOfComma++;
		if (count != nOfComma)
			return 0;
		return 1;
	}
	else
	{// not more than limit times
		if (count > limit )
			return 0;
		return 1;
	}

}