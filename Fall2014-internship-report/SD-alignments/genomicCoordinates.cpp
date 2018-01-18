#include "genomicCoordinates.h"

genomicCoordinates::genomicCoordinates(string sd_add , string align_both_add)
{
	genomic_SD_add = sd_add;
	align_both_directory_add = align_both_add;

}

void genomicCoordinates::getAlignAttributes (string &align_file, string &chr1, string &chr2 , int &chr1_s,
int &chr1_e, int &chr2_s,int  &chr2_e , string &strand , const string &line_tab)
{
	vector<string> res ;
	splitWithChar ('\t' , line_tab , res);
	chr1 = res[0];
	chr2 = res[6];
	strand = res[5];
	// reported start point in .tab file need to be incremented by one.
	chr1_s =  atoi(res[1].c_str()) + 1;
	chr1_e = atoi (res[2].c_str());
	chr2_s = atoi (res[7].c_str()) + 1;
	chr2_e = atoi (res[8].c_str());
	align_file = align_both_directory_add +res[16];
	cout << "working with file " << res[16] << endl; 
	res.clear();
}
void genomicCoordinates::findSeq ( const string &id ,string &res)
{
	for (int i=0 ; i< id.size() ; i++)
		if (id[i] == 'A' || id[i]=='C' || id[i]=='G' || id[i]=='T' || id[i]=='-')
		{
			res = id.substr (i , id.size()-i+1);
			break;
		}
}

void genomicCoordinates::reportIndel (string &align_pattern,string &first_seq, string &sec_seq,const string &chr1,
const string &chr2,const string &strand,const int &chr1_s,const int &chr2_s,const int &chr1_e , const int &chr2_e)
{
	ofstream report_insert , report_deletion;
	int first_gap_head =0 , sec_gap_head=0;
	report_insert.open ("res/report_insert.txt" , ofstream::app);
	report_deletion.open ("res/report_deletion.txt" , ofstream::app);
	for (int i=0 ; i<align_pattern.size() ; i++ )
	{ 

		if (first_seq[i]=='-')
		{
			first_gap_head++;
		}
		if (sec_seq[i]=='-')
		{
			sec_gap_head++;
		}

		if (align_pattern[i]==' ')
		{
			for (int j=i ; j<align_pattern.size() ; j++)
			{
				if (align_pattern[j]==' ')
					continue;
				else // finding how many consequative indels
				{
					// now the first seq has dels from i to j-1, the second seq has insertions from i to j-1
					if (first_seq[j-1]=='-')
					{
						//reporting deletions for first seq
						first_gap_head--;
						report_deletion << chr1 << "\t" << (chr1_s+(i-1))-first_gap_head << "\t";
						report_deletion << (chr1_s+(i-1))-first_gap_head+ 1 << "\t" << j-i << endl;
						first_gap_head += (j-i);
						// reporting insertions for second seq
						if (strand== "+")
						{
							report_insert<<chr2<<"\t"<<(chr2_s+i)-sec_gap_head<<"\t"<<(chr2_s+j)-sec_gap_head<<endl;
						}
						else
						{
							report_insert<<chr2 <<"\t"<<(chr2_e-(j-1))+sec_gap_head<<"\t"<<(chr2_e-(i-1))+sec_gap_head<< endl;
						}
					}
					if (sec_seq[j-1]=='-')
					{// first seq -> insertion , second deq -> deletions

						report_insert<<chr1<<"\t"<<(chr1_s+i)-first_gap_head<<"\t"<<(chr1_s+j)-first_gap_head<<endl;

						if (strand == "+")
						{ // first seq -> insert and + , sec seq -> del and + 
	
							sec_gap_head--;
							report_deletion << chr2 << "\t" << (chr2_s+(i-1))-sec_gap_head << "\t";
							report_deletion << (chr2_s+(i-1))-sec_gap_head + 1 << "\t" << j-i << endl;
							sec_gap_head += (j-i);
						}
						else
						{ // first seq -> insert and + , sec seq -> del and - 
							sec_gap_head--;
							report_deletion << chr2 << "\t" << (chr2_e-(i))+sec_gap_head << "\t";
							report_deletion << (chr2_e-(i-1))+sec_gap_head << "\t" << j-i << endl;
							sec_gap_head += (j-i);
						}
					}
					i=j-1;
					break;
				}	
			}
		}
		else
			continue;
	}
	align_pattern.clear();
	first_seq.clear();
	sec_seq.clear();
	report_insert.close();
	report_deletion.close();
}

void genomicCoordinates::reportMismatch ( string &align_pattern,string &first_seq, string &sec_seq,const string &chr1,
const string &chr2,const string &strand,const int &chr1_s,const int &chr2_s,const int &chr1_e , const int &chr2_e)
{
	ofstream report;
	int first_gap_head =0 , sec_gap_head=0;
	report.open ("res/report_mismatch.txt" , ofstream::app);

	for (int i=0 ; i<align_pattern.size() ; i++ )
	{ 
		if (align_pattern[i]=='*')
		{
			for (int j=i+1 ; j<align_pattern.size() ; j++)
			{
				if (align_pattern[j]=='*')
					continue;
				else // finding how many consequative mismatches
				{
					
					// now i to j-1 is mistmatch, we should write chr# i j
					report << chr1 << "\t" << (chr1_s+i)-first_gap_head << "\t" << (chr1_s+j)-first_gap_head << endl;
					// for the other chromosome, it depends on strand
					if (strand== "+")
					{
						report << chr2 << "\t" << (chr2_s+i)-sec_gap_head << "\t" << (chr2_s+j)-sec_gap_head << endl;
					}
					else
					{// as j-i * is exits, the ending position can be calculated relatively
						report << chr2 << "\t" << (chr2_e-(j-1))+sec_gap_head << "\t" << (chr2_e-(i-1))+sec_gap_head<< endl;
					}
					i=j-1;
					break;
				}
				
			}

		}
		else if (align_pattern[i]==' ')
		{
			if (first_seq[i]=='-')
				first_gap_head++;
			else
				sec_gap_head++;
		}
		else
			continue;
	}
	//align_pattern.clear();
	//first_seq.clear();
	//sec_seq.clear();
	report.close();
}

void genomicCoordinates::getSDalignInfFromFile (const string &file_name , string &allSeqAlignPatt , string &allFirstSeq, string &allSecSeq)
{
	ifstream file;
	int count , line_find;
	string line , first_seq , align_pat , sec_seq , first_seq_line , sec_seq_line;
	file.open (file_name.c_str());
	count =1 ;
	line_find = 1;
	while (!file.eof())
	{
		if (count == (6*line_find-1) )
		{
			getline (file , first_seq_line);
			findSeq ( first_seq_line , first_seq);

			getline (file , align_pat);

			getline (file , sec_seq_line);
			findSeq ( sec_seq_line , sec_seq);

			if (sec_seq.size() != first_seq.size())
			{// means that one of the sequences start with 'space'
				if (sec_seq.size() > first_seq.size())
				{
					first_seq = first_seq_line.substr (sec_seq_line.size() - sec_seq.size() , sec_seq.size());
				}
				else
				{
					sec_seq = sec_seq_line.substr (first_seq_line.size() - first_seq.size() , first_seq.size());
				}
			}

			// because of many spaces in align pattern, we cant findSeq it,
			// so we substr the sequence needed as we know the length of it from the tail
			allSeqAlignPatt += align_pat.substr ( align_pat.size()-first_seq.size() , first_seq.size());
			allFirstSeq += first_seq;
			allSecSeq += sec_seq;
			line_find++;
			count += 3;
		}
		else 
		{
			getline (file , line);
			count++;
		}
	
	}
	file.close();
}

void genomicCoordinates::splitWithChar (const char &identifier , const string &id ,vector<string> &res)
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

void genomicCoordinates::reportMismatchIndelCoordinatesOfSDalignments()
{
	int chr1_s , chr1_e , chr2_e, chr2_s;
	string allSeqAlignPatt ,allFirstSeq, allSecSeq, line_tab , chr1 , strand , chr2 , align_file;
	ifstream alignInf_file;
	alignInf_file.open (genomic_SD_add.c_str());
	system("mkdir \"res\"");
	remove ("res/report_mismatch.txt");
	remove ("res/report_deletion.txt");
	remove ("res/report_insert.txt");
	int count = 1;
	while (!alignInf_file.eof())
	{
		getline (alignInf_file , line_tab);
		if (count == 1) 
			{
				count++;
				continue;
			}
		

		if (count%2==0)
		{



			cout << "processing" << endl;
			getAlignAttributes (align_file, chr1, chr2 , chr1_s,chr1_e,chr2_s, chr2_e , strand , line_tab);
			//CHANGE
			

/*			chr1 = "chr1";
			chr2 = "chr2";
			chr2_s = 1;
			chr2_e = 1000;
			chr1_s = 1;
			chr1_e = 1100;
			strand = "-";
			align_file = "test2.txt";
			//CHANGE
*/			
			getSDalignInfFromFile (align_file ,  allSeqAlignPatt , allFirstSeq , allSecSeq);
			reportMismatch (allSeqAlignPatt,allFirstSeq,allSecSeq,chr1,chr2,strand,chr1_s,chr2_s,chr1_e,chr2_e);
			reportIndel (allSeqAlignPatt,allFirstSeq,allSecSeq,chr1,chr2,strand,chr1_s,chr2_s,chr1_e,chr2_e);
			//because two continues lines are giving us the same information 
			getline ( alignInf_file , line_tab);
			count+=2;
		}
	}
	
}