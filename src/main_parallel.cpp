#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <ctime>
#include <vector>
#include <stdlib.h>
#include <typeinfo>
#include <string>
//#include "hed.h"

#define NDEBUG
#include <cassert>

#define Pi 3.14159265
#define pi 3.14156

using namespace std;


struct trip{

    int first;
    int second;
    int third;

};

pair<int,int> count21=make_pair(0,0), count31=make_pair(0,0), count32=make_pair(0,0);

vector < vector <float> > corrected_rmaps;

vector < int > rmaps_to_correct;

vector< vector<float> > rmap;

vector< vector<int> > rmap_quantized;

vector < vector <int> > rmap_store;

vector < set <int> > rmap_relations;

vector < vector <int> > rmap_relations_filtered;

vector < vector <int> > rmap_relations_final;

vector <vector < vector < pair <int, int> > > > rmap_alignments;

vector <vector < vector < trip > > > multi_align_grid;

vector< string > rmap_name;

vector < vector < pair <int, int> > > kmer_store;

vector < vector < string > > concensus_res;

int buket_size=6;

map < string, int > kmer_map;

map < int, string > reverse_map;

vector < pair <int, int> > alignment;

//vector < pair <int, int> > ali;

int max_alignments=100;

int k_size=4;

int m_size=2;

int min_frags=15;

int align_size=3;

float threshold=0.85;

int min_alignments = 5;

int min_concensus = 3;

float bias_for_single=0.1;

int start_r,end_r,diff;

string nu;



float find_score( float a, float b){

    if (b>a)
        return(a/b);

    else
        return(b/a);
}

struct triplet{

    int pos1;
    int rmap_2;
    int pos2;

};


vector< vector < triplet > > complete_relations;

vector < pair<int ,int> > optimized_overlap_alignment(vector< float >& rmap1, vector< float >& rmap2, int p1, int p2);


float unit_bias=0.2;

trip make_trip(int a, int b, int c){

    trip t;
    t.first=a;
    t.second=b;
    t.third=c;
    return(t);
}

int check_relation(int rmap1, int pos1, int rmap2, int pos2){

//    return 1;

    int flag=1;

    //cout << endl << endl;

    //cout << "Aligning " << rmap1 << " with " << rmap2 <<endl ;
    float total_sc=0;

    while((pos1<rmap[rmap1].size()-1) && (pos2<rmap[rmap2].size()-1)){

        float sc=0;

        float score=0, new_score=0;

        int cp1=pos1, cp2=pos2;

        score= find_score(rmap[rmap1][pos1],rmap[rmap2][pos2])+bias_for_single; // 1-1
        sc=score/2;

//        if(cp1+1<rmap[rmap1].size() && cp2+1<rmap[rmap2].size())
//            score= score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][cp1+1],rmap[rmap2][cp2+1]);


        //cout << "compare " << rmap[rmap1][pos1] << "  AND  " << rmap[rmap2][pos2] << " Score : " << score << endl;
        if(pos1+1<rmap[rmap1].size()){

            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1],rmap[rmap2][pos2]); //2-1

            if(pos1+2<rmap[rmap1].size() && pos2+1<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+2],rmap[rmap2][pos2+1]);

            //cout << "compare " << rmap[rmap1][pos1]<< "+"<< rmap[rmap1][pos1+1] << "  AND  " << rmap[rmap2][pos2] << " Score : " << new_score<< endl;

            if (new_score>score){
                    score=new_score;
                    cp1=pos1+1;
                    cp2=pos2;

                    sc=score/3;


            }
        }


        if(pos1+2<rmap[rmap1].size())
        new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2],rmap[rmap2][pos2]); //3-1

        if(pos1+3<rmap[rmap1].size() && pos2+1<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+3],rmap[rmap2][pos2+1]);

        //cout << "compare " << rmap[rmap1][pos1]<< "+"<< rmap[rmap1][pos1+1] << "+"<< rmap[rmap1][pos1+2] << "  AND  " << rmap[rmap2][pos2] << " Score : " << new_score<< endl;


        if (new_score>score){
                score=new_score;
                cp1=pos1+2;
                cp2=pos2;
                sc=score/4;
        }

        if(pos2+1<rmap[rmap2].size())
        new_score= find_score(rmap[rmap1][pos1],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]); //1-2

        if(pos1+1<rmap[rmap1].size() && pos2+2<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+1],rmap[rmap2][pos2+2]);

        if (new_score>score){
                score=new_score;
                cp1=pos1;
                cp2=pos2+1;
                sc=score/3;
        }

        if(pos2+2<rmap[rmap2].size())
        new_score= find_score(rmap[rmap1][pos1],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]); //1-3

        if(pos1+1<rmap[rmap1].size() && pos2+3<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+1],rmap[rmap2][pos2+3]);

        if (new_score>score){
                score=new_score;
                cp1=pos1;
                cp2=pos2+2;
                sc=score/4;
        }

        if((pos1+1<rmap[rmap1].size()) && (pos2+1<rmap[rmap2].size()))
        new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]);//

        if(pos1+2<rmap[rmap1].size() && pos2+2<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+2],rmap[rmap2][pos2+2]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+1;
                    cp2=pos2+1;
                    sc=score/4;
        }

        if((pos1+1<rmap[rmap1].size()) && (pos2+2<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]);

        if(pos1+2<rmap[rmap1].size() && pos2+3<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+2],rmap[rmap2][pos2+3]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+1;
                    cp2=pos2+2;
                    sc=score/5;
        }

        if((pos1+2<rmap[rmap1].size()) && (pos2+1<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]);

        if(pos1+3<rmap[rmap1].size() && pos2+2<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+3],rmap[rmap2][pos2+2]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+2;
                    cp2=pos2+1;
                    sc=score/5;
        }


        if((pos1+2<rmap[rmap1].size()) && (pos2+2<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]);

        if(pos1+3<rmap[rmap1].size() && pos2+3<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+3],rmap[rmap2][pos2+3]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+2;
                    cp2=pos2+2;
                    sc=score/6;
        }

        if((pos1+3<rmap[rmap1].size()) && (pos2+2<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2]+rmap[rmap1][pos1+3],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]);

        if(pos1+4<rmap[rmap1].size() && pos2+3<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+4],rmap[rmap2][pos2+3]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+3;
                    cp2=pos2+2;
                    sc=score/7;
        }

        if((pos1+2<rmap[rmap1].size()) && (pos2+3<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]+rmap[rmap2][pos2+3]);

        if(pos1+3<rmap[rmap1].size() && pos2+4<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+3],rmap[rmap2][pos2+4]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+2;
                    cp2=pos2+3;
                    sc=score/7;
        }



        if (score>threshold){
            total_sc+=sc;

//            for(int i=pos1; i<=cp1; i++)
//                cout << i << " (" << rmap[rmap1][i] << ")   ";
//
//            cout << " aligns with ";
//            for(int i=pos2; i<=cp2; i++)
//                cout << i << " (" << rmap[rmap2][i] << ")   ";
//
//            cout << endl;
        }

        else{

//            cout << "Maps do not align ";
            flag=0;

            break;

        }

        pos1=cp1+1;
        pos2=cp2+1;

    }
    if(total_sc<2)
        flag=0;
    //if((pos1>=rmap[rmap1].size()) || (pos2>=rmap[rmap2].size()))
        //cout <<endl << " Alignment successful";
    //cout << "  FLAG: " << flag;

    if(flag==1){
        //cout<< endl << "Alignment success" << endl;
    }
        return(flag);
}



void read_data( char* filename){

    int rmap_count=0;

    cout<< "Reading rmaps....."<<endl;

    ifstream file(filename, ifstream::in);
	if (!file.is_open()) {
		cout << "Error opening file" << endl;
		exit(1);
	}

	int total_rmaps=0;
//	cout << "Number of Rmaps ? ";

	//cin >> total_rmaps;


	string str,str2;

	while(getline(file, str))
    {
        istringstream ss(str);

        vector < float > v1;
        vector < int > v2;

        ss >> str2;

        total_rmaps++;

        //cout << rmap_count << " ";

//        float first_num;
//        ss >> first_num;

        float num;
        while(ss >> num)
        {
        //cout << num << "    ";
         v1.push_back(num);
         v2.push_back(round(num/buket_size));
        }

        if(v1.size()<min_frags)
            continue;

        rmap_name.push_back(str2);
        rmap.push_back(v1);
        rmap_quantized.push_back(v2);
        rmap_count++;
    }
    //cout << endl;
    file.close();

    cout << "Total rmaps : " <<total_rmaps<<endl<<"Rmaps with 10+ fragments : "<< rmap_count<<endl;
    total_rmaps=rmap_count;

    rmap_store.resize(total_rmaps);

//	rmap_relations.resize(total_rmaps);
//	rmap_relations_filtered.resize(total_rmaps);
	rmap_relations_final.resize(total_rmaps);
	rmap_alignments.resize(total_rmaps);
	multi_align_grid.resize(total_rmaps);
	concensus_res.resize(total_rmaps);

//	complete_relations.resize(total_rmaps);

}


void find_kmers( int k ){

    int lastmapped=0;

    for( int i=0; i<rmap_quantized.size(); i++){


        for( int j=0; j<rmap_quantized[i].size()-1-k+1; j++){

            string kmer;
            ostringstream oss;
            for (int l=0; l<k; l++){
                oss << rmap_quantized[i][j+l] << "_";
            }


            istringstream iss(oss.str());
            iss >> kmer;

            if(kmer_map.find(kmer) == kmer_map.end()){
                kmer_map[kmer]=lastmapped;
                reverse_map[lastmapped]=kmer;
                lastmapped++;

            }



            kmer_store.resize(lastmapped);

            int value=kmer_map[kmer];
            rmap_store[i].push_back(value);

            pair < int, int> rmap_position;

            rmap_position=make_pair(i,j);

            kmer_store[value].push_back(rmap_position);


        }
        //cout << endl;
    }

}



void check_all_relations(){

   	vector < triplet >::iterator it_triplet;

//   	ofstream out;
//
//   	ofstream out2;
//
   	ofstream out3;
//
//   	ofstream out_data;

   	//out_data.open("see_alignments.txt");

//   	out.open("better_relations.txt");
//
//   	out2.open("review.txt");
   	out3.open("filter_stats.txt");

   	for (int i=0; i<rmap.size(); i++){
 //   for (int i=0; i<complete_relations.size(); i++){

//        out<<rmap_name[i]<< "    ";
        int first_fil_cnt=0;

        for (it_triplet=complete_relations[i].begin(); it_triplet!=complete_relations[i].end(); ++it_triplet){
            if(check_relation(i,it_triplet->pos1,it_triplet->rmap_2,it_triplet->pos2)){
//              out<<it_triplet->rmap_2<<"    ";

                first_fil_cnt++;
              alignment = optimized_overlap_alignment(rmap[i],rmap[it_triplet->rmap_2],it_triplet->pos1,it_triplet->pos2);

//              cout<<alignment.size()<<" "<<sizeof(alignment)<<endl;

              if(alignment[0].first>15 && alignment[0].second>0){
//                out<< it_triplet->rmap_2<< "   ";
                rmap_relations_final[i].push_back(it_triplet->rmap_2);
                rmap_alignments[i].push_back(alignment);
//                ali.clear();
//                out2<< " reference rmap : " << i << " target rmap : " << it_triplet->rmap_2 <<endl;
                for(int j=1; j<alignment.size();j++){
//                    out2 << "   ("<<alignment[j].first<<","<<alignment[j].second<<")";
//                    ali.push_back(make_pair(alignment[j].first,alignment[j].second));
                }

//                rmap_alignments[i].push_back(ali);
//                out2 << endl;
              }


//              rmap_relations_filtered[i].push_back(it_triplet->rmap_2);

            }

            }
//        cout << i << "  " << rmap_relations_filtered[i].size()<<endl;
//        out<< endl;

        out3 << rmap_alignments[i].size()<< " / "<<first_fil_cnt<<" /" <<complete_relations[i].size() << endl;
//        out2 << endl;

        multi_align_grid[i].resize(rmap_alignments[i].size());

        for(int k=0; k<multi_align_grid[i].size(); k++){
            multi_align_grid[i][k].resize(rmap[i].size(),make_trip(0,0,0));

        }

        if(i%100000==0)
            cout << endl<<i<<" "<<rmap_alignments[i].size()<< " / "<<first_fil_cnt<<" /" <<complete_relations[i].size();

        }

//    out.close();
//    out2.close();
    out3.close();
//    out_data.close();
//




}

void build_multi_align_grid(){


    long mag=0;
    for(int i=0; i<rmap_quantized.size(); i++){


        for(int j=0; j<rmap_alignments[i].size(); j++){

            mag+=rmap_alignments[i][j].size();
            for(int k=1; k<rmap_alignments[i][j].size()-1; k++){

                int ref_in=rmap_alignments[i][j][k].first;
                int tar_in=rmap_alignments[i][j][k].second;

//                if(k==rmap_alignments[i][j].size()-1){
//                    multi_align_grid[i][j][ref_in]=make_pair(ref_gap,tar_gap);
//
//                }

                int ref_gap=rmap_alignments[i][j][k+1].first-rmap_alignments[i][j][k].first;
                int tar_gap=rmap_alignments[i][j][k+1].second-rmap_alignments[i][j][k].second;

                while(ref_in<rmap_alignments[i][j][k+1].first){
                    multi_align_grid[i][j][ref_in]=make_trip(ref_gap,tar_gap,tar_in);
                    ref_in++;
                }


                }

            }

    }

}

map < string, int > mapit;
map < int, string > rev_mapit;
vector < int > keep_count;

void find_concensus(char* ch){



    int lastmapped=0;

    ofstream out;

    string st;
    string str(ch);

    st = "results/concensus"+str+".txt";
    out.open(st.c_str());

    for(int i=0; i<rmap_quantized.size(); i++){



    if(multi_align_grid[i].size()<min_alignments)
        continue;

    rmaps_to_correct.push_back(i);
    out<<rmap_name[i];

//    cout << multi_align_grid[i][0].size()<<endl;

        for(int k=0; k<multi_align_grid[i][0].size(); k++){


            fill(keep_count.begin(), keep_count.end(), 0);

            for(int j=0; j<multi_align_grid[i].size(); j++){

                if(multi_align_grid[i][j][k].first==0)
                    continue;

                string str;
                ostringstream oss;
                oss<<multi_align_grid[i][j][k].first<<","<<multi_align_grid[i][j][k].second;

                istringstream iss(oss.str());
                iss >> str;

                if(mapit.find(str)==mapit.end()){
                    mapit[str]=lastmapped;
                    rev_mapit[lastmapped]=str;
                    lastmapped++;
                    keep_count.resize(mapit.size());
                    keep_count[mapit[str]]++;
//                    keep_count.push_back(1);

                }
                else{
                    keep_count[mapit[str]]++;
                }

            }
            int maxel=0;
            for(int x=0; x<keep_count.size(); x++){
                if(keep_count[x]>keep_count[maxel])
                    maxel=x;
            }

            if(keep_count[maxel]<min_concensus){
                concensus_res[i].push_back("0,0");
                out << " (0,0) ";
            }
            else{
                concensus_res[i].push_back(rev_mapit[maxel]);
                out << " (" <<rev_mapit[maxel]<<") ";
            }



        }
        out<<endl;

    }

out.close();

}


bool check(string & str, int x, int y){

                char chara=str[0];
                char charb=str[2];

                int a=chara - '0';
                int b=charb - '0';

                if(a==x && b==y)
                    return true;
                else
                    return false;
}

void correct_rmaps(char* ch){

    ofstream out;

    ofstream out2;

    long tot_frag=0;
    long untouched=0;

    string st;

    string in(ch);

    st = "results/corrections"+in+".txt";

    out2.open(st.c_str());

    st=  "results/corrected_rmaps"+in+".txt";
    out.open(st.c_str());

    corrected_rmaps.resize(rmaps_to_correct.size());
    for(int i=0; i<rmaps_to_correct.size(); i++){

        if(i%1000==0){
            cout<<i<<endl;
        }

        int add_err=0;
        int del_err=0;
        int this_rmap=rmaps_to_correct[i];
        out<<rmap_name[this_rmap];
        out<<"\n enzyme enzyme";

        bool set_bit2=false;
        bool set_bit3=false;
        bool set_bit4=false;
        float counter;

        tot_frag+=concensus_res[this_rmap].size();

        for(int j=0; j<concensus_res[this_rmap].size(); j++){

            if (set_bit2==true){
                    set_bit2=false;
                    continue;
            }

            else if(set_bit3==true){
                set_bit3=false;
                set_bit2=true;
                continue;
            }

            else if(set_bit4==true){
                set_bit4=false;
                set_bit3=true;
                continue;
            }

            string str = concensus_res[this_rmap][j];

            char chara=str[0];
            char charb=str[2];

            int a=chara - '0';
            int b=charb - '0';

            if(a==1 && b==1){
                float tot11=0;
                counter=0;
                for(int k=0; k<multi_align_grid[this_rmap].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot11+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                    }
                }
                tot11+=rmap[this_rmap][j];
                counter++;

                corrected_rmaps[i].push_back(tot11/counter);
//                out<<"  "<<tot11/counter<<"("<<j<<","<<str<<"correction)";
                out<<"  "<<roundf((tot11/counter) * 1000) / 1000;


            }

            else if(a==1 && b==2){
                del_err++;
                float tot121=0;
                float tot122=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot121+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot122+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                    }

                }
                corrected_rmaps[i].push_back(tot121/counter);
                corrected_rmaps[i].push_back(tot122/counter);
//                out << "  "<<tot121/counter<<"  "<<tot122/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot121/counter) * 1000) / 1000<<"  "<<roundf((tot122/counter) * 1000) / 1000;
            }

            else if(a==2 && b==1){


                count21.first++;
                if(!check(concensus_res[this_rmap][j+1],2,1)){
                    count21.second++;
                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                    out << "  "<<rmap[this_rmap][j];
                    untouched++;
                    continue;
                }
                add_err++;
                float tot21=0;
                set_bit2=true;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot21+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                    }

                }
                corrected_rmaps[i].push_back(tot21/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot21/counter) * 1000) / 1000;

            }

            else if(a==2 && b==3){

            del_err++;

//                count21.first++;
//                if(!check(concensus_res[this_rmap][j+1],2,1)){
//                    count21.second++;
//                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
//                    out << "  "<<rmap[this_rmap][j];
//                    continue;
//                }
                float tot231=0;
                float tot232=0;
                float tot233=0;
                counter=0;
                set_bit2=true;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot231+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot232+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot233+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                    }

                }
                corrected_rmaps[i].push_back(tot231/counter);
                corrected_rmaps[i].push_back(tot232/counter);
                corrected_rmaps[i].push_back(tot233/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot231/counter) * 1000) / 1000<< "  "<<roundf((tot232/counter) * 1000) / 1000<< "  "<<roundf((tot233/counter) * 1000) / 1000;

            }

             else if(a==2 && b==4){
                del_err++;
                del_err++;

//                count21.first++;
//                if(!check(concensus_res[this_rmap][j+1],2,1)){
//                    count21.second++;
//                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
//                    out << "  "<<rmap[this_rmap][j];
//                    continue;
//                }
                float tot241=0;
                float tot242=0;
                float tot243=0;
                float tot244=0;
                set_bit2=true;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot241+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot242+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot243+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                        tot244+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+3];
                    }

                }
                corrected_rmaps[i].push_back(tot241/counter);
                corrected_rmaps[i].push_back(tot242/counter);
                corrected_rmaps[i].push_back(tot243/counter);
                corrected_rmaps[i].push_back(tot244/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot241/counter) * 1000) / 1000<< "  "<<roundf((tot242/counter) * 1000) / 1000<< "  "<<roundf((tot243/counter) * 1000) / 1000<< "  "<<roundf((tot244/counter) * 1000) / 1000;

            }

            else if(a==2 && b==2){

                if(!check(concensus_res[this_rmap][j+1],2,2)){
                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                    out << "  "<<rmap[this_rmap][j];
                    continue;
                }
                set_bit2=true;
                float tot221=0;
                float tot222=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot221+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot222+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                    }

                }
                corrected_rmaps[i].push_back(tot221/counter);
                corrected_rmaps[i].push_back(tot222/counter);
//                out << "  "<<tot221/counter<<"  "<<tot222/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot221/counter) * 1000) / 1000<<"  "<<roundf((tot222/counter) * 1000) / 1000;

            }

            else if(a==1 && b==3){
                del_err++;
                del_err++;
                float tot131=0;
                float tot132=0;
                float tot133=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot131+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot132+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot133+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                    }

                }
                corrected_rmaps[i].push_back(tot131/counter);
                corrected_rmaps[i].push_back(tot132/counter);
                corrected_rmaps[i].push_back(tot133/counter);
//                out << "  "<<tot121/counter<<"  "<<tot122/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot131/counter) * 1000) / 1000<<"  "<<roundf((tot132/counter) * 1000) / 1000<< "  " << roundf((tot133/counter) * 1000) / 1000;
            }

            else if(a==1 && b==4){
                del_err++;
                del_err++;
                del_err++;
                float tot141=0;
                float tot142=0;
                float tot143=0;
                float tot144=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot141+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot142+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot143+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                        tot144+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+3];
                    }

                }
                corrected_rmaps[i].push_back(tot141/counter);
                corrected_rmaps[i].push_back(tot142/counter);
                corrected_rmaps[i].push_back(tot143/counter);
                corrected_rmaps[i].push_back(tot144/counter);
//                out << "  "<<tot121/counter<<"  "<<tot122/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot141/counter) * 1000) / 1000<<"  "<<roundf((tot142/counter) * 1000) / 1000<< "  " << roundf((tot143/counter) * 1000) / 1000<< "  " << roundf((tot144/counter) * 1000) / 1000;
            }

            else if(a==3 && b==3){


                set_bit3=true;
                float tot331=0;
                float tot332=0;
                float tot333=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot331+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot332+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot333+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                    }

                }
                corrected_rmaps[i].push_back(tot331/counter);
                corrected_rmaps[i].push_back(tot332/counter);
                corrected_rmaps[i].push_back(tot333/counter);
//                out << "  "<<tot331/counter<<"  "<<tot332/counter<<"  "<<tot333/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot331/counter) * 1000) / 1000<<"  "<<roundf((tot332/counter) * 1000) / 1000<<"  "<<roundf((tot333/counter) * 1000) / 1000;

            }

            else if(a==3 && b==1){


                count31.first++;
                if(!(check(concensus_res[this_rmap][j+1],3,1))){
                    count31.second++;
                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                    out << "  "<<rmap[this_rmap][j];
                    untouched++;
                    continue;
                }
                add_err++;
                add_err++;

                set_bit3=true;
                float tot31=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot31+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                    }

                }
                corrected_rmaps[i].push_back(tot31/counter);
//                out << "  "<<tot321/counter<<"  "<<tot322/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot31/counter) * 1000) / 1000;

            }

            else if(a==3 && b==2){


                count32.first++;
                if(!(check(concensus_res[this_rmap][j+1],3,2))){
                    count32.second++;
                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                    out << "  "<<rmap[this_rmap][j];
                    untouched++;
                    continue;
                }
                add_err++;
                set_bit3=true;
                float tot321=0;
                float tot322=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot321+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot322+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                    }

                }
                corrected_rmaps[i].push_back(tot321/counter);
                corrected_rmaps[i].push_back(tot322/counter);
//                out << "  "<<tot321/counter<<"  "<<tot322/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot321/counter) * 1000) / 1000<<"  "<<roundf((tot322/counter) * 1000) / 1000;

            }

            else if(a==4 && b==1){

            add_err++;
            add_err++;
            add_err++;
//                count21.first++;
//                if(!check(concensus_res[this_rmap][j+1],2,1)){
//                    count21.second++;
//                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
//                    out << "  "<<rmap[this_rmap][j];
//                    continue;
//                }
                float tot41=0;
                counter=0;
                set_bit4=true;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot41+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];

                    }

                }
                corrected_rmaps[i].push_back(tot41/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot41/counter) * 1000) / 1000;

            }

            else if(a==4 && b==3){

            add_err++;
//                count21.first++;
//                if(!check(concensus_res[this_rmap][j+1],2,1)){
//                    count21.second++;
//                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
//                    out << "  "<<rmap[this_rmap][j];
//                    continue;
//                }
                float tot431=0;
                float tot432=0;
                float tot433=0;
                counter=0;
                set_bit4=true;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot431+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot432+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot433+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                    }

                }
                corrected_rmaps[i].push_back(tot431/counter);
                corrected_rmaps[i].push_back(tot432/counter);
                corrected_rmaps[i].push_back(tot433/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot431/counter) * 1000) / 1000<< "  "<<roundf((tot432/counter) * 1000) / 1000<< "  "<<roundf((tot433/counter) * 100) / 100;

            }



            else if(a==0 && b==0){
                corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                out << "  "<<rmap[this_rmap][j];
                untouched++;
            }

            else{
                corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                out << "  "<<rmap[this_rmap][j];
            }


        }

        out2<< rmap_name[this_rmap]<<"  "<<add_err<<"   "<< del_err<<endl;

        out<<"\n \n";

    }

    out.close();
    out2.close();

}



int main ( int argc, char *argv[] ) {


//    int num=atoi(argv[1]);
    if(argc<2){
        cout<<"Usage: "<<argv[0]<<" <Rmap_file> <Number_of_parallel_streams> <stream_index> \n For example, if total number of Rmaps is 1,000,000 and total number of streams is 1000. Then index=1 will error correct Rmaps 1 to 1000. index=2 will error correct Rmaps 1001 to 2000 and so on."<<endl;
        return 1;
    }
    int tot_streams=atoi(argv[2]);
    int index=atoi(argv[3]);

    time_t ostart = clock();

	read_data(argv[1]);

	time_t reading_t_end = clock();

	time_t kmer_finding_start = clock();
	cout << "Indexing kmers....."<<endl;
	find_kmers(k_size);
	time_t kmer_finding_end = clock();



	//print_kmer_store();

	//print_rmaps();
	//print_rmaps_quantized();
	//print_rmap_names();

	double t1=0;
	double t2=0;

	for (int i=0; i<kmer_store.size();i++){
        t1++;
        if(kmer_store[i].size()>1){
            t2++;
//            outfile<<reverse_map[i]<<"     "<< kmer_store[i].size()<<endl;
        }
//        if(i>=kmer_store.size()-1)
//            cout <<endl <<"last kmer reached." <<i<< endl;
	}

//	cout<<"Percentage of orphan kmers : "<<(1-t2/t1);

    float orphan_count=0;

//    ofstream rmap_file;

    cout<< "Finding related Rmaps....."<<endl;

    long cl=0;

    double timer_total=0;

    time_t checking_time_begin = clock();

    float tot_score=0;

    int tot_fil2=0;

    int start=rmap_quantized.size()*(index-1)/tot_streams;
    int fin=rmap_quantized.size()*(index)/tot_streams;

    cout<<"\n Rmap start index:"<< start<<" Rmap end index:"<<fin<<"\n";




    for (int i=start; i<fin; i++){

        vector<triplet> position_index;
        vector<int> rmap_count(rmap_name.size(),0);

        position_index.resize(rmap_name.size());

        time_t ts=clock();

        int fil_cnt1=0;
        int fil_cnt2=0;




        float avg_score=0;

        for (int j=0; j<rmap_store[i].size(); j++){

            for (int m=0; m<kmer_store[rmap_store[i][j]].size(); m++){

                if((kmer_store[rmap_store[i][j]][m].first)!=i)
                    rmap_count[kmer_store[rmap_store[i][j]][m].first]++;

                if(rmap_count[kmer_store[rmap_store[i][j]][m].first]==1){

                    triplet t;
                    t.pos1=j;
                    t.rmap_2=kmer_store[rmap_store[i][j]][m].first;
                    t.pos2=kmer_store[rmap_store[i][j]][m].second;

                    position_index[t.rmap_2]=t;
                }


                if(rmap_count[kmer_store[rmap_store[i][j]][m].first]==m_size){

                    fil_cnt1++;



//                    complete_relations[i].push_back(position_index[kmer_store[rmap_store[i][j]][m].first]);

                    triplet tem1=position_index[kmer_store[rmap_store[i][j]][m].first];

//                    related1<<" "<<tem1.rmap_2;

                    if(check_relation(i,tem1.pos1,tem1.rmap_2,tem1.pos2)){

//                        related2<<" "<<tem1.rmap_2;
                        fil_cnt2++;
                        alignment = optimized_overlap_alignment(rmap[i],rmap[tem1.rmap_2],tem1.pos1,tem1.pos2);

                        if(alignment[0].first>15 && alignment[0].second>0){
                            rmap_relations_final[i].push_back(tem1.rmap_2);
                            rmap_alignments[i].push_back(alignment);
                            cl+=rmap_alignments[i].size();
//                            related<<" "<<tem1.rmap_2;
                            avg_score+=alignment[0].first;
                        }
                    }
                 }

                 if(rmap_relations_final[i].size()==max_alignments)
                    break;

            }
            if(rmap_relations_final[i].size()==max_alignments)
                    break;
        }

        if(rmap_relations_final[i].size()>0)
            avg_score=avg_score/rmap_relations_final[i].size();

        tot_score+=avg_score;

        tot_fil2+=fil_cnt2;

//        related1<<endl;
//        related2<<endl;
//
//        related<<endl;
        multi_align_grid[i].resize(rmap_alignments[i].size());

        for(int k=0; k<multi_align_grid[i].size(); k++){
                multi_align_grid[i][k].resize(rmap[i].size(),make_trip(0,0,0));

        }

        time_t te=clock();

        timer_total+=double( te - ts) / CLOCKS_PER_SEC;

        float rmaps_left=rmap_quantized.size()-i-1;

        if(i%100==0)
                cout<<"Rmaps completed:"<<i+1<<" Time per Rmap: "<<timer_total/(i+1-start_r)<< " Approximate time remaining:"<< rmaps_left*(timer_total/(i+1-start_r))<<endl;



    }


    time_t checking_time_end = clock();

    cout <<endl<< "Multi alignment grid building begins....."<<endl;

    build_multi_align_grid();

    time_t mag_end = clock();

//    print_multi_align_grid();
    cout << "Finding consensus....."<<endl;
    find_concensus(argv[2]);

    time_t find_con_end = clock();

    cout<<endl<<"Rmap correction begins.....";
//    rmaps_to_corr();

    correct_rmaps(argv[2]);

    time_t corr_rmaps_end = clock();

    time_t oend = clock();
    cout << endl << "total running time:	" << double(oend - ostart) / CLOCKS_PER_SEC << endl;

}


