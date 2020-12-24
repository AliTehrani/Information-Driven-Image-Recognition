/*
 *  tGame.cpp
 *  HMMBrain
 *
 *  Created by Arend on 9/23/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "tGame.h"
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

#define activeFinalizeDecision false

//int auditIndex=0;
//int outputIndex=brainSize-1;
//int pos=0;
int windowSize = 2;

// read sequence of saccade positions from file
vector<pair<int, int>> readSaccadePositions(string fn){
    vector<pair<int, int>> xys;
    ifstream infile(fn);
    string line;
    while(getline(infile, line)){
	istringstream iss(line);
	int x,y;
	iss >> x >> y;
	xys.push_back(pair<int, int>{x, y});
    }
    return xys;
}

matrix tGame::readPix(ifstream &infile){
    matrix pix;

    string line;
    int p;
    getline(infile, line);
    while(line!=""){
      istringstream l(line);
      vector<int> row;
      while(l >> p){
	row.push_back(p);
      }
      pix.push_back(row);

      getline(infile, line);
    }
    return pix;
}

tGame::tGame(string mnistFn, string saccadeFn){//double rewPd, double rewSt, double rewNd){
    //readImageData(mnist, trials);

    // read x,y positions of the saccades
    xy = readSaccadePositions(saccadeFn);

    // read image data

    ifstream infile(mnistFn);
    
    int digit, rep;
    string line;
    while(getline(infile,line)){
        // read digit and rep
        auto pos=line.find('-');
	if (pos!=string::npos){
	  digit=stod(line.substr(0, pos));
	  rep=stod(line.substr(pos+1));
	  //cout << digit << '-' << rep << endl;
	}

	// read pixels
	auto image=readPix(infile);
    
	pair<int,int> key{digit,rep};

	vector<vector<int>> seq;
	for (auto p: xy){
	    int xInit=p.first*windowSize, xFin=(p.first+1)*windowSize, yInit=p.second*windowSize, yFin=(p.second+1)*windowSize;
	    vector<int> subimage;
	    for (int i=xInit; i<xFin; i++){
		for (int j=yInit; j<yFin; j++){
		    subimage.push_back(image[i][j]);
		//cout << trial->image[i][j]<< " ";
		}
	    }
	    seq.push_back(subimage);
	    //cout << endl;
	}
	sequences[key] = seq;

	auto trial = new tTrial(digit, rep, vector<vector<int>>{});//image);
	/*
	if (trials.find(key)!=trials.end())
	  cout << "image already exists" << endl;
	
	trials[key]=image;
	*/
	digitTrialIdx[digit].push_back( trials.size() );
	trials.push_back(trial);
    }
    datasetSize = trials.size();
}

tGame::~tGame(){
    for (int i=0; i<trials.size(); i++)
	delete trials[i];
    trials.clear();
}


void tGame::evaluateAgent(tAgent* agent, tTrial* trial) const{
    int outputSum=agent->getAction();
    auto outputs = agent->getOutputs();

    //cout << outputs.size() << '\t' << trial->digit << endl;
    //copy(outputs.begin(), outputs.end(), ostream_iterator<int>(cout,","));
    //int q; cin >> q;
    if (outputs[trial->digit]==1){
	//cout << outputs[trial->digit] << endl;
	agent->fitness += trial->rew/double(outputSum);
	if (outputSum==1){
	    agent->correct[trial->digit]++;
	}
    }
    
    //agent->fitness -= (outputSum-1)*trial->penalty;

}

string tGame::executeGame(tAgent* agent, bool test, bool report, tTrial* trial) const{
    //cout << trial->digit << '\t' << trial->rep << endl;
    agent->resetBrain();
    //for (auto p: xy){
    pair<int,int> key{trial->digit, trial->rep};
    for (const auto pix: sequences.at(key)){
	/*
	int xInit=p.first*windowSize, xFin=(p.first+1)*windowSize, yInit=p.second*windowSize, yFin=(p.second+1)*windowSize;
	vector<int> pix;
	for (int i=xInit; i<xFin; i++){
	    for (int j=yInit; j<yFin; j++){
		pix.push_back(trial->image[i][j]);
		//cout << trial->image[i][j]<< " ";
	    }
	    //cout << endl;
	}
	//cout << endl;
	*/
	
	agent->setSensors(pix);
	for (int bu=0; bu<agent->brainUpdates; bu++)
	    agent->updateStates();
	
	if (activeFinalizeDecision){
	    if (agent->isDecisionFinalized()){
		evaluateAgent(agent, trial);
		break;
	    }
	}
    }

    if (!activeFinalizeDecision)
	evaluateAgent(agent, trial);

    return "";
}
/*  
double tGame::mutualInformation(vector<int> A,vector<int>B){
	set<int> nrA,nrB;
	set<int>::iterator aI,bI;
	map<int,map<int,double> > pXY;
	map<int,double> pX,pY;
	int i,j;
	double c=1.0/(double)A.size();
	double I=0.0;
	for(i=0;i<A.size();i++){
		nrA.insert(A[i]);
		nrB.insert(B[i]);
		pX[A[i]]=0.0;
		pY[B[i]]=0.0;
	}
	for(aI=nrA.begin();aI!=nrA.end();aI++)
		for(bI=nrB.begin();bI!=nrB.end();bI++){
			pXY[*aI][*bI]=0.0;
		}
	for(i=0;i<A.size();i++){
		pXY[A[i]][B[i]]+=c;
		pX[A[i]]+=c;
		pY[B[i]]+=c;
	}
	for(aI=nrA.begin();aI!=nrA.end();aI++)
		for(bI=nrB.begin();bI!=nrB.end();bI++)
			if((pX[*aI]!=0.0)&&(pY[*bI]!=0.0)&&(pXY[*aI][*bI]!=0.0))
				I+=pXY[*aI][*bI]*log2(pXY[*aI][*bI]/(pX[*aI]*pY[*bI]));
	return I;
	
}

double tGame::entropy(vector<int> list){
	map<int, double> p;
	map<int,double>::iterator pI;
	int i;
	double H=0.0;
	double c=1.0/(double)list.size();
	for(i=0;i<list.size();i++)
		p[list[i]]+=c;
	for (pI=p.begin();pI!=p.end();pI++) {
			H+=p[pI->first]*log2(p[pI->first]);	
	}
	return -1.0*H;
}

double tGame::ei(vector<int> A,vector<int> B,int theMask){
	set<int> nrA,nrB;
	set<int>::iterator aI,bI;
	map<int,map<int,double> > pXY;
	map<int,double> pX,pY;
	int i,j;
	double c=1.0/(double)A.size();
	double I=0.0;
	for(i=0;i<A.size();i++){
		nrA.insert(A[i]&theMask);
		nrB.insert(B[i]&theMask);
		pX[A[i]&theMask]=0.0;
		pY[B[i]&theMask]=0.0;
	}
	for(aI=nrA.begin();aI!=nrA.end();aI++)
		for(bI=nrB.begin();bI!=nrB.end();bI++){
			pXY[*aI][*bI]=0.0;
		}
	for(i=0;i<A.size();i++){
		pXY[A[i]&theMask][B[i]&theMask]+=c;
		pX[A[i]&theMask]+=c;
		pY[B[i]&theMask]+=c;
	}
	for(aI=nrA.begin();aI!=nrA.end();aI++)
		for(bI=nrB.begin();bI!=nrB.end();bI++)
			if((pX[*aI]!=0.0)&&(pY[*bI]!=0.0)&&(pXY[*aI][*bI]!=0.0))
				I+=pXY[*aI][*bI]*log2(pXY[*aI][*bI]/(pY[*bI]));
	return -I;
}
double tGame::computeAtomicPhi(vector<int>A,int states){
	int i;
	double P,EIsystem;
	vector<int> T0,T1;
	T0=A;
	T1=A;
	T0.erase(T0.begin()+T0.size()-1);
	T1.erase(T1.begin());
	EIsystem=ei(T0,T1,(1<<states)-1);
	P=0.0;
	for(i=0;i<states;i++){
		double EIP=ei(T0,T1,1<<i);
//		cout<<EIP<<endl;
		P+=EIP;
	}
//	cout<<-EIsystem+P<<" "<<EIsystem<<" "<<P<<" "<<T0.size()<<" "<<T1.size()<<endl;
	return -EIsystem+P;
}



double tGame::computeR(vector<vector<int> > table,int howFarBack){
	double Iwh,Iws,Ish,Hh,Hs,Hw,Hhws,delta,R;
	int i;
	for(i=0;i<howFarBack;i++){
		table[0].erase(table[0].begin());
		table[1].erase(table[1].begin());
		table[2].erase(table[2].begin()+(table[2].size()-1));
	}
	table[4].clear();
	for(i=0;i<table[0].size();i++){
		table[4].push_back((table[0][i]<<14)+(table[1][i]<<10)+table[2][i]);
	}
	Iwh=mutualInformation(table[0],table[2]);
    Iws=mutualInformation(table[0],table[1]);
    Ish=mutualInformation(table[1],table[2]);
    Hh=entropy(table[2]);
    Hs=entropy(table[1]);
    Hw=entropy(table[0]);
    Hhws=entropy(table[4]);
    delta=Hhws+Iwh+Iws+Ish-Hh-Hs-Hw;
    R=Iwh-delta;
  	return R;
}

double tGame::computeOldR(vector<vector<int> > table){
	double Ia,Ib;
	Ia=mutualInformation(table[0], table[2]);
	Ib=mutualInformation(table[1], table[2]);
	return Ib-Ia;
}

double tGame::predictiveI(vector<int>A){
	vector<int> S,I;
	S.clear(); I.clear();
	for(int i=0;i<A.size();i++){
		S.push_back((A[i]>>12)&15);
		I.push_back(A[i]&3);
	}
	return mutualInformation(S, I);
}

double tGame::nonPredictiveI(vector<int>A){
	vector<int> S,I;
	S.clear(); I.clear();
	for(int i=0;i<A.size();i++){
		S.push_back((A[i]>>12)&15);
		I.push_back(A[i]&3);
	}
	return entropy(I)-mutualInformation(S, I);
}
double tGame::predictNextInput(vector<int>A){
	vector<int> S,I;
	S.clear(); I.clear();
	for(int i=0;i<A.size();i++){
		S.push_back((A[i]>>12)&15);
		I.push_back(A[i]&3);
	}
	S.erase(S.begin());
	I.erase(I.begin()+I.size()-1);
	return mutualInformation(S, I);
}

void tGame::computeAllMI(char *filename){
    set<int> target,source;
    int i;
    for(i=0;i<9;i++)
        source.insert(i);
    makeAllSets(filename,target,source);
}

void tGame::makeAllSets(char *filename,set<int> target,set<int> source){
    if(source.size()==0){
        double d;
        set<int>::iterator SI;
        FILE *F=fopen(filename,"a+t");
        for(SI=target.begin();SI!=target.end();SI++){
            cout<<(*SI)+1;
            fprintf(F,"%i",(*SI+1));
        }
        d=doInformationCombination(target);
        fprintf(F," %i  %f\n",target.size(),d);
        cout<<" "<<target.size()<<" "<<d<<endl;
        fclose(F);
    } else {
        set<int>::iterator SI;
        int i;
        SI=source.begin();
        i=(*SI);
        source.erase(SI);
        makeAllSets(filename,target, source);
        target.insert(i);
        makeAllSets(filename,target, source);
    }
}

double tGame::doInformationCombination(set<int> combo){
    double r=0.0;
    vector<int> A,B;
    int i,j;
    int a,b;
    set<int>::iterator SI;
    for(i=0;i<83;i++)
        for(j=0;j<83;j++)
            if(i!=j){
                if(i<j)
                    b=0;
                else
                    b=1;
                a=0;
                for(SI=combo.begin();SI!=combo.end();SI++){


                    a=(a<<1)+scoreTable[i][(*SI)];
                    a=(a<<1)+scoreTable[j][(*SI)];
                }
                A.push_back(a);
                B.push_back(b);
            }
    return mutualInformation(A, B);
}




*/
