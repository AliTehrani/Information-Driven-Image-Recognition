/*
 *  tHMM.cpp
 *  HMMBrain
 *
 *  Created by Arend on 9/16/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>
using std::pow;
#include <map>
using std::map;

#include "tHMM.h"
//#define feedbackON

tHMMU::tHMMU(){
}

tHMMU::~tHMMU(){
	hmm.clear();
	sums.clear();
	ins.clear();
	outs.clear();
	posLevelOfFB.clear();
	negLevelOfFB.clear();
	chosenInPos.clear();
	chosenInNeg.clear();
	chosenOutPos.clear();
	chosenOutNeg.clear();
}

void tHMMU::setup(vector<unsigned char> &genome, int start){
	int i,j,k;
	ins.clear();
	outs.clear();
	k=(start+2)%genome.size();

	_xDim=1+(genome[(k++)%genome.size()]&3);
	_yDim=1+(genome[(k++)%genome.size()]&3);
	posFBNode=genome[(k++)%genome.size()]&(brainSize-1);
	negFBNode=genome[(k++)%genome.size()]&(brainSize-1);
	nrPos=genome[(k++)%genome.size()]&3;
	nrNeg=genome[(k++)%genome.size()]&3;
	cout<<"setup "<<(int)genome[start+2]<<" "<<(int)_xDim<<" "<<(int)_yDim<<endl;
	ins.resize(_yDim);
	outs.resize(_xDim);
	posLevelOfFB.resize(nrPos);
	negLevelOfFB.resize(nrNeg);
	for(i=0;i<_yDim;i++)
		ins[i]=genome[(k+i)%genome.size()]&(brainSize-1);
	for(i=0;i<_xDim;i++)
		outs[i]=genome[(k+4+i)%genome.size()]&(brainSize-1);
	for(i=0;i<nrPos;i++)
		posLevelOfFB[i]=(int)(1+genome[(k+8+i)%genome.size()]);
	for(i=0;i<nrNeg;i++)
		negLevelOfFB[i]=(int)(1+genome[(k+12+i)%genome.size()]);
	chosenInPos.clear();
	chosenInNeg.clear();
	chosenOutPos.clear();
	chosenOutNeg.clear();
	
	k=k+16;
	hmm.resize(1<<_yDim);
	sums.resize(1<<_yDim);
	for(i=0;i<(1<<_yDim);i++){
		hmm[i].resize(1<<_xDim);
		for(j=0;j<(1<<_xDim);j++){
//			hmm[i][j]=(genome[(k+j+((1<<yDim)*i))%genome.size()]&1)*255;
			hmm[i][j]=genome[(k+j+((1<<_xDim)*i))%genome.size()];
			if(hmm[i][j]==0) hmm[i][j]=1;
			sums[i]+=hmm[i][j];
		}
	}
}

void tHMMU::setupQuick(vector<unsigned char> &genome, int start){
	int i,j,k;
	ins.clear();
	outs.clear();
	k=(start+2)%genome.size();
	
	_xDim=1+(genome[(k++)%genome.size()]&3);
	_yDim=1+(genome[(k++)%genome.size()]&3);
	posFBNode=genome[(k++)%genome.size()]&(brainSize-1);
	negFBNode=genome[(k++)%genome.size()]&(brainSize-1);
	nrPos=genome[(k++)%genome.size()]&3;
	nrNeg=genome[(k++)%genome.size()]&3;
	//cout<<"setup "<<(int)genome[start+2]<<" "<<(int)xDim<<" "<<(int)yDim<<endl;
	ins.resize(_yDim);
	outs.resize(_xDim);
	posLevelOfFB.resize(nrPos);
	negLevelOfFB.resize(nrNeg);
	for(i=0;i<_yDim;i++)
		ins[i]=genome[(k+i)%genome.size()]&(brainSize-1);
	for(i=0;i<_xDim;i++)
		outs[i]=genome[(k+4+i)%genome.size()]&(brainSize-1);
	for(i=0;i<nrPos;i++)
		posLevelOfFB[i]=(int)(1+genome[(k+8+i)%genome.size()]);
	for(i=0;i<nrNeg;i++)
		negLevelOfFB[i]=(int)(1+genome[(k+12+i)%genome.size()]);
	chosenInPos.clear();
	chosenInNeg.clear();
	chosenOutPos.clear();
	chosenOutNeg.clear();
	
	k=k+16;
	hmm.resize(1<<_yDim);
	sums.resize(1<<_yDim);
	for(i=0;i<(1<<_yDim);i++){
		hmm[i].resize(1<<_xDim);
		for(j=0;j<(1<<_xDim);j++)
			hmm[i][j]=0;
		hmm[i][genome[(k+j+((1<<_xDim)*i))%genome.size()]&((1<<_xDim)-1)]=255;
		sums[i]=255;
		//cout << (int)genome[(k+j+((1<<_xDim)*i))%genome.size()] << endl;
	}
	
}

void tHMMU::update(unsigned char *states,unsigned char *newStates,int *nrIncomingSignals){
	int I=0;
	int i,j,r;
	unsigned char mod;
#ifdef feedbackON
	if((nrPos!=0)&&(states[posFBNode]==1)){
		for(i=0;i<chosenInPos.size();i++){
			mod=(unsigned char)(rand()%(int)posLevelOfFB[i]);
			if((hmm[chosenInPos[i]][chosenOutPos[i]]+mod)<255){
				hmm[chosenInPos[i]][chosenOutPos[i]]+=mod;
				sums[chosenInPos[i]]+=mod;
			}
		}
	}
	if((nrNeg!=0)&&(states[negFBNode]==1)){
		for(i=0;i<chosenInNeg.size();i++){
			mod=(unsigned char)(rand()%(int)negLevelOfFB[i]);
			if((hmm[chosenInNeg[i]][chosenOutNeg[i]]-mod)>0){
				hmm[chosenInNeg[i]][chosenOutNeg[i]]-=mod;
				sums[chosenInNeg[i]]-=mod;
			}
		}
	}
#endif
	for(i=0;i<ins.size();i++)
		I=(I<<1)+((states[ins[i]])&1);
	r=1+(rand()%(sums[I]-1));
	j=0;
//	cout<<I<<" "<<(int)hmm.size()<<" "<<(int)hmm[0].size()<<endl;
	while(r>hmm[I][j]){
		r-=hmm[I][j];
		j++;
	}
	for(i=0;i<outs.size();i++){
	    nrIncomingSignals[outs[i]]++;
	    newStates[outs[i]]+=(j>>i)&1;
	    //newStates[outs[i]]|=(j>>i)&1; //commented out to implement synaptic inhibition
	    //newStates[outs[i]]=(j>>i)&1;
	}
#ifdef feedbackON
	chosenInPos.push_back(I);
	chosenInNeg.push_back(I);
	chosenOutPos.push_back(j);
	chosenOutNeg.push_back(j);
	while(chosenInPos.size()>nrPos) chosenInPos.pop_front();
	while(chosenOutPos.size()>nrPos) chosenOutPos.pop_front();
	while(chosenInNeg.size()>nrNeg) chosenInNeg.pop_front();
	while(chosenOutNeg.size()>nrNeg) chosenOutNeg.pop_front();
#endif
}

void tHMMU::show(void){
	int i,j;
	cout<<"INS: ";
	for(i=0;i<ins.size();i++)
		cout<<(int)ins[i]<<" ";
	cout<<endl;
	cout<<"OUTS: ";
	for(i=0;i<outs.size();i++)
		cout<<(int)outs[i]<<" ";
	cout<<endl;
	for(i=0;i<hmm.size();i++){
		for(j=0;j<hmm[i].size();j++)
			cout<<" "<<(double)hmm[i][j]/sums[i];
		cout<<endl;
	}
	cout<<endl;
	/*
	cout<<"posFB: "<<(int)posFBNode<<" negFB: "<<(int)negFBNode<<endl;
	cout<<"posQue:"<<endl;
	for(i=0;i<posLevelOfFB.size();i++)
		cout<<(int)posLevelOfFB[i]<<" ";
	cout<<endl;
	cout<<"negQue:"<<endl;
	for(i=0;i<negLevelOfFB.size();i++)
		cout<<(int)negLevelOfFB[i]<<" ";
	cout<<endl;
	*/
	// @AT
	if (ins.size()==2&&outs.size()==1)
	  cout<<'\n'<<determineGateType().first<<"\n\n";

/*
	for(i=0;i<hmm.size();i++){
		for(j=0;j<hmm[i].size();j++)
			cout<<(int)hmm[i][j]<<" ";
		cout<<endl;
	}
	*/
//	cout<<"------"<<endl;
}

// @AT ->
void tHMMU::save(ofstream &phenotypeFile){
	int i,j;
	phenotypeFile<<"INS: ";
	for(i=0;i<ins.size();i++)
		phenotypeFile<<(int)ins[i]<<" ";
	phenotypeFile<<endl;
	phenotypeFile<<"OUTS: ";
	for(i=0;i<outs.size();i++)
		phenotypeFile<<(int)outs[i]<<" ";
	
	phenotypeFile<<endl;
	for(i=0;i<hmm.size();i++){
		for(j=0;j<hmm[i].size();j++)
			phenotypeFile<<" "<<(double)hmm[i][j]/sums[i];
		phenotypeFile<<endl;
	}
	phenotypeFile<<endl;
	phenotypeFile<<"posFB: "<<(int)posFBNode<<" negFB: "<<(int)negFBNode<<endl;
	phenotypeFile<<"posQue:"<<endl;
	for(i=0;i<posLevelOfFB.size();i++)
		phenotypeFile<<(int)posLevelOfFB[i]<<" ";
	phenotypeFile<<endl;
	phenotypeFile<<"negQue:"<<endl;
	for(i=0;i<negLevelOfFB.size();i++)
		phenotypeFile<<(int)negLevelOfFB[i]<<" ";
	phenotypeFile<<endl;
	
	if (ins.size()==2&&outs.size()==1)
	  phenotypeFile<<'\n'<<determineGateType().first<<"\n\n";

}

/*
0 0000 ZERO
1 0001 NOR
2 0010 AND-NOT
3 0011 NOT
4 0100 AND-NOT
5 0101 NOT
6 0110 XOR
7 0111 NAND
8 1000 AND
9 1001 EQU
10 1010 COPY
11 1011 OR-NOT
12 1100 COPY
13 1101 OR-NOT
14 1110 OR
15 1111 ONE
*/
// 0: zero, 1: NOR, 2,4: AND-NOT, 3: COPY, 
int tHMMU::gateIndex(){
  unsigned int o[4]={hmm[0][1]/sums[0], hmm[1][1]/sums[1], hmm[2][1]/sums[2], hmm[3][1]/sums[3]};
  int indx=0;
  for (int i=0; i<4; i++)
    indx += o[i]*pow(2,i);
  return indx;
}

pair<string,int> tHMMU::determineGateType(){
  string gateType="";
  int input=-1;
  if (ins.size()!=2||outs.size()!=1)
    return make_pair(gateType, input);

  map<int, string> gateMap { {0, "ZERO"}, {1, "NOR"}, {2, "AND-NOT"}, {3,"NOT"}, {4, "AND-NOT"}, {5, "NOT"}, {6, "XOR"}, {7, "NAND"}, {8, "AND"}, {9, "EQU"}, {10, "COPY"}, {11, "OR-NOT"}, {12, "COPY"}, {13, "OR-NOT"}, {14, "OR"}, {15, "ONE"} };
  map<int, int> inputMap { {2, 0}, {3, 0}, {4, 1}, {5, 1}, {10, 1}, {11, 0}, {12, 0}, {13, 1} };
  
  if (gateMap.count(gateIndex()) != 0)
      gateType=gateMap.at(gateIndex());
  if (inputMap.count(gateIndex()) != 0)
      input = inputMap.at(gateIndex());
  
  return make_pair(gateType, input);
  //determine the type of 2-to-1 gates
  unsigned int o[4]={hmm[0][1]/sums[0], hmm[1][1]/sums[1], hmm[2][1]/sums[2], hmm[3][1]/sums[3]};
  //printf("truth table: %d %d %d %d\n", o[0], o[1], o[2], o[3]);
  int nrOnes=o[0]+o[1]+o[2]+o[3];
  switch(nrOnes){
  case 0 :
    gateType="ZERO";
    break;
  case 1:
    if (o[0]==1)
      gateType="NOR";
    else if (o[1]==1 || o[2]==1){
      gateType="AND-NOT";
      input=o[2]; //o[2]=1 means: (1 AND-NOT 0)=1, second input (input 1) is negated
    }
    else
      gateType="AND";
    break;
  case 2:
    if ((o[0]==1&&o[1]==1)||(o[0]==1&&o[2]==1)){
      gateType="NOT";
      input=o[2];//o[0] is 0 for both. o[2]=1 means: ~(10)=1, second input (input 1) is negated
    }
    else if ((o[2]==1&&o[3]==1)||(o[1]==1&&o[3]==1)){
      gateType="COPY";
      input=o[1];//o[1]=1 means: copy(01)=1, second input is copied
    }
    else if (o[0]==1&&o[3]==1)
      gateType="XNOR";
    else 
      gateType="XOR";
    break;
  case 3:
    if (o[0]==0)
      gateType="OR";
    else if (o[3]==0)
      gateType="NAND";
    else{
      gateType="OR-NOT";//TODO
      input=o[2];//o[2]=1 means: (1 OR-NOT 0)=1, second input is negated 
    }
    break;
  case 4:
    gateType="ONE";
    break;
  }
  return make_pair(gateType, input);
}
// <- @AT 
