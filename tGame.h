/*
 *  tGame.h
 *  HMMBrain
 *
 *  Created by Arend on 10/09/23.
 *  Adapted by Jory on 13/12/04
 *
 */
 
#ifndef _tGame_h_included_
#define _tGame_h_included_

#include "globalConst.h"
#include "tAgent.h"
#include "tTrial.h"
#include <vector>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


class tGame{
public:

    vector<tTrial*> trials;
    map<int, vector<int>> digitTrialIdx;
    int datasetSize;
	//map<pair<int,int>, tTrial*> trials;
	//void readImageData(string fn, map<pair<int,int>, tTrial*> &trials);
    matrix readPix(ifstream &);
    vector<pair<int,int>> xy;
    map< pair<int,int>, vector<vector<int>> > sequences;

    string executeGame(tAgent* agent, bool test, bool report, tTrial* trial) const;
    void evaluateAgent(tAgent* agent, tTrial* trial) const;
    tGame(string mnistFn, string saccadeFn);
    ~tGame();
    /*    
    double mutualInformation(vector<int> A,vector<int>B);
    double ei(vector<int> A,vector<int> B,int theMask);
    double computeAtomicPhi(vector<int>A,int states);
    double predictiveI(vector<int>A);
    double nonPredictiveI(vector<int>A);
    double predictNextInput(vector<int>A);
    double computeR(vector<vector<int> > table,int howFarBack);
    double computeOldR(vector<vector<int> > table);
    double entropy(vector<int> list);
		
    void computeAllMI(char *filename);
    void makeAllSets(char *filename,set<int> target,set<int> source);
    double doInformationCombination(set<int> combo);
    */

};
#endif
