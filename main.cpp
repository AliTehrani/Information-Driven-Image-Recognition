#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <map>
#include <math.h>
#include <time.h>
#include <iostream>
#include "globalConst.h"
#include "tHMM.h"
#include "tAgent.h"
#include "tGame.h"
#include <string.h>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <iomanip>

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

using namespace std;

int batchSize=100;
double maxScore=double(batchSize);//, rewPd, rewSt, rewNd;
string datasetFn="mnist.100.txt";
string saccadeFn="saccades.txt";
string runId="", outputsFn;

//double replacementRate=0.1;
double perSiteMutationRate=0.005;
//int update=0; // generation counter
int maxAgents=100; // population size
int totalGenerations=1000;
//double penaltyRatio=0.8;
int nrArg=1;
//int argi=11;
int rElitism=5;
int repeats;


// find maximum score of a subset of trials. to normalize fitness values with
double findMaxScore(vector<tTrial*> trials);

// calculate total fitness by suming up fitnesses vector and normalize by 100
//void sumFitnesses(vector<tAgent*> &agents);//, double maxScore);

// make a vector of batches, a batch is a vector of trials
vector<vector<tTrial*>> makeBatches(tGame* game, int batchSize);

// read output file
void readOutputFile(int argc, char* argv[]);

// run game for the population
void runGame(tGame* game, vector<tAgent*> &agents, vector<tTrial*> trials, bool test, bool report);

// calculate total fitness of each agent and find the best fitness
int findBestAgent(vector<tAgent*> &agents);

// apply elitism algorithm in selection with elitism rate rElitism
void elitism(vector<tAgent*> &agents, vector<tAgent*> &nextGen, tAgentSetter* agentSetter, int update);

// selection
void selection(vector<tAgent*> &agents, vector<tAgent*> &nextGen, tAgentSetter* agentSetter, double maxFitness, int update);

// write output files
void writeOutputFiles(tAgent* &bestAgent);

// write genome file of the best at generation gen
void writeGenomeFile(tAgent* &bestAgent, int gen);

// find best agent at the end of run and write outputfiles
//void findBestRun(tGame *game, vector<tAgent*> &agents, vector<tTrial*> trials);

// execute game (train and test) for one agent only and write output files
void testAgent(tGame *game, tAgent* &agent, vector<tTrial*> trials);

// perform knockout tests
void knockOutTest(tGame* game, tAgent* agent, vector<tTrial*> trials, string koFn);

// custom knockout
void customKO(tAgent* &agent, vector<int> gates);


int main(int argc, char *argv[]) { // ./acc -dataset (dataset filename) -batchsize (batchSize) -r (random seed) -bu (brain updates) -s (master agent) -test -gen (totalGenerations) -wg (write output gens) -o (outputListFn) -id (runID) 

        long long t0=clock();
	vector<tAgent*> agents;
	vector<tAgent*> nextGen; // write-buffer for population
	tAgent *masterAgent; // used to clone and start LoD
	int i,j;

	// read experiment parameters
	//char* trPatFn=argv[nrArg++];
	/*
	rewPd = atof(argv[nrArg++]);
	rewSt = atof(argv[nrArg++]);
	rewNd = atof(argv[nrArg++]);
	*/

	// setup game trials
	// input filename can be read from the command line
	// write a function to create patterns, trials, etc.
	if (nrArg<argc and strcmp(argv[nrArg], "-dataset")==0){
	    datasetFn=argv[++nrArg];
	    ++nrArg;
	}

	// saccade sequence files
	if (nrArg<argc and strcmp(argv[nrArg], "-xy")==0){
	    saccadeFn=argv[++nrArg];
	    ++nrArg;
	}
	tGame* game=new tGame(datasetFn, saccadeFn);
	/*
	cout << "dataset size:\t" << game->datasetSize << '\t' << batchSize << endl;	
	for (auto p: game->xy){
	    cout << p.first << "\t" << p.second << endl;
	}
	int q; cin >>q;
	*/
	// batch size
	if (nrArg<argc and strcmp(argv[nrArg], "-batchsize")==0){
	    batchSize=atoi(argv[++nrArg]);
	    ++nrArg;
	}

	// set random seed
	int randSeed;
	if (nrArg<argc and strcmp(argv[nrArg], "-r")==0){
	    randSeed=atoi(argv[++nrArg]);
	    ++nrArg;
	}
	else 
	    randSeed=getpid();

	srand(randSeed); // need different includes for windows XPLATFORM
	cout << "Random seed: " << randSeed << endl;

	// agent setter: setup that all agents share
	const int nrIns=4, nrOuts=10;
	int ins[nrIns]={0,1,2,3}, outs[nrOuts]={54,55,56,57,58,59,60,61,62,63};

	int brainUpdates=1;
	if (nrArg<argc and strcmp(argv[nrArg], "-bu")==0){
	    brainUpdates=atoi(argv[++nrArg]);
	    ++nrArg;
	}

	tAgentSetter *agentSetter=new tAgentSetter(brainUpdates, nrIns, nrOuts, ins, outs);

	/// setup initial population
	masterAgent=new tAgent(agentSetter); // can also loadAgent(specs...)
	
	// make initial population seeding master agent
	if (nrArg<argc and strcmp(argv[nrArg], "-s")==0){
	    char* savedAgent=argv[++nrArg];
	    masterAgent->loadAgent(savedAgent);
	    ++nrArg;
	    cout << "Master agent: " << savedAgent << endl;

	    // knockout tests
	    if (nrArg<argc and strcmp(argv[nrArg], "-ko")==0){
		string koFn=argv[++nrArg];
		knockOutTest(game, masterAgent, game->trials, koFn);
		return 0;
	    }
	    // test agent only
	    if (nrArg<argc and strcmp(argv[nrArg], "-test")==0){
		++nrArg;
		testAgent(game, masterAgent, game->trials);
		cout << "fitness: " << masterAgent->fitness << "\taccuracy: " << masterAgent->accuracy << endl;
		readOutputFile(argc, argv);
		writeOutputFiles(masterAgent);
		long long t1=clock();
		cout << "time to run this experiment: " << (t1-t0)*1e-6 <<endl;
		return 0;
	    }
	}
	else{
	    masterAgent->setupRandomAgent(10000, 20); // create with 5000 genes, 30 start codons
	    masterAgent->determinePhenotype();
	}

	/*
	cout << "before custom ko: " << masterAgent->hmmus.size() << endl;
	customKO(masterAgent, vector<int> {5});
	cout << "after custom ko: " << masterAgent->hmmus.size() << endl;
	*/

	// number of generations in evolution
	if (nrArg<argc and strcmp(argv[nrArg], "-gen")==0){
	    totalGenerations = atoi(argv[++nrArg]);
	    ++nrArg;
	}

	agents.resize(maxAgents);
	
	for(int i=0; i<agents.size(); i++){
	    agents[i]=new tAgent(agentSetter);
	    agents[i]->inherit(masterAgent, 0.01, 0); // {}(*from, mu, t);
	}
	nextGen.resize(agents.size());
	masterAgent->nrPointingAtMe--; // effectively kill first ancestor from whom all inherit

	// write output files (if genome file is saved it will be sufficient to generate further outputs) in these generations 
	vector<int> writeOutputGens;
	if (nrArg<argc and strcmp(argv[nrArg],"-wg")==0){
	    char* outGenFn=argv[++nrArg];
	    nrArg++;	

	    // read writing output generations
	    ifstream wog(outGenFn);
	    long g;
	    while(wog>>g)
	        writeOutputGens.push_back(g);
	}

	// read output filename and run id
	readOutputFile(argc, argv);

	//cout << "dataset:\t" << datasetFn << endl << "batch size:\t" << batchSize << endl;

	auto batchTrials = makeBatches(game, batchSize);
	cout << "Batches made\t" << batchTrials.size() << endl;

	cout.precision(2);
	cout << fixed;
	///  the main loop
	int update, epoch=0, datasetSize=game->trials.size();
	for (update=0; update<totalGenerations; update++){
	    // make a random subset of trials. usually from pre-made subsets, where each pre-made subset contains only of a particular outcome
	    // here I use all possible trials for evolution
	    if (batchTrials.size()==0){
		epoch++;
		batchTrials = makeBatches(game, batchSize);
	    }
	    auto trials = batchTrials.back();
	    batchTrials.pop_back();
	    /*
	    for (auto trial: trials)
		cout << trial->rep << " ";
	    cout << endl;
	    int q; cin >>q;
	    */
	    //cout << "Trials size\t"<< trials.size() << endl;
	    
	    // find maximum score of trials. This is necessary if agents are evaluated in a subset of trials that changes over generations
	    //double maxScore = findMaxScore(game->trials);

	    // execute game for population
	    bool is_test=false, is_report=false;
	    runGame(game, agents, trials, is_test, is_report);
		
	    // calculate total fitness and normalize by 100
	    //sumFitnesses(agents);//, maxScore);

	    // find the best fitness
	    int bestAgent = findBestAgent(agents);
	    double maxFitness = agents[bestAgent]->fitness, maxAcc = agents[bestAgent]->accuracy;

	    // print fitness and learned size
	    cout << "update: "<< update << " epoch " << epoch << " " << maxFitness << " " << maxAcc << endl;//" " ;//<< endl;
	    //copy(agent[bestAgent]->notLearned.begin(), agent[bestAgent]->notLearned.end(), ostream_iterator<string>(cout, ","));
	    /*
	      for (set<string>::iterator it=agents[bestAgent]->notLearned.begin(); it!=agents[bestAgent]->notLearned.end(); it++)
	      cout << (*it) << ",";
	      cout << endl;
	    */

	    // find elites and move them to next generation
	    elitism(agents, nextGen, agentSetter, update);

	    // selection
	    selection(agents, nextGen, agentSetter, maxFitness, update);

	    if (find(writeOutputGens.begin(), writeOutputGens.end(), update)!=writeOutputGens.end()){
		writeGenomeFile(agents[bestAgent], update);
	    }
	    //cout << update << " TODO" << endl;
	    //  findBestRun(game, agent, argc, argv, update);
	}
	
	// find best run
	cout << "find best agent\n";
	bool is_test=false, is_report=false;
	runGame(game, agents, game->trials, is_test, is_report);
	int bestAgent = findBestAgent(agents);
	double bestFitness = agents[bestAgent]->fitness, bestAcc = agents[bestAgent]->accuracy;
	cout << "fitness: "<< bestFitness << "\taccuracy: " << bestAcc << endl;
	writeOutputFiles(agents[bestAgent]);

	delete agentSetter;
	delete game;
	long long t1=clock();
	cout << "time to run this experiment: " << (t1-t0)*1e-6 <<endl;
	return 0;
}

// make a vector of batches, a batch is a vector of trials                                                                                                                                                      
vector<vector<tTrial*>> makeBatches(tGame* game, int batchSize){
    int nrBatches = game->datasetSize/batchSize + ((game->datasetSize%batchSize==0) ? 0 : 1);
    vector<vector<tTrial*>> batches(nrBatches);
    for (auto item: game->digitTrialIdx){
        // indices of trials for a particular digit                                                                                                                                                             
        auto indices = item.second;
        random_shuffle(indices.begin(), indices.end());
        int batchPerClass = indices.size()/nrBatches;// if number of images in each class is not uniform this would leave out some of the images                                                                
        for (int b=0; b<nrBatches; b++){
            for (auto i=b*batchPerClass; i<indices.size() and i<(b+1)*batchPerClass; i++){
                batches[b].push_back( game->trials[indices[i]] );
            }
        }
    }
    return batches;
}

// read output file
void readOutputFile(int argc, char* argv[]){
    if (nrArg<argc and strcmp(argv[nrArg], "-o")==0){
	outputsFn = argv[++nrArg];
    }
    nrArg++;
    if (nrArg<argc and strcmp(argv[nrArg], "-id")==0){
	runId = argv[++nrArg];
	runId.insert(0, "_");
    }
}

// find maximum score of a subset of trials. to normalize fitness values with
double findMaxScore(vector<tTrial*> trials){
    double maxScore=0;
    for (auto trial: trials)
        maxScore += trial->rew;

    return maxScore;
}
/*
// sum up fitnesses in fitnesses vector and normalize by 100
void sumFitnesses(vector<tAgent*> &agents){ //, double maxScore){
    for (int i=0; i<agents.size(); i++){
        double fitness=0.0;
	for (int j=0; j<agents[i]->fitnesses.size(); j++){
	    fitness += agents[i]->fitnesses[j];
	}
	agents[i]->fitness = 100.0 * fitness/maxScore;
    }
}
*/
// run game for the population
void runGame(tGame* game, vector<tAgent*> &agents, vector<tTrial*> trials, bool test, bool report){
    map<int,int> nrPerDigit;
    for (auto trial: trials){
	nrPerDigit[trial->digit]++;
    }
    for(int i=0; i<agents.size(); i++){
        // reset fitness, fitnesses
	//agents[i]->fitnesses.clear();
	agents[i]->correct.clear();
	agents[i]->accs.clear();
	agents[i]->fitness=0.0;
	agents[i]->accuracy=0.0;
	// execute game
	for(auto trial: trials){
	    // reset fitness, fitnesses
	    //agents[i]->fitness=0.0;
	    game->executeGame(agents[i], test, report, trial);
	    //agents[i]->fitnesses.push_back(agents[i]->fitness);
	}
	agents[i]->fitness *= 100.0/trials.size();//maxScore;
	int correctSum=0;
	for (auto dn: nrPerDigit){
	    int d=dn.first, n=dn.second;
	    correctSum += agents[i]->correct[d];
	    agents[i]->accs[d] = agents[i]->correct[d]/double(n);
	}
	agents[i]->accuracy = correctSum/double(trials.size());
	//if (test)
	//cout << agents[i]->fitness << '\t';
    }
    //if (test)
    //cout << endl;
}

// calculate total fitness of each agent and find the agent with best fitness
int findBestAgent(vector<tAgent*> &agents){
    double maxFitness=0.0;
    int bestAgent=0;
    for(int i=0; i<agents.size(); i++){
	//find best fitness
	if(agents[i]->fitness>maxFitness){
	    maxFitness=agents[i]->fitness;
	    bestAgent=i;
	}
    }
    return bestAgent;
}

// apply elitist selection with elitism rate rElitism
void elitism(vector<tAgent*> &agents, vector<tAgent*> &nextGen, tAgentSetter* agentSetter, int update){
    vector<int> elites;
    for (int r=0; r<rElitism; r++){
        double eliteFitness=0.0;
	int elite=0;
	for (int i=0; i<agents.size(); i++){
	    // if fitness>elites[-1].fitness then continue 
	    if (find(elites.begin(), elites.end(), i) != elites.end()) // if current i exists in elites then skip it
	      continue;

	    if(agents[i]->fitness>eliteFitness){
	        eliteFitness=agents[i]->fitness;
		elite=i;
	    }
	}
	elites.push_back(elite);

	//move elite r to next generation with no mutation
	tAgent *d=new tAgent(agentSetter);
	d->inherit(agents[elite], 0, update);
	nextGen[r]=d;
    }
}

// selection
void selection(vector<tAgent*> &agents, vector<tAgent*> &nextGen, tAgentSetter* agentSetter, double maxFitness, int update){
    //roulette wheel selection
    for(int i=rElitism; i<agents.size(); i++){
        tAgent *d;
	d=new tAgent(agentSetter);
	int j;
	do{
	    j=rand()%(int)agents.size();
	} while(randDouble>(pow(1.05, agents[j]->fitness) / pow(1.05, maxFitness)));
        
	d->inherit(agents[j], perSiteMutationRate, update);
	nextGen[i]=d;
    }

    // moves "nextGen" to current population
    for(int i=0; i<agents.size(); i++){
        agents[i]->retire();
	agents[i]->nrPointingAtMe--;
	/***************************** THIS IS JUST A TEST TO DETECT MEMORY LEAK [by ALI] *****************************/
	// KILL ANCESTOR NO MATTER WHAT
	if (agents[i]->nrPointingAtMe==0)
	    delete agents[i];

	//agents[i]=nextGen[i];
    }
    agents.clear();
    agents.resize(maxAgents);
    for (int i=0; i<agents.size(); i++){	
	agents[i]=nextGen[i];
	//agents[i]->ancestor = NULL;
    }
    
    agents=nextGen;
}


// write genome file at writeOutputGen
void writeGenomeFile(tAgent* &bestAgent, int gen){

    map<string, string> outList;
    ifstream infile(outputsFn);
    string sw, fn;
    while (infile >> sw >> fn){
        //outList.push_back(make_pair(sw, fn));                                                                                                                                                                 
        outList[sw] = fn;
    }
    if (outList.find("-g")!=outList.end()){
        ostringstream gFn;
        gFn << outList["-g"] << runId << "_" << gen;
        //FILE * genomeFile=fopen(gFn.str().c_str(), "w+t");
	cout << "write best agent genome ..." << endl;
        //bestAgent->saveGenome(genomeFile);
        bestAgent->saveGenome(gFn.str());
    }

}
// write output files
void writeOutputFiles(tAgent* &bestAgent){
    cout << "write outputs" << endl;

    map<string, string> outList;
    ifstream infile(outputsFn);
    string sw, fn;
    while (infile >> sw >> fn){
        //outList.push_back(make_pair(sw, fn));
        outList[sw] = fn;
    }
    if (outList.find("-l")!=outList.end()){
	ostringstream lodFn;
	lodFn << outList["-l"] << runId;
	//FILE *LOD=fopen(lodFn.str().c_str(), "w+t");
	ofstream LOD(lodFn.str().c_str());
	bestAgent->saveLOD(LOD);//,genomeFile);
    }
    if (outList.find("-g")!=outList.end()){
        ostringstream gFn;
	gFn << outList["-g"] << runId;
	//FILE * genomeFile=fopen(gFn.str().c_str(), "w+t");
        //bestAgent->saveGenome(genomeFile);
	bestAgent->saveGenome(gFn.str());
    }
    if (outList.find("-lt")!=outList.end()){
        ostringstream ltFn;
	ltFn << outList["-lt"] << runId;
	FILE* ltFile = fopen(ltFn.str().c_str(), "w+t");
	bestAgent->saveLogicTable(ltFile);
    }
    if (outList.find("-ph")!=outList.end()){
        ostringstream phFn;
	phFn << outList["-ph"] << runId;
	//FILE* ltFile = fopen(ltFn.str().c_str(), "w+t");
	bestAgent->savePhenotype(phFn.str().c_str());
    }
    if (outList.find("-lc")!=outList.end()){
        ostringstream lcFn;
	lcFn << outList["-lc"] << runId;
	//FILE* ltFile = fopen(ltFn.str().c_str(), "w+t");
	bestAgent->saveToDotFileDiagram(lcFn.str());
    }
        
    /*
    vector<int> transitions=agent[bestAgent]->attention();
    int nonZeroTr=0;
    for (vector<int>::iterator tr=transitions.begin(); tr!=transitions.end(); tr++)
    if (*tr>1)
    nonZeroTr++;
    cout << transitions.size() << '\t' << nonZeroTr << endl;
    */
}

// perform knockout tests
void knockOutTest(tGame* game, tAgent* agent, vector<tTrial*> trials, string koFn){
    cout << "ko\n";

    // calculate fitness of the agent in the absence of knockouts
    vector<tAgent*> agents{agent};
    runGame(game, agents, trials, false, false);
    maxScore = findMaxScore(trials);
    agent->fitness *= 100.0/maxScore;
    double fitnessOrig=agent->fitness;

    ofstream koFile(koFn);
    koFile << ",gateIndex,fitnessOrig,fitnessLoss,secondFitnessLoss\n";
    for (int i=0; i<agent->hmmus.size(); i++){
        tAgent* koAgent = new tAgent(*agent);

	int gi= koAgent->hmmus[i]->gateIndex();
	koAgent->knockOutGate(i);

	// execute game for the ko agent
	bool is_test=true, is_report=true;
	vector<tAgent*> koAgents{koAgent};
	runGame(game, koAgents, trials, is_test, is_report);

	// normalize fitness
	maxScore = findMaxScore(trials);
	koAgent->fitness *= 100.0/maxScore;
	double firstKoFitness = koAgent->fitness;
	
	// second knockouts
	vector<double> secondKoFitnesses;
	for (int j=0; j<agent->hmmus.size(); j++){
	    if (i==j){
	      secondKoFitnesses.push_back(fitnessOrig - firstKoFitness);
	      continue;
	    }
	        
	    tAgent* secondKoAgent = new tAgent(*agent);
	    
	    int gj= secondKoAgent->hmmus[j]->gateIndex();

	    secondKoAgent->knockOutGate(i);
	    int k = j - int(j>i);
	    secondKoAgent->knockOutGate(k);

	    // execute game for the ko agent
	    bool is_test=true, is_report=true;
	    vector<tAgent*> secondKoAgents{secondKoAgent};
	    runGame(game, secondKoAgents, trials, is_test, is_report);

	    // normalize fitness
	    maxScore = findMaxScore(trials);
	    secondKoAgent->fitness *= 100.0/maxScore;
	    
	    secondKoFitnesses.push_back(fitnessOrig - secondKoAgent->fitness);
	}

	koFile << i << "," << gi <<  "," << fitnessOrig << "," << fitnessOrig - firstKoFitness << ",";
	ostringstream fits;
	copy(secondKoFitnesses.begin(), secondKoFitnesses.end(), ostream_iterator<double>(fits, ","));
	koFile  << "\"" << fits.str().substr(0, fits.str().size()-1) << "\"" << endl;
    }    
}

void customKO(tAgent* &agent, vector<int> gates){
    for (auto g: gates)
        agent->knockOutGate(g);
}

// test agent only
void testAgent(tGame *game, tAgent* &agent, vector<tTrial*> trials){
    cout << "test agent\n";
    
    vector<tAgent*> testAgent{agent};

    // execute game for the test Agent
    bool is_test=true, is_report=true;
    runGame(game, testAgent, trials, is_test, is_report);
    
    //maxScore = findMaxScore(trials);
    //agent->fitness *= 100.0/maxScore;

    //writeOutputFiles(argc, argv, agent);

    /*
    vector<int> transitions=agent->attention();
    int nonZeroTr=0;
    for (vector<int>::iterator tr=transitions.begin(); tr!=transitions.end(); tr++)
      if (*tr>1)
        nonZeroTr++;
    cout << transitions.size() << '\t' << nonZeroTr << endl;
    */
}
