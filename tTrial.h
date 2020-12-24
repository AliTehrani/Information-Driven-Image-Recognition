#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <utility>
#include <map>

using namespace std;

using matrix=vector<vector<int>>;

class tTrial{
public:
  //int index;
  //string category, label;
    
    int digit, rep;
    vector<vector<int>> image;

    double rew=1, penalty=1.0/9.0;

    tTrial(int dig, int rep, vector<vector<int>> pix);
    tTrial();
    ~tTrial();

    //void printTrial(ostream &out);
    //int MD();
    //void doIOI(int pos);
    //void doOB(int dec);
};
