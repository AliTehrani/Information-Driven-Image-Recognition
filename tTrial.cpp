#include "tTrial.h"
#include <sstream>
#include <iostream>
using namespace std;

tTrial::tTrial(){
}

tTrial::~tTrial(){
}

tTrial::tTrial(int dig, int rep, vector<vector<int>> pix): digit(dig), rep(rep), image(pix){}
