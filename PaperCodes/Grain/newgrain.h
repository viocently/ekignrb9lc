#ifndef __NEW_GRAIN_H__
#define  __NEW_GRAIN_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

void newupdate(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& b, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep);

STATUS newMidSolutionCounter(int rounds, bitset<96>& cube,
    bitset<256>& last, map<bitset<256>, int, CMPS<256>>& counterMap, float
    time, int thread, int midround);
STATUS newMidSolutionCounter(int rounds, bitset<96>& cube,
    bitset<256>& last, int& sols, float
    time, int thread);



#endif