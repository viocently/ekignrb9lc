#ifndef __NEW_TRIVIUM_H__
#define  __NEW_TRIVIUM_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"

void newupdate(GRBModel& model, vector<GRBVar>& s,  bitset<288>& nonZeroTrackRep, bitset<288>& varsTrackRep);
void setRoundFlags(int r, bitset<80>& cube, bitset<288>& nonZeroFlag, bitset<288>& varsFlag, vector<bitset<288>>& nonZeroTrackReps, vector<bitset<288>>& varsTrackReps);
STATUS newMidSolutionCounter(int rounds, bitset<80>& cube, bitset<288>& last, map<bitset<288>,
    int, CMPS<288>> &counterMap, float time, int threads, int midround);

STATUS newMidSolutionCounter(int rounds, bitset<80>& cube, const bitset<288>& last,
    int& sols, float time, int threads);
#endif