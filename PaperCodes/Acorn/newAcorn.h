#ifndef __NEW_ACORN_H__
#define __NEW_ACORN_H__
#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"

void newupdate256(GRBModel& model, vector<GRBVar>& s, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep);
void setRoundFlags(int r, bitset<293>& nonZeroFlag, bitset<293>& varsFlag, vector<bitset<293>>& nonZeroTrackReps, vector<bitset<293>>& varsTrackReps);
STATUS newAcornForwardExpand(int rounds, bitset<128>& cube, vector<bitset<293>>& terms, float time, int thread);
STATUS newMidSolutionCounter(int rounds, bitset<293>& startcube, bitset<293 + 256>& last, map<bitset<293 + 256>,
    int, CMPS<293 + 256>> &counterMap, float time, int threads, int midround);
STATUS newMidSolutionCounter(int rounds, bitset<293>& startcube, bitset<293>& last, int& solcnt, float time, int threads);
#endif