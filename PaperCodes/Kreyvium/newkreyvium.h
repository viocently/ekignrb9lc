#ifndef __NEW_KREYVIUM_H__
#define  __NEW_KREYVIUM_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

void newupdate(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& key, vector<GRBVar>& iv, bitset<288 + 256>& nonZeroTrackRep, bitset<288 + 256>& varsTrackRep);
void setRoundFlags(int r, bitset<128>& cube, bitset<288 + 256>& nonZeroFlag, bitset<288 + 256>& varsFlag, vector<bitset<288 + 256>>& nonZeroTrackReps, vector<bitset<288 + 256>>& varsTrackReps);
STATUS newMidSolutionCounter(int rounds, bitset<128>& cube, bitset<544>& last, map<bitset<288 + 256>,
    vector<bitset<288 + 256>>, CMPS<288 + 256>> &counterMap, float time, int threads, int midround);

STATUS newMidSolutionCounter(int rounds, bitset<128>& cube, bitset<544>& last,
    int& sols, float time, int threads);
#endif