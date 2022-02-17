#ifndef __GRAIN_CALLBACK_H__
#define  __GRAIN_CALLBACK_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

STATUS uMidSolutionCallBack(int rounds, bitset<96>& cube,
    bitset<256>& last, int midround, vector<bitset<256>>& sols, float
    time, int thread, int coRound);
#endif