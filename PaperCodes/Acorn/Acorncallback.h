#ifndef __ACORN_CALLBACK_H__
#define  __ACORN_CALLBACK_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

STATUS uMidSolutionCallBack(int rounds, bitset<293>& startcube, bitset<293 + 256>& last, vector<bitset<293 + 256>>& sols, float time, int threads, int midround, int coRound);

#endif