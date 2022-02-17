#ifndef __TRIVIUM_CALLBACK_H__
#define  __TRIVIUM_CALLBACK_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
STATUS uMidSolutionCallBack(int rounds, bitset<80>& cube, const bitset<288>& last, vector<bitset<288> >& sols, float time, int threads, int midround);

#endif