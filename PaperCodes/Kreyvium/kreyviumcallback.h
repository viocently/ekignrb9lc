#ifndef __KREYVIUM_CALLBACK_H__
#define  __KREYVIUM_CALLBACK_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
STATUS MidSolutionCallBack(int rounds, bitset<128>& cube, const bitset<544>& last, vector<bitset<288 + 256> >& sols, float time, int threads, int midround);
STATUS uMidSolutionCallBack(int rounds, bitset<128>& cube, bitset<544>& last,
    vector<bitset<288 + 256>>& sols, float time, int threads, int midround);

#endif