#ifndef __GRAIN_TRACK_H__
#define  __GRAIN_TRACK_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

void grainTrack(bitset<96>& cube, vector<bitset<256>>& trackReps);
void grainVarsTrack(bitset<96>& cube, vector<bitset<256>>& trackReps, vector<bitset<256>>& nonZeros);
#endif