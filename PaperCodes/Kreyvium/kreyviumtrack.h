#ifndef __KREYVIUM_TRACK_H__
#define  __KREYVIUM_TRACK_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

void kreyviumTrack(bitset<128>& cube, vector<bitset<288 + 256>>& trackReps);
void kreyviumVarsTrack(bitset<128>& cube, vector<bitset<288 + 256>>& trackReps, vector<bitset<288 + 256>>& nonZeros);
#endif