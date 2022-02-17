#ifndef __TRIVIUM_TRACK_H__
#define  __TRIVIUM_TRACK_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"

void triviumTrack(bitset<80>& cube, vector<bitset<288>>& trackReps);
void triviumVarsTrack(bitset<80>& cube, vector<bitset<288>>& trackReps, vector<bitset<288>>& nonZeros);
#endif