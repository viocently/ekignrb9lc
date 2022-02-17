#ifndef __ACORN_TRACK_H__
#define __ACORN_TRACK_H__
#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"

void xorFTrack(bitset<293>& nonZeroTrack, int i, int j, int k);
void xorFVarsTrack(bitset<293>& varsTrack, int i, int j, int k);
int acornNonZeroTrack(vector<bitset<293>>& nonZeroTrackReps);
int acornVarsTrack0(vector<bitset<293>>& varsTrackReps, bitset<128>& cube, vector<bitset<293> >& nonZeroTrackReps);
int acornVarsTrack256(vector<bitset<293>>& varsTrackReps, bitset<293>& startcube, vector<bitset<293> >& nonZeroTrackReps);


#endif