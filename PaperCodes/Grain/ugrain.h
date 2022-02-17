#ifndef __U_GRAIN_H__
#define  __U_GRAIN_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"


void uupdate(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& b, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep);

#endif