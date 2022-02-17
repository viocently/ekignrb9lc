#ifndef __U_KREYVIUM_H__
#define  __U_KREYVIUM_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"


void uupdate(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& key, vector<GRBVar>& iv, bitset<288 + 256>& nonZeroTrackRep, bitset<288 + 256>& varsTrackRep);

#endif