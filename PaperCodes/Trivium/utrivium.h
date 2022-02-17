#ifndef __U_TRIVIUM_H__
#define  __U_TRIVIUM_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"


void uupdate(GRBModel& model, vector<GRBVar>& s, bitset<288>& nonZeroTrackRep, bitset<288>& varsTrackRep);

#endif