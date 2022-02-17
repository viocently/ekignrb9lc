#ifndef __U_ACORN_H__
#define __U_ACORN_H__
#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"

void uupdate(GRBModel& model, vector<GRBVar>& s, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep);

#endif