#ifndef __GRAIN_H__
#define  __GRAIN_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"

void filterMap(map<bitset<256>, int, CMPS<256>>& mp);
int BackExpandPolynomial(int rounds, vector<bitset<256>>& term);

int SecondBackExpandPolynomial(int rounds, const bitset<256>& last,
	vector<bitset<256> >& term, int threads);
STATUS MidSolutionCounter(int rounds, const bitset<96>& cube, const
	bitset<256>& last, map<bitset<128>, int, CMPS<128>>& counterMap, float
	time, int thread);
STATUS MidSolutionCounter(int start, int end,  const bitset<256>& startstate,
	const bitset<256>& last, int& sols, float
	time, int thread);

void update(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& b);

#endif