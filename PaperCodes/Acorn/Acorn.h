#ifndef __ACORN_H__
#define __ACORN_H__
#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"

void update(GRBModel& model, vector<GRBVar>& x, vector<GRBVar>& k,
	vector<GRBVar>& v, int r);
int BackExpandPolynomial(int current, int rounds, vector<bitset<549>>& term);
int SecondBackExpandPolynomial(int current, int rounds, const bitset<549>& last,
	vector<bitset<549>>& term, int threads);
STATUS MidSolutionCounter(int midround, int rounds, const bitset<549>& start, const
	bitset<549>& last, int& solcnt, float
	time, int thread);
STATUS MidSolutionCounter(int midround, int rounds, const bitset<293>& startcube, const
	bitset<549>& last, map<bitset<549>, int, CMPS<549>>& counterMap, float
	time, int thread);
STATUS MidSolutionCounter(int rounds, const bitset<128>& cube, const
	bitset<549>& last, map<bitset<128>, int, CMPS<128>>& counterMap, float
	time, int thread);
STATUS MidSolutionCounter2(int rounds, const bitset<128>& cube, const
	bitset<549>& last, map<bitset<128>, int, CMPS<128>>& counterMap, float
	time, int thread);
#endif












