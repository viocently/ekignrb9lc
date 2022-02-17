#ifndef __KREYVIUM_H__
#define  __KREYVIUM_H__

#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h" 

void update(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& key, vector<GRBVar>& iv);

void filterMap(map<bitset<544>, int, CMPS<544>>& mp);
void filterTerms(const bitset<128>& cube, int rounds, vector<bitset<544>>& terms);
void BackExpandPolynomial( int rounds, vector<bitset<544>> & terms );

void SecondBackExpandPolynomial( int rounds, const bitset<288 + 256> & last,
		vector<bitset<288 + 256>> & terms, int threads );

STATUS MidSolutionCounter( int rounds, bitset<128> cube, const bitset<544> & last, map<bitset<128>, 
		int, CMPS<128>> & counterMap,  float time, int threads );

STATUS MidSolutionCounter(int start, int end, const bitset<256 + 288> startstate,
	const bitset<544>& last, int& sols, float time, int threads);

#endif
