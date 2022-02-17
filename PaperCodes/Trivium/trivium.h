#pragma once
#ifndef __TRIVIUM_H__
#define __TRIVIUM_H__
#include<vector>
#include<bitset>
#include<map>
#include"node.h"
#include"gurobi_c++.h"


void triviumCore(GRBModel&, vector<GRBVar>&, int, int, int, int, int);
void update(GRBModel& model, vector<GRBVar>& s);
void filterMap(map<bitset<288>, int, CMPS<288>>& mp);
void filterTerms(const bitset<80>& cube, int rounds, vector<bitset<288>>& terms);
int SecondBackExpandPolynomial(int, const bitset<288>&, vector<bitset<288> >&, int);

STATUS MidSolutionCounter(int, const bitset<80>&, const bitset<288>&,
    map<bitset<80>, int, CMPS<80>>&, float, int);

STATUS MidSolutionCounter(int start, int end, const bitset<288>& startstate, const
    bitset<288>& laststate, int& solcnt, float
    time, int thread);
int BackExpandPolynomial(int, vector<bitset<288> >&);

#endif