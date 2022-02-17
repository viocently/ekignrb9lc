#include<iostream>
#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include<set>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"thread_pool.h"
#include"Acorntrack.h"


const int MAX = 200000000;

using namespace std;
using namespace thread_pool;

GRBVar utap(GRBModel& model, GRBVar& x)
{
    GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(x == y + z);
    x = z;
    return y;
}

void uAnd2(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _out)
{
    GRBVar vars[2] = { _in1,_in2 };
    model.addGenConstrOr(_out, vars, 2);
}

void uXor(GRBModel& model, vector<GRBVar> _ins, GRBVar _out)
{
    GRBLinExpr exp = 0;
    for (auto& x : _ins)
        exp += x;
    model.addConstr(exp == _out);
}

void vAnd(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _out)
{
    GRBVar vars[2] = { _in1,_in2 };
    model.addGenConstrAnd(_out, vars, 2);
}

// start macro definition
#define AndOpr1(x1) \
{ \
    GRBVar s##x1 = utap(model,s[x1]); \
    vars.emplace_back(s##x1); \
}

#define AndOpr2(x1,x2) \
{ \
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);\
    GRBVar s##x1 = utap(model,s[x1]); \
    GRBVar s##x2 = utap(model,s[x2]); \
    uAnd2(model,s##x1,s##x2,o); \
    \
    vAnd(model,o,varMap[x1],s##x1); \
    vAnd(model,o,varMap[x2],s##x2); \
    \
    vars.emplace_back(o); \
}
    /*// vAnd(model,o,us##x1,s##x1); \
    // vAnd(model,o,us##x2,s##x2); */

#define IsNonZero1(x1) (nonZeroTrackRep[x1])
#define IsNonZero2(x1,x2) IsNonZero1(x1) && IsNonZero1(x2)


static set<int> QuadExtractIndex(int i, int j, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    set<int> resSet;
    if (IsNonZero2(i, j))
    {
        if (varsTrackRep[i])
            resSet.emplace(i);
        if (varsTrackRep[j])
            resSet.emplace(j);
    }
    return resSet;
}

static set<int> LinearExtractIndex(int i, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    set<int> resSet;
    if (IsNonZero1(i))
        if (varsTrackRep[i])
            resSet.emplace(i);

    return resSet;
}


GRBVar uMaj(GRBModel& model, vector<GRBVar>& s, int i, int j, int k, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    map<int, GRBVar> varMap;
    varMap[i] = s[i];
    varMap[j] = s[j];
    varMap[k] = s[k];

    set<set<int>> rSets;
    set<int> tmpSet = QuadExtractIndex(i, j, nonZeroTrackRep, varsTrackRep);
    if (tmpSet.size() > 0) rSets.emplace(tmpSet);
    tmpSet = QuadExtractIndex(i, k, nonZeroTrackRep, varsTrackRep);
    if (tmpSet.size() > 0) rSets.emplace(tmpSet);
    tmpSet = QuadExtractIndex(j, k, nonZeroTrackRep, varsTrackRep);
    if (tmpSet.size() > 0) rSets.emplace(tmpSet);

    // remove redundance
    map<int, int> degMap;
    for (auto& x : rSets)
        for (auto& y : x)
            if (degMap[y] < x.size())
                degMap[y] = x.size();

    int size0 = rSets.size();
    set<set<int>> toDelete;
    for (auto& x : rSets)
        for (auto& y : x)
            if (degMap[y] > x.size())
            {
                toDelete.emplace(x);
                break;
            }

    for (auto& x : toDelete)
        rSets.erase(x);
    int size1 = rSets.size();
    cout << __func__ << " : " << size0 << "   " << size1 << endl;

    // build model
    vector<GRBVar> vars;
    for (auto& x : rSets)
        if (x.size() == 1)
        {
            int i = *x.begin();
            AndOpr1(i);
        }
        else // x == 2
        {
            int i = *x.begin();
            int j = *x.rbegin();
            GRBVar usi = s[i];
            GRBVar usj = s[j];
            AndOpr2(i, j);
        }
    /*vector<GRBVar> vars;
    if (IsNonZero2(i, j))
        AndOpr2(i, j);

    if (IsNonZero2(i, k))
        AndOpr2(i, k);

    if (IsNonZero2(k, j))
        AndOpr2(k, j);*/

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

GRBVar uCh(GRBModel& model, vector<GRBVar>& s, int i, int j, int k, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    map<int, GRBVar> varMap;
    varMap[i] = s[i];
    varMap[j] = s[j];
    varMap[k] = s[k];

    set<set<int>> rSets;
    set<int> tmpSet = QuadExtractIndex(i, j, nonZeroTrackRep, varsTrackRep);
    if (tmpSet.size() > 0) rSets.emplace(tmpSet);
    tmpSet = QuadExtractIndex(i, k, nonZeroTrackRep, varsTrackRep);
    if (tmpSet.size() > 0) rSets.emplace(tmpSet);
    tmpSet = LinearExtractIndex(k, nonZeroTrackRep, varsTrackRep);
    if (tmpSet.size() > 0) rSets.emplace(tmpSet);

    // remove redundance
    map<int, int> degMap;
    for (auto& x : rSets)
        for (auto& y : x)
            if (degMap[y] < x.size())
                degMap[y] = x.size();

    int size0 = rSets.size();
    set<set<int>> toDelete;
    for (auto& x : rSets)
        for (auto& y : x)
            if (degMap[y] > x.size())
            {
                toDelete.emplace(x);
                break;
            }

    for (auto& x : toDelete)
        rSets.erase(x);
    int size1 = rSets.size();
    cout << __func__ << " : " << size0 << "   " << size1 << endl;

    // build model
    vector<GRBVar> vars;
    for (auto& x : rSets)
        if (x.size() == 1)
        {
            int i = *x.begin();
            AndOpr1(i);
        }
        else // x == 2
        {
            int i = *x.begin();
            int j = *x.rbegin();
            AndOpr2(i, j);
        }
    /*vector<GRBVar> vars;
    if (IsNonZero2(i, j))
        AndOpr2(i, j);

    if (IsNonZero2(i, k))
        AndOpr2(i, k);

    if (IsNonZero1(k))
        AndOpr1(k);*/

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

void uXorF(GRBModel& model, vector<GRBVar>& s, int i, int j, int k, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    vector<GRBVar> vars;
    
    vars.emplace_back(s[i]);
    if (IsNonZero1(k))
        AndOpr1(k);
    if (IsNonZero1(j))
        AndOpr1(j);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, O);
    else
        model.addConstr(O == 0);

    s[i] = O;

    xorFTrack(nonZeroTrackRep, i, j, k);
    xorFVarsTrack(varsTrackRep, i, j, k);
}

GRBVar uKsg(GRBModel& model, vector<GRBVar>& s, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    vector<GRBVar> vars;
    if (IsNonZero1(12))
        AndOpr1(12);

    if (IsNonZero1(154))
        AndOpr1(154);

    GRBVar c = uMaj(model, s, 235, 61, 193, nonZeroTrackRep,varsTrackRep);
    GRBVar d = uCh(model, s, 230, 111, 66, nonZeroTrackRep,varsTrackRep);
    vars.emplace_back(c);
    vars.emplace_back(d);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

GRBVar uFbk(GRBModel& model, vector<GRBVar>& s, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    vector<GRBVar> vars;
    if (IsNonZero1(0))
        AndOpr1(0);

    if (IsNonZero1(107))
        AndOpr1(107);

    if (IsNonZero1(196))
        AndOpr1(196);

    GRBVar c = uMaj(model, s, 244, 23, 160, nonZeroTrackRep,varsTrackRep);
    vars.emplace_back(c);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

// we dont define genM here because we only use this after round 256
void uupdate(GRBModel& model, vector<GRBVar>& s, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    uXorF(model, s, 289, 235, 230, nonZeroTrackRep,varsTrackRep);
    uXorF(model, s, 230, 196, 193, nonZeroTrackRep, varsTrackRep);
    uXorF(model, s, 193, 160, 154, nonZeroTrackRep, varsTrackRep);
    uXorF(model, s, 154, 111, 107, nonZeroTrackRep, varsTrackRep);
    uXorF(model, s, 107, 66, 61, nonZeroTrackRep, varsTrackRep);
    uXorF(model, s, 61, 23, 0, nonZeroTrackRep, varsTrackRep);

    GRBVar e = uKsg(model, s, nonZeroTrackRep, varsTrackRep);
    GRBVar b = uFbk(model, s, nonZeroTrackRep, varsTrackRep);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(O == e + b);

    model.addConstr(s[0] == 0);

    for (int i = 0; i < 292; i++)
        s[i] = s[i + 1];
    s[292] = O;
}