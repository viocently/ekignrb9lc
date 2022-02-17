#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<algorithm>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"

using namespace std;

const int MAX = 200000000; // the maximum value of PoolSearchMode, P625

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
    GRBVar us##x1 = s[x1];\
    GRBVar us##x2 = s[x2];\
    GRBVar s##x1 = utap(model,s[x1]); \
    GRBVar s##x2 = utap(model,s[x2]); \
    uAnd2(model,s##x1,s##x2,o); \
    \
    vAnd(model,o,us##x1,s##x1); \
    vAnd(model,o,us##x2,s##x2); \
    \
    vars.emplace_back(o); \
}

#define IsNonZero1(x1) (nonZeroTrackRep[x1])
#define IsNonZero2(x1,x2) IsNonZero1(x1) && IsNonZero1(x2)

void uTriviumCore(GRBModel& model, vector<GRBVar>& s, int i1, int i2, int i3, int i4, int i5, bitset<288>& nonZeroTrackRep, bitset<288>& varsTrackRep)
{
    vector<GRBVar> vars;

    if (IsNonZero1(i1))
        AndOpr1(i1);

    if (IsNonZero1(i2))
        AndOpr1(i2);

    if (IsNonZero2(i3, i4))
        AndOpr2(i3, i4);

    if (IsNonZero1(i5))
        vars.emplace_back(s[i5]);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, O);
    else
        model.addConstr(O == 0);
    s[i5] = O;
}

void uupdate(GRBModel& model, vector<GRBVar>& s, bitset<288>& nonZeroTrackRep, bitset<288>& varsTrackRep)
{



    uTriviumCore(model, s, 65, 170, 90, 91, 92, nonZeroTrackRep, varsTrackRep);
    uTriviumCore(model, s, 161, 263, 174, 175, 176, nonZeroTrackRep, varsTrackRep);
    uTriviumCore(model, s, 242, 68, 285, 286, 287, nonZeroTrackRep, varsTrackRep);

    vector<GRBVar> tmp = s;
    for (int i = 0; i < 288; i++)
        s[(i + 1) % 288] = tmp[i];
}