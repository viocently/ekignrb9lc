#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<algorithm>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"newkreyvium.h"

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

GRBVar uLFSR(GRBModel& model, vector <GRBVar>& x)
{
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);

    model.addConstr(x[0] == a + b);

    vector<GRBVar> t = x;
    for (int i = 0; i < 127; i++)
        x[i] = t[i + 1];
    x[127] = a;
    return b;
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

void uKreyviumCore(GRBModel& model, vector<GRBVar>& s, int i1, int i2, int i3, int i4, int i5, bitset<288 + 256>& nonZeroTrackRep, bitset<288 + 256>& varsTrackRep)
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

void uupdate(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& key, vector<GRBVar>& iv, bitset<288 + 256>& nonZeroTrackRep, bitset<288 + 256>& varsTrackRep)
{
    // k dont need to LFSR because its all 0 ,so the k in the register has no contribution to the propagation
    GRBVar genIV;
    if (nonZeroTrackRep[288] == 1)
        genIV = uLFSR(model, iv);
    else
    {
        genIV = model.addVar(0, 1, 0, GRB_BINARY);
        model.addConstr(genIV == 0);
        vector<GRBVar> tmpiv = iv;
        for (int i = 0; i < 128; i++)
            iv[i] = tmpiv[(i + 1) % 128];
    }

    vector<GRBVar> tmpkey = key;
    for (int i = 0; i < 128; i++)
        key[i] = tmpkey[(i + 1) % 128];

    uKreyviumCore(model, s, 65, 170, 90, 91, 92, nonZeroTrackRep, varsTrackRep);
    uKreyviumCore(model, s, 161, 263, 174, 175, 176, nonZeroTrackRep, varsTrackRep);
    uKreyviumCore(model, s, 242, 68, 285, 286, 287, nonZeroTrackRep, varsTrackRep);

    GRBVar t1 = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(t1 == s[92] + genIV);

    s[92] = t1;
    vector<GRBVar> tmp = s;
    for (int i = 0; i < 288; i++)
        s[(i + 1) % 288] = tmp[i];
}

