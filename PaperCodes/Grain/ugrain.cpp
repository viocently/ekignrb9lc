#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"grain.h"
#include"graintrack.h"

using namespace std;
const int MAX = 200000000; // the maximum value of PoolSearchMode, P625

GRBVar utap(GRBModel& model, GRBVar& x) {

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

void uAnd3(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _in3, GRBVar& _out)
{
    GRBVar vars[3] = { _in1,_in2,_in3 };
    model.addGenConstrOr(_out, vars, 3);
}

void uAnd4(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _in3, GRBVar& _in4, GRBVar& _out)
{
    GRBVar vars[4] = { _in1,_in2,_in3,_in4 };
    model.addGenConstrOr(_out, vars, 4);
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

#define AndOpr1(s1,x1) \
{ \
    GRBVar s1##x1 = utap(model,s1[x1]); \
    vars.emplace_back(s1##x1); \
}

#define AndOpr2(s1,x1,s2,x2) \
{ \
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);\
    GRBVar u##s1##x1 = s1[x1];\
    GRBVar u##s2##x2 = s2[x2];\
    GRBVar s1##x1 = utap(model,s1[x1]); \
    GRBVar s2##x2 = utap(model,s2[x2]); \
    uAnd2(model,s1##x1,s2##x2,o); \
    \
    vAnd(model,o,u##s1##x1,s1##x1); \
    vAnd(model,o,u##s2##x2,s2##x2); \
    \
    vars.emplace_back(o); \
} 

#define AndOpr3(s1,x1,s2,x2,s3,x3) \
{ \
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);\
    GRBVar u##s1##x1 = s1[x1];\
    GRBVar u##s2##x2 = s2[x2];\
    GRBVar u##s3##x3 = s3[x3];\
    GRBVar s1##x1 = utap(model,s1[x1]); \
    GRBVar s2##x2 = utap(model,s2[x2]); \
    GRBVar s3##x3 = utap(model,s3[x3]); \
    uAnd3(model,s1##x1,s2##x2,s3##x3,o); \
    \
    vAnd(model,o,u##s1##x1,s1##x1); \
    vAnd(model,o,u##s2##x2,s2##x2); \
    vAnd(model,o,u##s3##x3,s3##x3); \
    \
    vars.emplace_back(o); \
} 

#define AndOpr4(s1,x1,s2,x2,s3,x3,s4,x4) \
{ \
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);\
    GRBVar u##s1##x1 = s1[x1];\
    GRBVar u##s2##x2 = s2[x2];\
    GRBVar u##s3##x3 = s3[x3];\
    GRBVar u##s4##x4 = s4[x4];\
    GRBVar s1##x1 = utap(model,s1[x1]); \
    GRBVar s2##x2 = utap(model,s2[x2]); \
    GRBVar s3##x3 = utap(model,s3[x3]); \
    GRBVar s4##x4 = utap(model,s4[x4]); \
    uAnd4(model,s1##x1,s2##x2,s3##x3,s4##x4,o); \
    \
    vAnd(model,o,u##s1##x1,s1##x1); \
    vAnd(model,o,u##s2##x2,s2##x2); \
    vAnd(model,o,u##s3##x3,s3##x3); \
    vAnd(model,o,u##s4##x4,s4##x4); \
    \
    vars.emplace_back(o); \
} 

#define IsNonZero1(s1,x1) (nonZeroTrackRep[#s1 == "b" ? x1 : x1+128])
#define IsNonZero2(s1,x1,s2,x2) IsNonZero1(s1,x1) && IsNonZero1(s2,x2)
#define IsNonZero3(s1,x1,s2,x2,s3,x3) IsNonZero1(s1,x1) && IsNonZero2(s2,x2,s3,x3)
#define IsNonZero4(s1,x1,s2,x2,s3,x3,s4,x4) IsNonZero2(s1,x1,s2,x2) && IsNonZero2(s3,x3,s4,x4) 

GRBVar ufuncH(GRBModel& model, vector<GRBVar>& b, vector<GRBVar>& s, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
{
    vector<GRBVar> vars;

    if (IsNonZero2(b, 12, s, 8))
        AndOpr2(b, 12, s, 8);

    if (IsNonZero2(s, 13, s, 20))
        AndOpr2(s, 13, s, 20);

    if (IsNonZero2(b, 95, s, 42))
        AndOpr2(b, 95, s, 42);

    if (IsNonZero2(s, 60, s, 79))
        AndOpr2(s, 60, s, 79);

    if (IsNonZero3(b, 12, b, 95, s, 94))
        AndOpr3(b, 12, b, 95, s, 94);

    GRBVar H = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, H);
    else
        model.addConstr(H == 0);

    return H;
}

GRBVar ufuncO(GRBModel& model, vector<GRBVar>& b, vector<GRBVar>& s, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
{
    vector<GRBVar> vars;

    if (IsNonZero1(s, 93))
        AndOpr1(s, 93);

    if (IsNonZero1(b, 2))
        AndOpr1(b, 2);

    if (IsNonZero1(b, 15))
        AndOpr1(b, 15);

    if (IsNonZero1(b, 36))
        AndOpr1(b, 36);

    if (IsNonZero1(b, 45))
        AndOpr1(b, 45);

    if (IsNonZero1(b, 64))
        AndOpr1(b, 64);

    if (IsNonZero1(b, 73))
        AndOpr1(b, 73);

    if (IsNonZero1(b, 89))
        AndOpr1(b, 89);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

GRBVar ufuncF(GRBModel& model, vector<GRBVar>& s, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
{
    vector<GRBVar> vars;

    if (IsNonZero1(s, 0))
        AndOpr1(s, 0);

    if (IsNonZero1(s, 7))
        AndOpr1(s, 7);

    if (IsNonZero1(s, 38))
        AndOpr1(s, 38);

    if (IsNonZero1(s, 70))
        AndOpr1(s, 70);

    if (IsNonZero1(s, 81))
        AndOpr1(s, 81);

    if (IsNonZero1(s, 96))
        AndOpr1(s, 96);

    GRBVar F = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, F);
    else
        model.addConstr(F == 0);

    return F;
}

GRBVar ufuncG(GRBModel& model, vector<GRBVar>& b, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
{
    vector<GRBVar> vars;

    if (IsNonZero1(b, 0))
        vars.emplace_back(b[0]);

    if (IsNonZero1(b, 26))
        AndOpr1(b, 26);

    if (IsNonZero1(b, 56))
        AndOpr1(b, 56);

    if (IsNonZero1(b, 91))
        AndOpr1(b, 91);

    if (IsNonZero1(b, 96))
        AndOpr1(b, 96);

    if (IsNonZero2(b, 3, b, 67))
        AndOpr2(b, 3, b, 67);

    if (IsNonZero2(b, 11, b, 13))
        AndOpr2(b, 11, b, 13);

    if (IsNonZero2(b, 17, b, 18))
        AndOpr2(b, 17, b, 18);

    if (IsNonZero2(b, 27, b, 59))
        AndOpr2(b, 27, b, 59);

    if (IsNonZero2(b, 40, b, 48))
        AndOpr2(b, 40, b, 48);

    if (IsNonZero2(b, 61, b, 65))
        AndOpr2(b, 61, b, 65);

    if (IsNonZero2(b, 68, b, 84))
        AndOpr2(b, 68, b, 84);

    if (IsNonZero4(b, 88, b, 92, b, 93, b, 95))
        AndOpr4(b, 88, b, 92, b, 93, b, 95);

    if (IsNonZero3(b, 22, b, 24, b, 25))
        AndOpr3(b, 22, b, 24, b, 25);

    if (IsNonZero3(b, 70, b, 78, b, 82))
        AndOpr3(b, 70, b, 78, b, 82);

    GRBVar G = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        uXor(model, vars, G);
    else
        model.addConstr(G == 0);

    return G;
}

void uupdate(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& b, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
{
    GRBVar s0 = s[0];
    GRBVar h = ufuncH(model, b, s, nonZeroTrackRep, varsTrackRep);
    GRBVar o = ufuncO(model, b, s, nonZeroTrackRep, varsTrackRep);

    GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(z == h + o);

    GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar tmpVar[2] = { z1, z2 };
    model.addConstr(z == z1 + z2);

    GRBVar f = ufuncF(model, s, nonZeroTrackRep, varsTrackRep);
    GRBVar g = ufuncG(model, b, nonZeroTrackRep, varsTrackRep);

    GRBVar us = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(us == z1 + f);

    GRBVar ub = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(ub == z2 + g + s[0]);


    for (int i = 0; i < 127; i++)
    {
        b[i] = b[i + 1];
        s[i] = s[i + 1];
    }
    b[127] = ub;
    s[127] = us;

    // remove (s[r][0], z[r]) = (1,1) 
    model.addConstr((1 - s0) + (1 - z) >= 1);
}