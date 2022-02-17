#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include<set>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"grain.h"
#include"graintrack.h"

using namespace std;
const int MAX = 200000000; // the maximum value of PoolSearchMode, P625

GRBVar newtap(GRBModel& model, GRBVar& x) {

    GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar tmp[2] = { y,z };
    model.addGenConstrOr(x, tmp, 2);
    x = z;
    return y;
}

void newAnd2(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _out)
{
    GRBVar vars[2] = { _in1,_in2 };
    model.addGenConstrOr(_out, vars, 2);
}

void newAnd3(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _in3,GRBVar& _out)
{
    GRBVar vars[3] = { _in1,_in2,_in3 };
    model.addGenConstrOr(_out, vars, 3);
}

void newAnd4(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _in3, GRBVar& _in4,GRBVar& _out)
{
    GRBVar vars[4] = { _in1,_in2,_in3,_in4 };
    model.addGenConstrOr(_out, vars, 4);
}

void newXor(GRBModel& model, vector<GRBVar> _ins, GRBVar _out)
{
    GRBLinExpr exp = 0;
    for (auto& x : _ins)
        exp += x;
    model.addConstr(exp == _out);
}


// start macro definition
#define AndOpr1(s1,x1) \
{ \
    GRBVar s1##x1 = newtap(model,s1[x1]); \
    vars.emplace_back(s1##x1); \
}

#define AndOpr2(s1,x1,s2,x2) \
{ \
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);\
    GRBVar s1##x1 = newtap(model,s1[x1]); \
    GRBVar s2##x2 = newtap(model,s2[x2]); \
    newAnd2(model,s1##x1,s2##x2,o); \
    \
    if(varsTrackRep[#s1 == "b" ? x1 : x1+128]) model.addConstr(o==s1##x1); \
    if(varsTrackRep[#s2 == "b" ? x2 : x2+128]) model.addConstr(o==s2##x2); \
    \
    vars.emplace_back(o); \
} 

#define AndOpr3(s1,x1,s2,x2,s3,x3) \
{ \
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);\
    GRBVar s1##x1 = newtap(model,s1[x1]); \
    GRBVar s2##x2 = newtap(model,s2[x2]); \
    GRBVar s3##x3 = newtap(model,s3[x3]); \
    newAnd3(model,s1##x1,s2##x2,s3##x3,o); \
    \
    if(varsTrackRep[#s1 == "b" ? x1 : x1+128]) model.addConstr(o==s1##x1); \
    if(varsTrackRep[#s2 == "b" ? x2 : x2+128]) model.addConstr(o==s2##x2); \
    if(varsTrackRep[#s3 == "b" ? x3 : x3+128]) model.addConstr(o==s3##x3); \
    \
    vars.emplace_back(o); \
} 

#define AndOpr4(s1,x1,s2,x2,s3,x3,s4,x4) \
{ \
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);\
    GRBVar s1##x1 = newtap(model,s1[x1]); \
    GRBVar s2##x2 = newtap(model,s2[x2]); \
    GRBVar s3##x3 = newtap(model,s3[x3]); \
    GRBVar s4##x4 = newtap(model,s4[x4]); \
    newAnd4(model,s1##x1,s2##x2,s3##x3,s4##x4,o); \
    \
    if(varsTrackRep[#s1 == "b" ? x1 : x1+128]) model.addConstr(o==s1##x1); \
    if(varsTrackRep[#s2 == "b" ? x2 : x2+128]) model.addConstr(o==s2##x2); \
    if(varsTrackRep[#s3 == "b" ? x3 : x3+128]) model.addConstr(o==s3##x3); \
    if(varsTrackRep[#s4 == "b" ? x4 : x4+128]) model.addConstr(o==s4##x4); \
    \
    vars.emplace_back(o); \
} 

#define IsNonZero1(s1,x1) (nonZeroTrackRep[#s1 == "b" ? x1 : x1+128])
#define IsNonZero2(s1,x1,s2,x2) IsNonZero1(s1,x1) && IsNonZero1(s2,x2)
#define IsNonZero3(s1,x1,s2,x2,s3,x3) IsNonZero1(s1,x1) && IsNonZero2(s2,x2,s3,x3)
#define IsNonZero4(s1,x1,s2,x2,s3,x3,s4,x4) IsNonZero2(s1,x1,s2,x2) && IsNonZero2(s3,x3,s4,x4) 

GRBVar newfuncH(GRBModel& model, vector<GRBVar>& b, vector<GRBVar>& s, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
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
        newXor(model, vars, H);
    else
        model.addConstr(H == 0);

    return H;
}

GRBVar newfuncO(GRBModel& model, vector<GRBVar>& b, vector<GRBVar>& s, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
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
        newXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

GRBVar newfuncF(GRBModel& model, vector<GRBVar>& s, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
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
        newXor(model, vars, F);
    else
        model.addConstr(F == 0);

    return F;
}

GRBVar newfuncG(GRBModel& model, vector<GRBVar>& b, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
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
        newXor(model, vars, G);
    else
        model.addConstr(G == 0);

    return G;
}

void newupdate(GRBModel& model, vector<GRBVar>& s, vector<GRBVar>& b, bitset<256>& nonZeroTrackRep, bitset<256>& varsTrackRep)
{
    GRBVar s0 = s[0];
    GRBVar h = newfuncH(model, b, s, nonZeroTrackRep, varsTrackRep);
    GRBVar o = newfuncO(model, b, s, nonZeroTrackRep, varsTrackRep);

    GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(z == h + o);

    GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar tmpVar[2] = { z1, z2 };
    model.addGenConstrOr(z, tmpVar, 2);

    GRBVar f = newfuncF(model, s, nonZeroTrackRep, varsTrackRep);
    GRBVar g = newfuncG(model, b, nonZeroTrackRep, varsTrackRep);

    GRBVar news = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(news == z1 + f);

    GRBVar newb = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(newb == z2 + g + s[0]);


    for (int i = 0; i < 127; i++)
    {
        b[i] = b[i + 1];
        s[i] = s[i + 1];
    }
    b[127] = newb;
    s[127] = news;

    // remove (s[r][0], z[r]) = (1,1) 
    model.addConstr((1 - s0) + (1 - z) >= 1);
}



STATUS newMidSolutionCounter(int rounds, bitset<96>& cube,
    bitset<256>& last, map<bitset<256>, int, CMPS<256>>& counterMap, float
    time, int thread, int midround)
{
    if (midround > rounds || rounds < 0 || midround < 0)
    {
        logger(__func__ + string(":parameter error."));
        exit(-1);
    }


    bitset<256> flag1;
    flag1.set();

    counterMap.clear();


    vector<bitset<256> > nonZeroTrackReps;
    grainTrack(cube, nonZeroTrackReps);


    vector<bitset<256> > varsTrackReps;
    grainVarsTrack(cube, varsTrackReps, nonZeroTrackReps);

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    //env.set(GRB_IntParam_Presolve, 2);
    //env.set(GRB_IntParam_MIPFocus, 1);
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> s(128);
    vector<GRBVar> b(128);
    for (int i = 0; i < 128; i++)
    {
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
        b[i] = model.addVar(0, 1, 0, GRB_BINARY);
    }

   // constraints
    for (int i = 0; i < 96; i++)
        if (cube[i] == 1)
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);

    for (int i = 96; i < 128; i++)
        model.addConstr(s[i] == 0);

    for (int i = 0; i < 128; i++)
        model.addConstr(b[i] == 0);

    // Round function
    vector<GRBVar> workb = b;
    vector<GRBVar> works = s;

    for (int r = 0; r < midround; r++)
    {
        newupdate(model, works, workb, r < nonZeroTrackReps.size() ? nonZeroTrackReps[r] : flag1, r < varsTrackReps.size() ? varsTrackReps[r] : flag1);
    }

    vector<GRBVar> midworks(128);
    vector<GRBVar> midworkb(128);

    vector<GRBVar> midworks2 = works;
    vector<GRBVar> midworkb2 = workb;

    // mid constraints
    bitset<256> midNonZeroFlag = (midround < nonZeroTrackReps.size()) ? nonZeroTrackReps[midround] : flag1;
    bitset<256> midVarsFlag = (midround < varsTrackReps.size()) ? varsTrackReps[midround] : flag1;
    
    GRBLinExpr nk = 0;
    for (int i = 0; i < 128; i++)
    {
        midworkb[i] = model.addVar(0, 1, 0, GRB_BINARY);
        if (midVarsFlag[i])
            model.addConstr(midworkb[i] == workb[i]);
        else
            if (midNonZeroFlag[i])
            {
                model.addConstr(midworkb[i] >= workb[i]);
                nk += midworkb[i] - workb[i];
            }
            else
                model.addConstr(midworkb[i] == 0);
    }

    for (int i = 0; i < 128; i++)
    {
        midworks[i] = model.addVar(0, 1, 0, GRB_BINARY);
        if (midVarsFlag[i + 128])
            model.addConstr(midworks[i] == works[i]);
        else
            if (midNonZeroFlag[i + 128])
            {
                model.addConstr(midworks[i] >= works[i]);
                nk += midworks[i] - works[i];
            }
            else
                model.addConstr(midworks[i] == 0);
    }

    works = midworks;
    workb = midworkb;

    // second Phase
    for (int r = midround; r < rounds; r++)
        update(model, works, workb);

    // output constraints
    for (int i = 0; i < 128; i++)
        if (last[i] == 1)
            model.addConstr(workb[i] == 1);
        else
            model.addConstr(workb[i] == 0);

    for (int i = 0; i < 128; i++)
        if (last[i + 128] == 1)
            model.addConstr(works[i] == 1);
        else
            model.addConstr(works[i] == 0);

    if (time > 0)
        model.set(GRB_DoubleParam_TimeLimit, time);

    model.setObjective(nk, GRB_MAXIMIZE);

    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
    {
        logger(__func__ + string(" : ") + to_string(rounds) +
            string(" | EXPAND "));
        return EXPAND;
    }
    else
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
        if (solCount >= MAX)
        {
            cerr << "solCount value  is too big !" << endl;
            exit(-1);
        }

        double time = model.get(GRB_DoubleAttr_Runtime);
        logger(__func__ + string(" : ") + to_string(rounds) + string(" | "
        ) + to_string(time) + string(" | "
        ) + to_string(solCount));


        if (solCount > 0)
        {
            bitset<256> midState;
            bitset<256> midState2;
            for (int i = 0; i < solCount; i++)
            {
                model.set(GRB_IntParam_SolutionNumber, i);
                for (int j = 0; j < 128; j++)
                {
                    if (round(midworkb[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState[j] = 1;
                    else
                        midState[j] = 0;
                    if (round(midworks[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState[j+128] = 1;
                    else
                        midState[j+128] = 0;
                }

                counterMap[midState]++;
            }
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }
}

STATUS newMidSolutionCounter(int rounds, bitset<96>& cube,
    bitset<256>& last, int & sols, float
    time, int thread)
{


    bitset<256> flag1;
    flag1.set();



    vector<bitset<256> > nonZeroTrackReps;
    grainTrack(cube, nonZeroTrackReps);


    vector<bitset<256> > varsTrackReps;
    grainVarsTrack(cube, varsTrackReps, nonZeroTrackReps);

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_MIPFocus, 1);
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> s(128);
    vector<GRBVar> b(128);
    for (int i = 0; i < 128; i++)
    {
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
        b[i] = model.addVar(0, 1, 0, GRB_BINARY);
    }

    // constraints
    for (int i = 0; i < 96; i++)
        if (cube[i] == 1)
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);

    for (int i = 96; i < 128; i++)
        model.addConstr(s[i] == 0);

    for (int i = 0; i < 128; i++)
        model.addConstr(b[i] == 0);

    // Round function
    vector<GRBVar> workb = b;
    vector<GRBVar> works = s;

    for (int r = 0; r < rounds; r++)
    {
        newupdate(model, works, workb, r < nonZeroTrackReps.size() ? nonZeroTrackReps[r] : flag1, r < varsTrackReps.size() ? varsTrackReps[r] : flag1);
    }

    bitset<256> lastNonZeroFlag, lastVarsFlag;
    lastNonZeroFlag = rounds < nonZeroTrackReps.size() ? nonZeroTrackReps[rounds] : flag1;
    lastVarsFlag = rounds < varsTrackReps.size() ? varsTrackReps[rounds] : flag1;

    // output constraints
    for (int i = 0; i < 128; i++)
        if (last[i] == 1 && lastVarsFlag[i])
            model.addConstr(workb[i] == 1);
        else
            model.addConstr(workb[i] == 0);

    for (int i = 0; i < 128; i++)
        if (last[i + 128] == 1 && lastVarsFlag[i+128])
            model.addConstr(works[i] == 1);
        else
            model.addConstr(works[i] == 0);

    if (time > 0)
        model.set(GRB_DoubleParam_TimeLimit, time);


    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
    {
        return EXPAND;
    }
    else
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
        if (solCount >= MAX)
        {
            cerr << "solCount value  is too big !" << endl;
            exit(-1);
        }

        sols = solCount;

        if (solCount > 0)
        {
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }
}



