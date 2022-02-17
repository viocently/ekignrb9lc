#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<algorithm>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"triviumtrack.h"
#include"trivium.h"

using namespace std;

const int MAX = 200000000; // the maximum value of PoolSearchMode, P625

void newAnd2(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _out)
{
    GRBVar vars[2] = { _in1,_in2 };
    model.addGenConstrOr(_out, vars, 2);
}

void newXor(GRBModel& model, vector<GRBVar> _ins, GRBVar _out)
{
    GRBLinExpr exp = 0;
    for (auto& x : _ins)
        exp += x;
    model.addConstr(exp == _out);
}

GRBVar newtap(GRBModel& model, GRBVar& x) {

    GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar tmp[2] = { y,z };
    model.addGenConstrOr(x, tmp, 2);
    x = z;
    return y;
}

// start macro definition
#define AndOpr1(x1) \
{ \
    GRBVar s##x1 = newtap(model,s[x1]); \
    vars.emplace_back(s##x1); \
}

#define AndOpr2(x1,x2) \
{ \
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);\
    GRBVar s##x1 = newtap(model,s[x1]); \
    GRBVar s##x2 = newtap(model,s[x2]); \
    newAnd2(model,s##x1,s##x2,o); \
    \
    if(varsTrackRep[x1]) model.addConstr(o==s##x1); \
    if(varsTrackRep[x2]) model.addConstr(o==s##x2); \
    \
    vars.emplace_back(o); \
}

#define IsNonZero1(x1) (nonZeroTrackRep[x1])
#define IsNonZero2(x1,x2) IsNonZero1(x1) && IsNonZero1(x2)

void newTriviumCore(GRBModel& model, vector<GRBVar>& s, int i1, int i2, int i3, int i4, int i5, bitset<288>& nonZeroTrackRep, bitset<288>& varsTrackRep)
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
        newXor(model, vars, O);
    else
        model.addConstr(O == 0);
    s[i5] = O;
}

void newupdate(GRBModel& model, vector<GRBVar>& s, bitset<288>& nonZeroTrackRep, bitset<288>& varsTrackRep)
{

    newTriviumCore(model, s, 65, 170, 90, 91, 92, nonZeroTrackRep, varsTrackRep);
    newTriviumCore(model, s, 161, 263, 174, 175, 176, nonZeroTrackRep, varsTrackRep);
    newTriviumCore(model, s, 242, 68, 285, 286, 287, nonZeroTrackRep, varsTrackRep);

    vector<GRBVar> tmp = s;
    for (int i = 0; i < 288; i++)
        s[(i + 1) % 288] = tmp[i];
}

void setRoundFlags(int r, bitset<80>& cube, bitset<288>& nonZeroFlag, bitset<288>& varsFlag, vector<bitset<288>>& nonZeroTrackReps, vector<bitset<288>>& varsTrackReps)
{
    nonZeroFlag.reset();
    varsFlag.reset();

    if (r < nonZeroTrackReps.size())
        nonZeroFlag = nonZeroTrackReps[r];
    else
    {
        for (int i = 0; i < 288; i++)
            nonZeroFlag[i] = 1;
    }

    if (r < varsTrackReps.size())
        varsFlag = varsTrackReps[r];
    else
    {
        for (int i = 0; i < 288; i++)
            varsFlag[i] = 1;
    }
}

STATUS newMidSolutionCounter(int rounds, bitset<80>& cube, bitset<288>& last, map<bitset<288>,
    int, CMPS<288>> &counterMap, float time, int threads, int midround)
{
    if (midround > rounds || rounds < 0 || midround < 0)
    {
        logger(__func__ + string(":parameter error."));
        exit(-1);
    }

    counterMap.clear();


    vector<bitset<288>> nonZeroTrackReps;
    triviumTrack(cube, nonZeroTrackReps);

    vector<bitset<288>> varsTrackReps;
    triviumVarsTrack(cube, varsTrackReps, nonZeroTrackReps);


    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, MAX);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_MIPFocus, 1);

    // Create the model
    GRBModel model = GRBModel(env);

    // register
    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
   

    // constraints
    for (int i = 0; i < 80; i++)
        if (cube[i] == 0)
            model.addConstr(s[i + 93] == 0);
        else
            model.addConstr(s[i + 93] == 1);

    for (int i = 0; i < 93; i++)
        model.addConstr(s[i] == 0);
    for (int i = 93 + 80; i < 288; i++)
        model.addConstr(s[i] == 0);

    // update
    for (int r = 0; r < midround; r++)
    {
        bitset<288> nonZeroFlag;
        bitset<288> varsFlag;
        setRoundFlags(r, cube, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);
        newupdate(model, s, nonZeroFlag, varsFlag);
    }

    // mid constraints
    bitset<288> midNonZeroFlag;
    bitset<288> midVarsFlag;
    setRoundFlags(midround, cube, midNonZeroFlag, midVarsFlag, nonZeroTrackReps, varsTrackReps);


    vector<GRBVar> mids2 = s;

    GRBLinExpr nk = 0;
    vector<GRBVar> mids(288);
    for (int i = 0; i < 288; i++)
        mids[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for (int i = 0; i < 288; i++)
        if (midVarsFlag[i])
            model.addConstr(mids[i] == s[i]);
        else
            if (midNonZeroFlag[i])
            {
                model.addConstr(mids[i] >= s[i]);
                nk += (mids[i] - s[i]);
            }
            else
                model.addConstr(mids[i] == 0);

    s = mids;

    // second Phase
    for (int r = midround; r < rounds; r++)
        update(model, s);

    // output constraints
    // key + iv + 288
    for (int i = 0; i < 288; i++)
        if (last[i])
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);

    // objective function
    model.setObjective(nk, GRB_MAXIMIZE);

    if (time > 0)
        model.set(GRB_DoubleParam_TimeLimit, time);

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
            // k + iv +288
            bitset<288> midState;
            for (int i = 0; i < solCount; i++)
            {
                model.set(GRB_IntParam_SolutionNumber, i);

                // midstate
                for (int j = 0; j < 288; j++)
                    if (round(mids[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState[j] = 1;
                    else
                        midState[j] = 0;

                counterMap[midState]++;
            }
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }


}

STATUS newMidSolutionCounter(int rounds, bitset<80>& cube, const bitset<288>& last,
    int& sols, float time, int threads)
{

  
    vector<bitset<288>> nonZeroTrackReps;
    triviumTrack(cube, nonZeroTrackReps);


    vector<bitset<288>> varsTrackReps;
    triviumVarsTrack(cube, varsTrackReps, nonZeroTrackReps);


    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, MAX);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_MIPFocus, 1);

    // Create the model
    GRBModel model = GRBModel(env);

    // register
    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
   

    // constraints
    for (int i = 0; i < 80; i++)
        if (cube[i] == 0)
            model.addConstr(s[i + 93] == 0);
        else
            model.addConstr(s[i + 93] == 1);

    for (int i = 0; i < 93; i++)
        model.addConstr(s[i] == 0);
    for (int i = 93 + 80; i < 288; i++)
        model.addConstr(s[i] == 0);

    // update
    for (int r = 0; r < rounds; r++)
    {
        bitset<288> nonZeroFlag;
        bitset<288> varsFlag;
        setRoundFlags(r, cube, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);
        newupdate(model, s, nonZeroFlag, varsFlag);
    }

    // output constraints
    // key + iv + 288
    bitset<288> finalNonZeroFlag;
    bitset<288> finalVarsFlag;
    setRoundFlags(rounds, cube, finalNonZeroFlag, finalVarsFlag, nonZeroTrackReps, varsTrackReps);

    for (int i = 0; i < 288; i++)
        if (finalVarsFlag[i] && last[i] == 1)
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);


    if (time > 0)
        model.set(GRB_DoubleParam_TimeLimit, time);

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

        sols = solCount;

        if (solCount > 0)
        {
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }


}
