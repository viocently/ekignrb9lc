#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include<array>
#include<thread>
#include<set>
#include<algorithm>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"thread_pool.h"
#include"Acorntrack.h"
#include"Acorn.h"

const int MAX = 200000000;

using namespace std;
using namespace thread_pool;

void newAnd2(GRBModel& model, GRBVar& _in1, GRBVar& _in2, GRBVar& _out)
{
    GRBVar vars[2] = { _in1,_in2 };
    model.addGenConstrOr(_out, vars, 2);
}

void newXor(GRBModel& model, vector<GRBVar> & _ins, GRBVar &_out)
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

static set<int> QuadExtractIndex(int i, int j,bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
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



GRBVar newMaj(GRBModel& model, vector<GRBVar>& s, int i, int j, int k, bitset<293>& nonZeroTrackRep,bitset<293> & varsTrackRep)
{
    set<set<int>> rSets;
    set<int> tmpSet = QuadExtractIndex(i, j, nonZeroTrackRep, varsTrackRep);
    if(tmpSet.size() > 0) rSets.emplace(tmpSet);
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
    //cout << __func__ << " : " << size0 << "   " << size1 << endl;

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

    if (IsNonZero2(k, j))
        AndOpr2(k, j);*/

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        newXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

GRBVar newCh(GRBModel& model, vector<GRBVar>& s, int i, int j, int k, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
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
    // cout << __func__ << " : " << size0 << "   " << size1 << endl;

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
        newXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

void newXorF(GRBModel& model, vector<GRBVar>& s, int i, int j, int k, bitset<293>& nonZeroTrackRep,bitset<293> & varsTrackRep)
{
    vector<GRBVar> vars;
    
    vars.emplace_back(s[i]);
    if (IsNonZero1(k))
        AndOpr1(k);
    if (IsNonZero1(j))
        AndOpr1(j);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        newXor(model, vars, O);
    else
        model.addConstr(O == 0);

    s[i] = O;

    xorFTrack(nonZeroTrackRep, i, j, k);
    xorFVarsTrack(varsTrackRep, i, j, k);
}

GRBVar newKsg(GRBModel& model, vector<GRBVar>& s, bitset<293>& nonZeroTrackRep,bitset<293> & varsTrackRep)
{
    vector<GRBVar> vars;
    if (IsNonZero1(12))
        AndOpr1(12);

    if (IsNonZero1(154))
        AndOpr1(154);

    GRBVar c = newMaj(model, s, 235, 61, 193, nonZeroTrackRep,varsTrackRep);
    GRBVar d = newCh(model, s, 230, 111, 66, nonZeroTrackRep,varsTrackRep);
    vars.emplace_back(c);
    vars.emplace_back(d);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        newXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

GRBVar newFbk(GRBModel& model, vector<GRBVar>& s, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    vector<GRBVar> vars;
    if (IsNonZero1(0))
        AndOpr1(0);

    if (IsNonZero1(107))
        AndOpr1(107);

    if (IsNonZero1(196))
        AndOpr1(196);

    GRBVar c = newMaj(model, s, 244, 23, 160, nonZeroTrackRep,varsTrackRep);
    vars.emplace_back(c);

    GRBVar O = model.addVar(0, 1, 0, GRB_BINARY);
    if (vars.size() > 0)
        newXor(model, vars, O);
    else
        model.addConstr(O == 0);

    return O;
}

GRBVar newGenM(GRBModel& model, bitset<128>& cube, vector<GRBVar>& v, int r)
{
    if (r < 128 || r >= 256)
    {
        logger(__func__ + string(" : parameter rounds is out of range."));
        exit(-1);
    }

    if (cube[r - 128] == 0)
    {
        logger(__func__ + string(" : error func call."));
        exit(-2);
    }


    GRBVar a = newtap(model, v[r - 128]);
    return a;
}

// starting from round 0
void newupdate0(GRBModel& model, bitset<128>& cube, vector<GRBVar>& s, vector<GRBVar>& v, int r, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    newXorF(model, s, 289, 235, 230, nonZeroTrackRep,varsTrackRep);
    newXorF(model, s, 230, 196, 193, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 193, 160, 154, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 154, 111, 107, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 107, 66, 61, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 61, 23, 0, nonZeroTrackRep, varsTrackRep);

    vector<GRBVar> vars;
    GRBVar e = newKsg(model, s, nonZeroTrackRep, varsTrackRep);
    GRBVar b = newFbk(model,s, nonZeroTrackRep, varsTrackRep);
    vars.emplace_back(e);
    vars.emplace_back(b);

    if (r >= 128 && r < 256 && (cube[r - 128] == 1))
    {
        GRBVar m = newGenM(model, cube, v, r);
        vars.emplace_back(m);
    }

    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);
    newXor(model, vars, o);

    model.addConstr(s[0] == 0);
    for (int i = 0; i < 292; i++)
        s[i] = s[i + 1];
    s[292] = o;
}

// starting from 256
void newupdate256(GRBModel& model, vector<GRBVar>& s, bitset<293>& nonZeroTrackRep, bitset<293>& varsTrackRep)
{
    newXorF(model, s, 289, 235, 230, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 230, 196, 193, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 193, 160, 154, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 154, 111, 107, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 107, 66, 61, nonZeroTrackRep, varsTrackRep);
    newXorF(model, s, 61, 23, 0, nonZeroTrackRep, varsTrackRep);

    GRBVar e = newKsg(model, s, nonZeroTrackRep, varsTrackRep);
    GRBVar b = newFbk(model, s, nonZeroTrackRep, varsTrackRep);

    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(o == e + b);

    model.addConstr(s[0] == 0);
    for (int i = 0; i < 292; i++)
        s[i] = s[i + 1];
    s[292] = o;
}

void setRoundFlags(int r, bitset<293>& nonZeroFlag, bitset<293>& varsFlag, vector<bitset<293>>& nonZeroTrackReps, vector<bitset<293>>& varsTrackReps)
{
    if (r < nonZeroTrackReps.size())
        nonZeroFlag = nonZeroTrackReps[r]; 
    else
        nonZeroFlag.set();

    if (r < varsTrackReps.size())
        varsFlag = varsTrackReps[r];
    else
        varsFlag.set();
}


STATUS newAcornForwardExpand(int rounds, bitset<128>& cube, vector<bitset<293>>& terms, float time, int thread)
{
    vector<bitset<293>> nonZeroTrackReps;
    vector<bitset<293>> varsTrackReps;

    acornNonZeroTrack(nonZeroTrackReps);
    acornVarsTrack0(varsTrackReps, cube,nonZeroTrackReps);

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_MIPFocus, 1);
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    // input constraints
    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for (int i = 0; i < 293; i++)
        model.addConstr(is[i] == 0);

    vector<GRBVar> iV(128);
    for (int i = 0; i < 128; i++)
        iV[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for (int i = 0; i < 128; i++)
        if (cube[i])
            model.addConstr(iV[i] == 1);
        else
            model.addConstr(iV[i] == 0);

    vector<GRBVar> s = is;
    vector<GRBVar> V = iV;

    for (int r = 0; r < rounds; r++)
    {
        bitset<293> nonZeroFlags, varsFlags;
        setRoundFlags(r, nonZeroFlags, varsFlags, nonZeroTrackReps, varsTrackReps);
        newupdate0(model, cube, s, V, r, nonZeroFlags, varsFlags);
    }

    for (int i = 0; i < (rounds - 128 > 128 ? 128 : rounds - 128); i++)
        model.addConstr(V[i] == 0);

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
            map<bitset<293>, int, CMPS<293>> counterMap;
            bitset<293> last;
            for (int i = 0; i < solCount; i++)
            {
                model.set(GRB_IntParam_SolutionNumber, i);
                for (int j = 0; j < 293; j++)
                    if (round(s[j].get(GRB_DoubleAttr_Xn)) == 1)
                        last[j] = 1;
                    else
                        last[j] = 0;
                counterMap[last]++;
            }

            for (auto& x : counterMap)
                terms.emplace_back(x.first);
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }
}

STATUS newMidSolutionCounter(int rounds, bitset<293> & startcube, bitset<293+256> & last, map<bitset<293+256>,
    int, CMPS<293+256>> &counterMap, float time, int threads, int midround)
{
    if (midround > rounds || rounds < 0 || midround < 0 || rounds < 256 || midround <256)
    {
        logger(__func__ + string(":parameter error."));
        exit(-1);
    }

    counterMap.clear();

    vector<bitset<293>> nonZeroTrackReps;
    vector<bitset<293>> varsTrackReps;

    acornNonZeroTrack(nonZeroTrackReps);
    acornVarsTrack256(varsTrackReps, startcube,nonZeroTrackReps);

    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, MAX);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_Presolve, 2);

    // Create the model
    GRBModel model = GRBModel(env);

    // input constraints
    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for (int i = 0; i < 293; i++)
        if (startcube[i])
            model.addConstr(is[i] == 1);
        else
            model.addConstr(is[i] == 0);

    vector<GRBVar> s = is;
    // update
    for (int r = 256; r < midround; r++)
    {
        bitset<293> nonZeroFlag;
        bitset<293> varsFlag;
        setRoundFlags(r, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);
        newupdate256(model, s, nonZeroFlag, varsFlag);
    }

    // mid constraints
    bitset<293> midNonZeroFlag;
    bitset<293> midVarsFlag;
    setRoundFlags(midround, midNonZeroFlag, midVarsFlag, nonZeroTrackReps, varsTrackReps);

    vector<GRBVar> mids2 = s;


    vector<GRBVar> mids(293);
    for (int i = 0; i < 293; i++)
        mids[i] = model.addVar(0, 1, 0, GRB_BINARY);

    GRBLinExpr nk = 0;
    for (int i = 0; i < 293; i++)
        if (midVarsFlag[i])
            model.addConstr(mids[i] == s[i]);
        else
            if (midNonZeroFlag[i])
            {
                model.addConstr(mids[i] >= s[i]);
                nk += mids[i] - s[i];
            }
            else
                model.addConstr(mids[i] == 0);

    s = mids;
    vector<GRBVar> midK(128);
    for (int i = 0; i < 128; i++)
    {
        midK[i] = model.addVar(0, 1, 0, GRB_BINARY);
        nk += midK[i];
    }
    vector<GRBVar> midV(128);

    vector<GRBVar> K = midK;
    vector<GRBVar> V = midV;

    // second Phase
    for (int r = midround; r < rounds; r++)
        update(model, s, K, V, r);

    // output constraints
    // key + iv + 288
    for (int i = 0; i < 293; i++)
        if (last[i+256])
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);
    for (int i = 0; i < 128; i++)
        if (last[i])
            model.addConstr(K[i] == 1);
        else
            model.addConstr(K[i] == 0);


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
            bitset<293+256> midState;
            bitset<293> midState2;
            for (int i = 0; i < solCount; i++)
            {
                model.set(GRB_IntParam_SolutionNumber, i);

                // midstate
                for (int j = 0; j < 293; j++)
                    if (round(mids[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState[j+256] = 1;
                    else
                        midState[j+256] = 0;
                for (int j = 0; j < 128; j++)
                    if (round(midK[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState[j] = 1;
                    else
                        midState[j] = 0;

                counterMap[midState]++;
            }
            return SOLUTION;
        }
        else
        {
            return NOSOLUTION;
        }
    }


}

STATUS newMidSolutionCounter(int rounds, bitset<293>& startcube, bitset<293>& last, int & solcnt, float time, int threads)
{
    vector<bitset<293>> nonZeroTrackReps;
    vector<bitset<293>> varsTrackReps;

    acornNonZeroTrack(nonZeroTrackReps);
    acornVarsTrack256(varsTrackReps, startcube,nonZeroTrackReps);

    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, MAX);
    env.set(GRB_IntParam_MIPFocus, 1);
    env.set(GRB_IntParam_Threads, threads);

    // Create the model
    GRBModel model = GRBModel(env);

    // input constraints
    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for (int i = 0; i < 293; i++)
        if (startcube[i])
            model.addConstr(is[i] == 1);
        else
            model.addConstr(is[i] == 0);

    vector<GRBVar> s = is;
    // update
    for (int r = 256; r < rounds; r++)
    {
        bitset<293> nonZeroFlag;
        bitset<293> varsFlag;
        setRoundFlags(r, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);
        newupdate256(model, s, nonZeroFlag, varsFlag);
    }

    bitset<293> lastNonZeroFlag, lastVarsFlag;
    setRoundFlags(rounds, lastNonZeroFlag, lastVarsFlag, nonZeroTrackReps, varsTrackReps);

    // output constraints
    for (int i = 0; i < 293; i++)
        if (last[i] && lastVarsFlag[i])
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);

    // objective function

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

        solcnt = solCount;

        if (solCount > 0)
        {
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }


}