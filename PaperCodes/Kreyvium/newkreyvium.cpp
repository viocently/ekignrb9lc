#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<algorithm>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"kreyviumtrack.h"
#include"kreyvium.h"

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

GRBVar newLFSR(GRBModel& model, vector <GRBVar>& x)
{
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);

    // copy
    model.addConstr(x[0] <= a + b);
    model.addConstr(x[0] >= a);
    model.addConstr(x[0] >= b);

    vector<GRBVar> t = x;
    for (int i = 0; i < 127; i++)
        x[i] = t[i + 1];
    x[127] = a;
    return b;
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

void newKreyviumCore(GRBModel& model, vector<GRBVar>& s, int i1, int i2, int i3, int i4, int i5,bitset<288+256> & nonZeroTrackRep,bitset<288+256> & varsTrackRep)
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

void newupdate(GRBModel& model, vector<GRBVar>& s, vector<GRBVar> & key,vector<GRBVar> & iv, bitset<288 + 256>& nonZeroTrackRep, bitset<288 + 256>& varsTrackRep)
{
    // k dont need to LFSR because its all 0 ,so the k in the register has no contribution to the propagation
    GRBVar genIV;
    if (nonZeroTrackRep[288] == 1)
        genIV = newLFSR(model, iv);
    else
    {
        genIV = model.addVar(0, 1, 0, GRB_BINARY);
        model.addConstr(genIV == 0);
        vector<GRBVar> tmpiv = iv;
        for(int i =0;i<128;i++)
            iv[i] = tmpiv[(i + 1) % 128];
    }

    vector<GRBVar> tmpkey = key;
    for (int i = 0; i < 128; i++)
        key[i] = tmpkey[(i + 1)%128];

    newKreyviumCore(model, s, 65, 170, 90, 91, 92, nonZeroTrackRep,varsTrackRep);
    newKreyviumCore(model, s, 161, 263, 174, 175, 176, nonZeroTrackRep, varsTrackRep);
    newKreyviumCore(model, s, 242, 68, 285, 286, 287, nonZeroTrackRep, varsTrackRep);

    GRBVar t1 = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(t1 == s[92] + genIV);

    s[92] = t1;
    vector<GRBVar> tmp = s;
    for (int i = 0; i < 288; i++)
        s[(i + 1) % 288] = tmp[i];
}

void setRoundFlags(int r, bitset<128>& cube, bitset<288 + 256>& nonZeroFlag, bitset<288 + 256>& varsFlag, vector<bitset<288 + 256>>& nonZeroTrackReps, vector<bitset<288 + 256>>& varsTrackReps)
{
    nonZeroFlag.reset();
    varsFlag.reset();

    if (r < nonZeroTrackReps.size())
        nonZeroFlag = nonZeroTrackReps[r];
    else
    {
        for (int i = 0; i < 288; i++)
            nonZeroFlag[i] = 1;
        for (int i = 0; i < 128; i++)
            nonZeroFlag[i + 288 + 128] = 1;

        for (int i = 0; i < 128; i++)
            if (cube[ 127 - ((i +r)%128)])
                nonZeroFlag[i+288] = 1;
    }

    if (r < varsTrackReps.size())
        varsFlag = varsTrackReps[r];
    else
    {
        for (int i = 0; i < 288; i++)
            varsFlag[i] = 1;

        for (int i = 0; i < 128; i++)
            if (cube[127 - ((i + r) % 128)])
                varsFlag[i+288] = 1;
    }
}

STATUS newMidSolutionCounter(int rounds, bitset<128>& cube, bitset<544>& last, map<bitset<288 + 256>,
    vector<bitset<288+256>>, CMPS<288 + 256>> &counterMap, float time, int threads, int midround)
{
    if (midround > rounds || rounds < 0 || midround < 0)
    {
        logger(__func__ + string(":parameter error."));
        exit(-1);
    }

    counterMap.clear();


    vector<bitset<544>> nonZeroTrackReps;
    kreyviumTrack(cube, nonZeroTrackReps);


    vector<bitset<544>> varsTrackReps;
    kreyviumVarsTrack(cube, varsTrackReps,nonZeroTrackReps);


    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, MAX);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_MIPFocus, 3);

    // Create the model
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> KEY(128);
    for (int i = 0; i < 128; i++)
        KEY[i] = model.addVar(0, 1, 0, GRB_BINARY);
    vector<GRBVar> IV(128);
    for (int i = 0; i < 128; i++)
        IV[i] = model.addVar(0, 1, 0, GRB_BINARY);

    // set cube
    for (int i = 0; i < 128; i++)
        if (cube[i] == 1)
            model.addConstr(IV[i] == 1);
        else
            model.addConstr(IV[i] == 0);

    // set key
    for (int i = 0; i < 128; i++)
        model.addConstr(KEY[i] == 0);

    // register
    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
    vector<GRBVar> key(128);
    for (int i = 0; i < 128; i++)
        key[i] = model.addVar(0, 1, 0, GRB_BINARY);
    vector<GRBVar> iv(128);
    for (int i = 0; i < 128; i++)
        iv[i] = model.addVar(0, 1, 0, GRB_BINARY);

    // distribute key
    for (int i = 0; i < 93; i++)
        model.addConstr(s[i] == 0);
    for (int i = 0; i < 128; i++)
        model.addConstr(key[127 - i] == 0);
    // distribute iv
    for (int i = 0; i < 128; i++)
    {
        model.addConstr(IV[i] <= iv[127 - i] + s[93 + i]);
        model.addConstr(IV[i] >= iv[127 - i]);
        model.addConstr(IV[i] >= s[93 + i]);
    }

    // constant 
    for (int i = 93 + 128; i < 288; i++)
        model.addConstr(s[i] == 0);

    // update
    for (int r = 0; r < midround; r++)
    {
        bitset<288 + 256> nonZeroFlag;
        bitset<288 + 256> varsFlag;
        setRoundFlags(r, cube, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);
        newupdate(model, s, key,iv, nonZeroFlag, varsFlag);
    }

    // mid constraints
    bitset<288 + 256> midNonZeroFlag;
    bitset<288 + 256> midVarsFlag;
    setRoundFlags(midround, cube, midNonZeroFlag, midVarsFlag, nonZeroTrackReps, varsTrackReps);


    vector<GRBVar> mids2 = s;
    vector<GRBVar> midkey2 = key;
    vector<GRBVar> midiv2 = iv;


    vector<GRBVar> mids(288);
    for (int i = 0; i < 288; i++)
        mids[i] = model.addVar(0, 1, 0, GRB_BINARY);
    vector<GRBVar> midkey(128);
    for (int i = 0; i < 128; i++)
        midkey[i] = model.addVar(0, 1, 0, GRB_BINARY);
    vector<GRBVar> midiv(128);
    for (int i = 0; i < 128; i++)
        midiv[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for (int i = 0; i < 288; i++)
        if (midVarsFlag[i])
            model.addConstr(mids[i] == s[i]);
        else
            if (midNonZeroFlag[i])
                model.addConstr(mids[i] >= s[i]);
            else
                model.addConstr(mids[i] == 0);

    for (int i = 0; i < 128; i++)
        if (midVarsFlag[i + 288])
            model.addConstr(midiv[i] == iv[i]);
        else
            model.addConstr(midiv[i] == 0);

    for (int i = 0; i < 128; i++)
        model.addConstr(midkey[i] >= key[i]);

    key = midkey;
    s = mids;
    iv = midiv;

    // second Phase
    for (int r = midround; r < rounds; r++)
        update(model, s, key, iv);

    // output constraints
    // key + iv + 288
    for (int i = 0; i < 128; i++)
        if (last[i])
            model.addConstr(key[i] == 1);
        else
            model.addConstr(key[i] == 0);
    for (int i = 0; i < 128; i++)
        if (last[i + 128])
            model.addConstr(iv[i] == 1);
        else
            model.addConstr(iv[i] == 0);
    for (int i = 0; i < 288; i++)
        if(last[i+256])
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);

    // objective function
    GRBLinExpr nk = 0;
    for (int i = 0; i < 128; i++)
        nk += midkey[i];

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
            bitset<256+288> midState;
            bitset<256+288> midState2;
            for (int i = 0; i < solCount; i++)
            {
                model.set(GRB_IntParam_SolutionNumber, i);

                // midstate
                for (int j = 0; j < 128; j++)
                {
                    if (round(midkey[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState[j] = 1;
                    else
                        midState[j] = 0;
                    if (round(midiv[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState[j + 128] = 1;
                    else
                        midState[j + 128] = 0;
                }
                for(int j=0;j<288;j++)
                    if (round(mids[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState[j + 256] = 1;
                    else
                        midState[j + 256] = 0;

                // midstate2
                for (int j = 0; j < 128; j++)
                {
                    if (round(midkey2[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState2[j] = 1;
                    else
                        midState2[j] = 0;
                    if (round(midiv2[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState2[j + 128] = 1;
                    else
                        midState2[j + 128] = 0;
                }
                for (int j = 0; j < 288; j++)
                    if (round(mids2[j].get(GRB_DoubleAttr_Xn)) == 1)
                        midState2[j + 256] = 1;
                    else
                        midState2[j + 256] = 0;


              

                counterMap[midState].emplace_back(midState2);
            }
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }


}

STATUS newMidSolutionCounter(int rounds, bitset<128>& cube, bitset<544>& last,
    int & sols, float time, int threads)
{


    vector<bitset<544>> nonZeroTrackReps;
    kreyviumTrack(cube, nonZeroTrackReps);


    vector<bitset<544>> varsTrackReps;
    kreyviumVarsTrack(cube, varsTrackReps, nonZeroTrackReps);


    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, MAX);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_MIPFocus, 1);

    // Create the model
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> KEY(128);
    for (int i = 0; i < 128; i++)
        KEY[i] = model.addVar(0, 1, 0, GRB_BINARY);
    vector<GRBVar> IV(128);
    for (int i = 0; i < 128; i++)
        IV[i] = model.addVar(0, 1, 0, GRB_BINARY);

    // set cube
    for (int i = 0; i < 128; i++)
        if (cube[i] == 1)
            model.addConstr(IV[i] == 1);
        else
            model.addConstr(IV[i] == 0);

    // set key
    for (int i = 0; i < 128; i++)
        model.addConstr(KEY[i] == 0);

    // register
    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);
    vector<GRBVar> key(128);
    for (int i = 0; i < 128; i++)
        key[i] = model.addVar(0, 1, 0, GRB_BINARY);
    vector<GRBVar> iv(128);
    for (int i = 0; i < 128; i++)
        iv[i] = model.addVar(0, 1, 0, GRB_BINARY);

    // distribute key
    for (int i = 0; i < 128; i++)
    {
        if (i < 93)
        {
            model.addConstr(KEY[i] <= s[i] + key[127 - i]);
            model.addConstr(KEY[i] >= s[i]);
            model.addConstr(KEY[i] >= key[127 - i]);
        }
        else
        {
            model.addConstr(KEY[i] == key[127 - i]);
        }
    }
    // distribute iv
    for (int i = 0; i < 128; i++)
    {
        model.addConstr(IV[i] <= iv[127 - i] + s[93 + i]);
        model.addConstr(IV[i] >= iv[127 - i]);
        model.addConstr(IV[i] >= s[93 + i]);
    }

    // constant 
    for (int i = 93 + 128; i < 288; i++)
        model.addConstr(s[i] == 0);

    // update
    for (int r = 0; r < rounds; r++)
    {
        bitset<288 + 256> nonZeroFlag;
        bitset<288 + 256> varsFlag;
        setRoundFlags(r, cube, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);
        newupdate(model, s, key, iv, nonZeroFlag, varsFlag);
    }

    // output constraints
    // key + iv + 288
    for (int i = 0; i < 128; i++)
        if (last[i])
            model.addConstr(key[i] == 1);
        else
            model.addConstr(key[i] == 0);
    for (int i = 0; i < 128; i++)
        if (last[i + 128])
            model.addConstr(iv[i] == 1);
        else
            model.addConstr(iv[i] == 0);
    for (int i = 0; i < 288; i++)
        if (last[i + 256])
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
