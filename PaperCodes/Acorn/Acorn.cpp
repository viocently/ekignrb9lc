#include<iostream>
#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"

const int MAX = 2000000000;

GRBVar maj(GRBModel& model, vector<GRBVar>& x, int i, int j, int k)
{
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar c = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x0 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x2 = model.addVar(0, 1, 0, GRB_BINARY);

    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);

    // copy
    model.addConstr(x[i] <= a + b + x0);
    model.addConstr(x[i] >= a);
    model.addConstr(x[i] >= b);
    model.addConstr(x[i] >= x0);
    x[i] = x0;

    model.addConstr(x[j] <= a + c + x1);
    model.addConstr(x[j] >= a);
    model.addConstr(x[j] >= c);
    model.addConstr(x[j] >= x1);
    x[j] = x1;

    model.addConstr(x[k] <= b + c + x2);
    model.addConstr(x[k] >= b);
    model.addConstr(x[k] >= c);
    model.addConstr(x[k] >= x2);
    x[k] = x2;

    model.addConstr(o == a + b + c);

    return o;
}

GRBVar ch(GRBModel& model, vector<GRBVar>& x, int i, int j, int k)
{
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar c = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x0 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x2 = model.addVar(0, 1, 0, GRB_BINARY);

    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);

    // copy
    model.addConstr(x[i] <= a + b + x0);
    model.addConstr(x[i] >= a);
    model.addConstr(x[i] >= b);
    model.addConstr(x[i] >= x0);
    x[i] = x0;

    model.addConstr(x[j] <= a + x1);
    model.addConstr(x[j] >= a);
    model.addConstr(x[j] >= x1);
    x[j] = x1;

    model.addConstr(x[k] <= b + c + x2);
    model.addConstr(x[k] >= b);
    model.addConstr(x[k] >= c);
    model.addConstr(x[k] >= x2);
    x[k] = x2;

    model.addConstr(o == a + b + c);

    return o;
}

void xorF(GRBModel& model, vector<GRBVar>& x, int i, int j, int k)
{
    GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar c = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x2 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);

    // copy
    model.addConstr(x[j] <= b + x1);
    model.addConstr(x[j] >= b);
    model.addConstr(x[j] >= x1);
    x[j] = x1;

    // copy
    model.addConstr(x[k] <= c + x2);
    model.addConstr(x[k] >= c);
    model.addConstr(x[k] >= x2);
    x[k] = x2;

    model.addConstr(o == x[i] + b + c);

    x[i] = o;
}

GRBVar ksg(GRBModel& model, vector<GRBVar>& x)
{
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x0 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x1 = model.addVar(0, 1, 0, GRB_BINARY);

    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);

    GRBVar c = maj(model, x, 235, 61, 193);
    GRBVar d = ch(model, x, 230, 111, 66);

    // copy
    model.addConstr(x[12] <= a + x0);
    model.addConstr(x[12] >= a);
    model.addConstr(x[12] >= x0);
    x[12] = x0;

    model.addConstr(x[154] <= b + x1);
    model.addConstr(x[154] >= b);
    model.addConstr(x[154] >= x1);
    x[154] = x1;

    model.addConstr(o == a + b + c + d);

    return o;
}

GRBVar fbk(GRBModel& model, vector<GRBVar>& x)
{
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x0 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar d = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x2 = model.addVar(0, 1, 0, GRB_BINARY);

    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);

    // ksg

    // maj
    GRBVar c = maj(model, x, 244, 23, 160);


    // copy
    model.addConstr(x[0] <= a + x0);
    model.addConstr(x[0] >= a);
    model.addConstr(x[0] >= x0);
    x[0] = x0;

    model.addConstr(x[107] <= b + x1);
    model.addConstr(x[107] >= b);
    model.addConstr(x[107] >= x1);
    x[107] = x1;

    model.addConstr(x[196] <= d + x2);
    model.addConstr(x[196] >= d);
    model.addConstr(x[196] >= x2);
    x[196] = x2;

    // 
    model.addConstr(o >= a + b + c + d);

    return o;
}

GRBVar genM(GRBModel& model, vector<GRBVar>& k, vector<GRBVar>& v,
    int r)
{
    GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar x = model.addVar(0, 1, 0, GRB_BINARY);

    if (r < 128)
    {
        model.addConstr(k[r] <= a + x);
        model.addConstr(k[r] >= a);
        model.addConstr(k[r] >= x);
        k[r] = x;
        return a;
    }
    else if (r < 256)
    {
        model.addConstr(v[r - 128] <= a + x);
        model.addConstr(v[r - 128] >= a);
        model.addConstr(v[r - 128] >= x);
        v[r - 128] = x;
        return a;
    }
    else
    {
        model.addConstr(k[r % 128] <= a + x);
        model.addConstr(k[r % 128] >= a);
        model.addConstr(k[r % 128] >= x);
        k[r % 128] = x;
        if (r != 256)
            return a;
        else
        {
            GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);
            model.addConstr(o >= a);
            return o;
        }
    }
}

void update(GRBModel& model, vector<GRBVar>& x, vector<GRBVar>& k,
    vector<GRBVar>& v, int r)
{
    xorF(model, x, 289, 235, 230);
    xorF(model, x, 230, 196, 193);
    xorF(model, x, 193, 160, 154);
    xorF(model, x, 154, 111, 107);
    xorF(model, x, 107, 66, 61);
    xorF(model, x, 61, 23, 0);

    GRBVar e = ksg(model, x);
    auto b = fbk(model, x);
    auto m = genM(model, k, v, r);

    GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);
    model.addConstr(o == e + m + b);

    model.addConstr(x[0] == 0);
    for (int i = 0; i < 292; i++)
        x[i] = x[i + 1];
    x[292] = o;
}

int BackExpandPolynomial(int current, int rounds, vector<bitset<549>>& term)
{
    if (term.size() > 0)
    {
        cerr << __func__ << " Term is not empty" << endl;
        exit(-1);
    }
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, MAX);

    GRBModel model = GRBModel(env);

    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> iK(128);
    for (int i = 0; i < 128; i++)
        iK[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> iV(128);
    for (int i = 0; i < 128; i++)
        iV[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> s = is;
    vector<GRBVar> K = iK;
    vector<GRBVar> V = iV;

    for (int r = current - rounds; r < current; r++)
        update(model, s, K, V, r);

    xorF(model, s, 289, 235, 230);
    xorF(model, s, 230, 196, 193);
    xorF(model, s, 193, 160, 154);
    xorF(model, s, 154, 111, 107);
    xorF(model, s, 107, 66, 61);
    xorF(model, s, 61, 23, 0);

    auto o = ksg(model, s);

    model.addConstr(o == 1);

    for (int i = 0; i < 128; i++)
        model.addConstr(K[i] == 0);

    for (int i = 0; i < 128; i++)
        model.addConstr(V[i] == 0);

    for (int i = 0; i < 293; i++)
        model.addConstr(s[i] == 0);

    model.optimize();

    map<bitset<549>, int, CMPS<549>> counterMap;

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
        if (solCount >= MAX)
        {
            cerr << "solCount value  is too big !" << endl;
            exit(-1);
        }
        bitset<549> start;
        for (int i = 0; i < solCount; i++)
        {
            model.set(GRB_IntParam_SolutionNumber, i);

            for (int j = 0; j < 128; j++)
                if (round(iK[j].get(GRB_DoubleAttr_Xn)) == 1)
                    start[j] = 1;
                else
                    start[j] = 0;

            for (int j = 0; j < 128; j++)
                if (round(iV[j].get(GRB_DoubleAttr_Xn)) == 1)
                    start[j + 128] = 1;
                else
                    start[j + 128] = 0;

            for (int j = 0; j < 293; j++)
                if (round(is[j].get(GRB_DoubleAttr_Xn)) == 1)
                    start[j + 256] = 1;
                else
                    start[j + 256] = 0;

            counterMap[start]++;
        }
    }
    else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        exit(-2);
    }
    else
    {
        exit(-1);
    }

    for (auto& it : counterMap)
        if (it.second % 2 == 1)
            term.push_back(it.first);

    logger(string(__func__) + string(":") + string("Expand Terms: ") + to_string(term.size()));

    return 0;
}




int SecondBackExpandPolynomial(int current, int rounds, const bitset<549>& last,
    vector<bitset<549>>& term, int threads)
{
    if (term.size() > 0)
    {
        cerr << __func__ << " Term is not empty" << endl;
        exit(-1);
    }
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_PoolSolutions, MAX);

    GRBModel model = GRBModel(env);

    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> iK(128);
    for (int i = 0; i < 128; i++)
        iK[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> iV(128);
    for (int i = 0; i < 128; i++)
        iV[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> s = is;
    vector<GRBVar> K = iK;
    vector<GRBVar> V = iV;

    for (int r = current - rounds; r < current; r++)
        update(model, s, K, V, r);

    for (int i = 0; i < 128; i++)
        if (last[i] == 0)
            model.addConstr(K[i] == 0);
        else
            model.addConstr(K[i] == 1);

    for (int i = 0; i < 128; i++)
        if (last[i + 128] == 0)
            model.addConstr(V[i] == 0);
        else
            model.addConstr(V[i] == 1);

    for (int i = 0; i < 293; i++)
        if (last[i + 256] == 0)
            model.addConstr(s[i] == 0);
        else
            model.addConstr(s[i] == 1);

    model.optimize();

    map<bitset<549>, int, CMPS<549>> counterMap;

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        int solCount = model.get(GRB_IntAttr_SolCount);
        if (solCount >= MAX)
        {
            cerr << "solCount value  is too big !" << endl;
            exit(-1);
        }
        bitset<549> start;
        for (int i = 0; i < solCount; i++)
        {
            model.set(GRB_IntParam_SolutionNumber, i);

            for (int j = 0; j < 128; j++)
                if (round(iK[j].get(GRB_DoubleAttr_Xn)) == 1)
                    start[j] = 1;
                else
                    start[j] = 0;

            for (int j = 0; j < 128; j++)
                if (round(iV[j].get(GRB_DoubleAttr_Xn)) == 1)
                    start[j + 128] = 1;
                else
                    start[j + 128] = 0;

            for (int j = 0; j < 293; j++)
                if (round(is[j].get(GRB_DoubleAttr_Xn)) == 1)
                    start[j + 256] = 1;
                else
                    start[j + 256] = 0;

            counterMap[start]++;
        }
    }
    else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        exit(-2);
    }
    else
    {
        exit(-1);
    }

    for (auto& it : counterMap)
        if (it.second % 2 == 1)
            term.push_back(it.first);

    logger(string(__func__) + string(":") + string("Expand Terms: ") + to_string(term.size()));
    return 0;
}

STATUS MidSolutionCounter(int midround, int rounds, const bitset<549>& start, const
    bitset<549>& last, int& solcnt, float
    time, int thread)
{
    if (midround < 256)
    {
        logger(__func__ + string(" : midround should be larger than 256."));
        exit(-1);
    }

    if (midround > rounds)
    {
        logger(__func__ + string(" : midround should not be larger than rounds."));
        exit(-2);
    }


    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_MIPFocus, 1);
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);
    for (int i = 0; i < 293; i++)
        if (start[i + 256])
            model.addConstr(is[i] == 1);
        else
            model.addConstr(is[i] == 0);


    vector<GRBVar> iK(128);
    for (int i = 0; i < 128; i++)
        iK[i] = model.addVar(0, 1, 0, GRB_BINARY);
    for (int i = 0; i < 128; i++)
        if (start[i])
            model.addConstr(iK[i] == 1);
        else
            model.addConstr(iK[i] == 0);


    vector<GRBVar> iV(128);

    vector<GRBVar> s = is;
    vector<GRBVar> K = iK;
    vector<GRBVar> V = iV;

    for (int r = midround; r < rounds; r++)
        update(model, s, K, V, r);

    for (int i = 0; i < 293; i++)
        if (last[256+i])
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);

    for (int i = 0; i < 128; i++)
        if (last[i])
            model.addConstr(K[i] == 1);
        else
            model.addConstr(K[i] == 0);

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
        if (solcnt == 0)
            return NOSOLUTION;
        else
            return SOLUTION;
    }
}

STATUS MidSolutionCounter(int rounds, const bitset<128>& cube, const
    bitset<549>& last, map<bitset<128>, int, CMPS<128>>& counterMap, float
    time, int thread)
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for (int i = 0; i < 293; i++)
        model.addConstr(is[i] == 0);

    vector<GRBVar> iK(128);
    for (int i = 0; i < 128; i++)
        iK[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> iV(128);
    for (int i = 0; i < 128; i++)
        iV[i] = model.addVar(0, 1, 0, GRB_BINARY);

    // cube
    for (int i = 0; i < 128; i++)
        if (cube[i] == 0)
            model.addConstr(iV[i] == 0);
        else
            model.addConstr(iV[i] == 1);

    vector<GRBVar> s = is;
    vector<GRBVar> K = iK;
    vector<GRBVar> V = iV;

    for (int r = 0; r < rounds; r++)
        update(model, s, K, V, r);

    for (int i = 0; i < 128; i++)
        if (last[i] == 0)
            model.addConstr(K[i] == 0);
        else
            model.addConstr(K[i] == 1);

    for (int i = 0; i < 128; i++)
        if (last[i + 128] == 0)
            model.addConstr(V[i] == 0);
        else
            model.addConstr(V[i] == 1);

    for (int i = 0; i < 293; i++)
        if (last[i + 256] == 0)
            model.addConstr(s[i] == 0);
        else
            model.addConstr(s[i] == 1);

    if (time > 0)
        model.set(GRB_DoubleParam_TimeLimit, time);

    GRBLinExpr nk = 0;
    for (int i = 0; i < 128; i++)
        nk += iK[i];
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
        /*logger(__func__ + string(" : ") + to_string(rounds) + string(" | "
        ) + to_string(time) + string(" | "
        ) + to_string(solCount));*/


        if (solCount > 0)
        {
            bitset<128> start;
            for (int i = 0; i < solCount; i++)
            {
                model.set(GRB_IntParam_SolutionNumber, i);
                for (int j = 0; j < 128; j++)
                    if (round(iK[j].get(GRB_DoubleAttr_Xn)) == 1)
                        start[j] = 1;
                    else
                        start[j] = 0;
                counterMap[start]++;
            }
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }
}

// this MidSolutionCounter is to compute the expression of y^v in K with bits in cube set 1
STATUS MidSolutionCounter2(int rounds, const bitset<128>& cube, const
    bitset<549>& last, map<bitset<128>, int, CMPS<128>>& counterMap, float
    time, int thread)
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);

    for (int i = 0; i < 293; i++)
        model.addConstr(is[i] == 0);

    vector<GRBVar> iK(128);
    for (int i = 0; i < 128; i++)
        iK[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> iV(128);
    for (int i = 0; i < 128; i++)
        iV[i] = model.addVar(0, 1, 0, GRB_BINARY);

    // cube
    for (int i = 0; i < 128; i++)
        if (cube[i] == 0)
            model.addConstr(iV[i] == 0);

    vector<GRBVar> s = is;
    vector<GRBVar> K = iK;
    vector<GRBVar> V = iV;

    for (int r = 0; r < rounds; r++)
        update(model, s, K, V, r);

    for (int i = 0; i < 128; i++)
        if (last[i] == 0)
            model.addConstr(K[i] == 0);
        else
            model.addConstr(K[i] == 1);

    for (int i = 0; i < 128; i++)
        if (last[i + 128] == 0)
            model.addConstr(V[i] == 0);
        else
            model.addConstr(V[i] == 1);

    for (int i = 0; i < 293; i++)
        if (last[i + 256] == 0)
            model.addConstr(s[i] == 0);
        else
            model.addConstr(s[i] == 1);

    if (time > 0)
        model.set(GRB_DoubleParam_TimeLimit, time);

    GRBLinExpr nk = 0;
    for (int i = 0; i < 128; i++)
        nk += iK[i];
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
        /*logger(__func__ + string(" : ") + to_string(rounds) + string(" | "
        ) + to_string(time) + string(" | "
        ) + to_string(solCount));*/


        if (solCount > 0)
        {
            bitset<128> start;
            for (int i = 0; i < solCount; i++)
            {
                model.set(GRB_IntParam_SolutionNumber, i);
                for (int j = 0; j < 128; j++)
                    if (round(iK[j].get(GRB_DoubleAttr_Xn)) == 1)
                        start[j] = 1;
                    else
                        start[j] = 0;
                counterMap[start]++;
            }
            return SOLUTION;
        }
        else
            return NOSOLUTION;
    }
}


STATUS MidSolutionCounter(int midround, int rounds, const bitset<293>& startcube, const
    bitset<549>& last, map<bitset<549>,int,CMPS<549>> & counterMap, float
    time, int thread)
{
    if (midround < 256)
    {
        logger(__func__ + string(" : midround should be larger than 256."));
        exit(-1);
    }

    if (midround > rounds)
    {
        logger(__func__ + string(" : midround should not be larger than rounds."));
        exit(-2);
    }


    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread);
    env.set(GRB_IntParam_PoolSearchMode, 2);//focus on finding n best solutions 
    env.set(GRB_IntParam_PoolSolutions, MAX); // try to find 2000000
    GRBModel model = GRBModel(env);

    vector<GRBVar> is(293);
    for (int i = 0; i < 293; i++)
        is[i] = model.addVar(0, 1, 0, GRB_BINARY);
    for (int i = 0; i < 293; i++)
        if (startcube[i])
            model.addConstr(is[i] == 1);

    for (int i = 0; i < 293 - 256; i++)
        model.addConstr(is[i] == 0);


    vector<GRBVar> iK(128);
    for (int i = 0; i < 128; i++)
        iK[i] = model.addVar(0, 1, 0, GRB_BINARY);


    vector<GRBVar> iV(128);

    vector<GRBVar> s = is;
    vector<GRBVar> K = iK;
    vector<GRBVar> V = iV;

    for (int r = midround; r < rounds; r++)
        update(model, s, K, V, r);

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

    GRBLinExpr nk = 0;
    for (int i = 0; i < 128; i++)
        nk += iK[i];
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
        /*logger(__func__ + string(" : ") + to_string(rounds) + string(" | "
        ) + to_string(time) + string(" | "
        ) + to_string(solCount));*/

        if (solCount == 0)
            return NOSOLUTION;
        else
        {
            bitset<549> tmp;
            for (int i = 0; i < solCount; i++)
            {
                model.set(GRB_IntParam_SolutionNumber, i);
                for (int j = 0; j < 128; j++)
                    if (round(iK[j].get(GRB_DoubleAttr_Xn)) == 1)
                        tmp[j] = 1;
                    else
                        tmp[j] = 0;

                for (int j = 0; j < 293; j++)
                {
                    if (round(is[j].get(GRB_DoubleAttr_Xn)) == 1)
                        tmp[j + 256] = 1;
                    else
                        tmp[j + 256] = 0;
                }

                counterMap[tmp]++;
            }

            return SOLUTION;
        }
    }
}

