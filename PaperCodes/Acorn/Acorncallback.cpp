#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"Acorn.h"
#include"Acorncallback.h"
#include"newAcorn.h"
#include"Acorntrack.h"
#include"uAcorn.h"

using namespace std;
const int MAX = 200000000;

class backcallback : public GRBCallback
{
public:
    vector<bitset<256+293> >* psolutions = NULL;
    vector<GRBVar> solvars;
    bitset<256+293> last;
    int rounds;
    int midround;

    backcallback(vector<GRBVar>& solvars, vector<bitset<256+293> >& solutions, int rounds, int midround, const bitset<256+293>& last)
    {
        this->solvars = solvars;
        this->psolutions = &solutions;
        this->last = last;
        this->rounds = rounds;
        this->midround = midround;
    }

protected:
    void callback()
    {
        try {
            if (where == GRB_CB_MIPSOL)
            {
                // read solution
                bitset<256+293> sol;
                for (int i = 0; i < 128; i++)
                    if (round(getSolution(solvars[i])) == 1)
                        sol[i] = 1;

                for (int i = 0; i < 293; i++)
                    if (round(getSolution(solvars[i + 256])) == 1)
                        sol[i + 256] = 1;

                // add lazy constr
                GRBLinExpr excludeCon = 0;
                for (int i = 0; i < 128; i++)
                    if (sol[i] == 1)
                        excludeCon += (1 - solvars[i]);
                    else
                        excludeCon += solvars[i];

                for (int i = 0; i < 293; i++)
                    if (sol[i+256] == 1)
                        excludeCon += (1 - solvars[i+256]);
                    else
                        excludeCon += solvars[i+256];

                addLazy(excludeCon >= 1);

                // check its parity
                int pathcnt = 0;
                STATUS status;
                status = MidSolutionCounter(midround, rounds, sol, last, pathcnt, 120, 1);
                if (status == EXPAND)
                {
                    logger("In CallBack : parity check failed.");
                    exit(-1);
                }

                if (pathcnt % 2 == 1)
                {
                    // save solution
                    cout << "back find one solution : " << sol << endl;
                    (*psolutions).emplace_back(sol);

                }

            }
        }
        catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        }
        catch (...) {
            cout << "Error during callback" << endl;
        }
    }
};


STATUS uMidSolutionCallBack(int rounds, bitset<293>& startcube, bitset<293 + 256>& last, vector<bitset<293 + 256>>& sols, float time, int threads, int midround, int coRound)
{
    if (midround > rounds || rounds < 0 || midround < 0 || rounds < 256 || midround < 256 || coRound < 256)
    {
        logger(__func__ + string(":parameter error."));
        exit(-1);
    }

    sols.clear();

    vector<bitset<293>> nonZeroTrackReps;
    vector<bitset<293>> varsTrackReps;

    acornNonZeroTrack(nonZeroTrackReps);
    acornVarsTrack256(varsTrackReps, startcube,nonZeroTrackReps);

    GRBEnv env = GRBEnv();

    // close standard output
    // env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_LazyConstraints, 1);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_Presolve, 2);
    //env.set(GRB_IntParam_PoolSearchMode, 2);
    //env.set(GRB_IntParam_PoolSolutions, MAX);
    //env.set(GRB_IntParam_MIPFocus, 1);

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
    for (int r = 256; r < coRound; r++)
    {
        bitset<293> nonZeroFlag;
        bitset<293> varsFlag;
        setRoundFlags(r, nonZeroFlag, varsFlag, nonZeroTrackReps, varsTrackReps);
        uupdate(model, s, nonZeroFlag,varsFlag);

    }

    // mid constraints
    bitset<293> midNonZeroFlag;
    bitset<293> midVarsFlag;
    setRoundFlags(coRound, midNonZeroFlag, midVarsFlag, nonZeroTrackReps, varsTrackReps);

    vector<GRBVar> mids2 = s;


    vector<GRBVar> mids(293);
    for (int i = 0; i < 293; i++)
        mids[i] = model.addVar(0, 1, 0, GRB_BINARY);

    GRBLinExpr nk = 0;
    for (int i = 0; i < 293; i++)
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
    vector<GRBVar> tmps;
    vector<GRBVar> tmpK;
    for (int r = coRound; r < rounds; r++)
    {
        if (r == midround)
        {
            tmps = s;
            tmpK = K;
        }
        for (int i = 0; i < 128; i++)
            K[i].set(GRB_DoubleAttr_VarHintVal, 0);
        update(model, s, K, V, r);
    }

    // output constraints
    // key + iv + 288
    for (int i = 0; i < 293; i++)
        if (last[i + 256])
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

    vector<GRBVar> midstate;
    for (int i = 0; i < 128; i++)
        midstate.emplace_back(tmpK[i]);
    for (int i = 0; i < 128; i++)
    {
        GRBVar paddingVar;
        midstate.emplace_back(paddingVar);
    }

    for (int i = 0; i < 293; i++)
        midstate.emplace_back(tmps[i]);
    backcallback cb = backcallback(midstate, sols, rounds, midround, last);
    model.setCallback(&cb);

    for (int i = 0; i < 293; i++)
        tmps[i].set(GRB_IntAttr_BranchPriority, 2);

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
            return SOLUTION;
        }
        else
            return NOSOLUTION;

    }
}