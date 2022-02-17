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
#include"newgrain.h"
#include"ugrain.h"

using namespace std;
const int MAX = 200000000;

class mycallback : public GRBCallback
{
public:
	vector<bitset<256> >* psolutions = NULL;
	vector<GRBVar> solvars;
	bitset<256> last;
	int rounds;
	int midround;

    mycallback(vector<GRBVar>& solvars, vector<bitset<256> >& solutions, int rounds, int midround, const bitset<256>& last)
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
                bitset<256> sol;
                for (int i = 0; i < 256; i++)
                    if (round(getSolution(solvars[i])) == 1)
                        sol[i] = 1;

                // add lazy constr
                GRBLinExpr excludeCon = 0;
                for (int i = 0; i < 256; i++)
                    if (sol[i] == 1)
                        excludeCon += (1 - solvars[i]);
                    else
                        excludeCon += solvars[i];
                addLazy(excludeCon >= 1);

                // check its parity
                int pathcnt = 0;
                STATUS status;
                status = MidSolutionCounter(midround, rounds, sol, last, pathcnt,120, 1);
                if (status == EXPAND)
                {
                    logger("In CallBack : parity check failed.");
                    exit(-1);
                }

                if (pathcnt % 2 == 1)
                {
                    // save solution
                    cout << "back find solution:" << sol << endl;
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


STATUS uMidSolutionCallBack(int rounds, bitset<96>& cube,
    bitset<256>& last, int midround, vector<bitset<256>>& sols, float
    time, int thread, int coRound)
{
    if (midround > rounds || rounds < 0 || midround < 0)
    {
        logger(__func__ + string(":parameter error."));
        exit(-1);
    }

 
    bitset<256> flag1;
    flag1.set();

    sols.clear();


    vector<bitset<256> > nonZeroTrackReps;
    grainTrack(cube, nonZeroTrackReps);


    vector<bitset<256> > varsTrackReps;
    grainVarsTrack(cube, varsTrackReps, nonZeroTrackReps);

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, thread);
    env.set(GRB_IntParam_LazyConstraints, 1);
    // env.set(GRB_IntParam_MIPFocus, 1);
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

    for (int r = 0; r < coRound; r++)
    {
        uupdate(model, works, workb, r < nonZeroTrackReps.size() ? nonZeroTrackReps[r] : flag1, r < varsTrackReps.size() ? varsTrackReps[r] : flag1);
    }

    vector<GRBVar> midworks(128);
    vector<GRBVar> midworkb(128);

    vector<GRBVar> midworks2 = works;
    vector<GRBVar> midworkb2 = workb;

    // mid constraints
    bitset<256> midNonZeroFlag = (coRound < nonZeroTrackReps.size()) ? nonZeroTrackReps[coRound] : flag1;
    bitset<256> midVarsFlag = (coRound < varsTrackReps.size()) ? varsTrackReps[coRound] : flag1;

    GRBLinExpr nk = 0;
    for (int i = 0; i < 128; i++)
    {
        midworkb[i] = model.addVar(0, 1, 0, GRB_BINARY);
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
    vector<GRBVar> tmps;
    vector<GRBVar> tmpb;
    for (int r = coRound; r < rounds; r++)
    {
        if (r == midround)
        {
            tmps = works;
            tmpb = workb;
        }
        update(model, works, workb);
    }

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

    vector<GRBVar> midstate(tmpb);
    for (auto& x : tmps)
        midstate.emplace_back(x);

    mycallback cb = mycallback(midstate, sols, rounds, midround, last);
    model.setCallback(&cb);

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
