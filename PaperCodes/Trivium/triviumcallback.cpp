#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"trivium.h"
#include"triviumtrack.h"
#include"newtrivium.h"
#include"utrivium.h"

using namespace std;
const int MAX = 200000000;

// for back expand
class backcallback : public GRBCallback
{
public:
    vector<bitset<288> >* psolutions = NULL;
    vector<GRBVar> solvars;
    bitset<288> last;
    int rounds;
    int midround;

    backcallback(vector<GRBVar>& solvars, vector<bitset<288> >& solutions, int rounds, int midround, const bitset<288>& last)
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
                bitset<288> sol;
                for (int i = 0; i < 288; i++)
                    if (round(getSolution(solvars[i])) == 1)
                        sol[i] = 1;

                // add lazy constr
                GRBLinExpr excludeCon = 0;
                for (int i = 0; i < 288; i++)
                    if (sol[i] == 1)
                        excludeCon += (1 - solvars[i]);
                    else
                        excludeCon += solvars[i];
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

STATUS uMidSolutionCallBack(int rounds, bitset<80>& cube, const bitset<288>& last, vector<bitset<288> >& sols, float time, int threads, int midround)
{


    vector<bitset<288>> nonZeroTrackReps;
    triviumTrack(cube, nonZeroTrackReps);


    vector<bitset<288>> varsTrackReps;
    triviumVarsTrack(cube, varsTrackReps, nonZeroTrackReps);


    // close standard output
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_LazyConstraints, 1);
    env.set(GRB_IntParam_MIPFocus, 1);
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
        uupdate(model, s, nonZeroFlag, varsFlag);
    }

    // mid constraints
    bitset<288> midNonZeroFlag;
    bitset<288> midVarsFlag;
    setRoundFlags(midround, cube, midNonZeroFlag, midVarsFlag, nonZeroTrackReps, varsTrackReps);

    vector<GRBVar> mids(288);
    for (int i = 0; i < 288; i++)
        mids[i] = model.addVar(0, 1, 0, GRB_BINARY);

    GRBLinExpr nk = 0;
    for (int i = 0; i < 288; i++)
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

    for (int i = 0; i < 288; i++)
        if (last[i])
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);

    model.setObjective(nk, GRB_MAXIMIZE);


    // set callback function
    backcallback cb = backcallback(mids, sols, rounds, midround, last);
    model.setCallback(&cb);

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

