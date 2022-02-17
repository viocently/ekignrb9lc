#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"
#include"kreyvium.h"
#include"kreyviumtrack.h"
#include"newkreyvium.h"
#include"ukreyvium.h"

using namespace std;
const int MAX = 200000000;

class mycallback : public GRBCallback
{
public:
    vector<bitset<288+256> >* psolutions = NULL;
    vector<GRBVar> solvars;
    bitset<288+256> last;
    int rounds;
    int midround;

    mycallback(vector<GRBVar>& solvars, vector<bitset<256+288> >& solutions, int rounds, int midround, const bitset<256+288>& last)
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
                bitset<256+288> sol;
                for (int i = 0; i < 256+288; i++)
                    if (round(getSolution(solvars[i])) == 1)
                        sol[i] = 1;
                cout << sol << endl;

                // add lazy constr
                GRBLinExpr excludeCon = 0;
                for (int i = 0; i < 256+288; i++)
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


STATUS MidSolutionCallBack(int rounds, bitset<128>& cube,const bitset<544> & last, vector<bitset<288 + 256> >& sols, float time, int threads, int midround)
{
    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_LazyConstraints, 1);
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
    model.addConstr(s[287] == 0);

    vector<GRBVar> mids;
    vector<GRBVar> midkey;
    vector<GRBVar> midiv;
    // update
    for (int i = 0; i < rounds; i++)
    {
        update(model, s, key, iv);
        if (i == midround)
        {
            mids = s;
            midkey = key;
            midiv = iv;
        }
    }

    // key + iv + 288
    for (int i = 0; i < 128; i++)
        model.addConstr(key[i] == last[i]);
    for (int i = 0; i < 128; i++)
        model.addConstr(iv[i] == last[i + 128]);
    for (int i = 0; i < 288; i++)
        model.addConstr(s[i] == last[i + 256]);

    // objective function
    /*GRBLinExpr nk = 0;
    for (int i = 0; i < 128; i++)
        nk += midkey[i];
    model.setObjective(nk, GRB_MAXIMIZE);*/
    GRBLinExpr nv = 0;
    for (int i = 0; i < 128; i++)
        nv += KEY[i];
    model.setObjective(nv, GRB_MAXIMIZE);

    if (time > 0)
        model.set(GRB_DoubleParam_TimeLimit, time);
    vector<GRBVar> midstate;
    for (int i = 0; i < 128; i++)
        midstate.emplace_back(midkey[i]);
    for (int i = 0; i < 128; i++)
        midstate.emplace_back(midiv[i]);
    for (int i = 0; i < 288; i++)
        midstate.emplace_back(mids[i]);

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

STATUS uMidSolutionCallBack(int rounds, bitset<128>& cube, bitset<544>& last, 
    vector<bitset<288 + 256>> & sols, float time, int threads, int midround)
{


    vector<bitset<544>> nonZeroTrackReps;
    kreyviumTrack(cube, nonZeroTrackReps);

    vector<bitset<544>> varsTrackReps;
    kreyviumVarsTrack(cube, varsTrackReps, nonZeroTrackReps);


    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads, threads);
    env.set(GRB_IntParam_LazyConstraints, 1);
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
        model.addConstr(IV[i] == iv[127 - i] + s[93 + i]);
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
        uupdate(model, s, key, iv, nonZeroFlag, varsFlag);
    }

    // mid constraints
    bitset<288 + 256> midNonZeroFlag;
    bitset<288 + 256> midVarsFlag;
    setRoundFlags(midround, cube, midNonZeroFlag, midVarsFlag, nonZeroTrackReps, varsTrackReps);


    vector<GRBVar> mids2 = s;
    vector<GRBVar> midkey2 = key;
    vector<GRBVar> midiv2 = iv;

    GRBLinExpr nk = 0;
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
        if (midNonZeroFlag[i])
        {
            model.addConstr(mids[i] >= s[i]);
            nk += mids[i] - s[i];
        }
        else
            model.addConstr(mids[i] == 0);

    for (int i = 0; i < 128; i++)
        if (midNonZeroFlag[i + 288])
        {
            model.addConstr(midiv[i] == iv[i]);
        }
        else
            model.addConstr(midiv[i] == 0);

    for (int i = 0; i < 128; i++)
    {
        model.addConstr(midkey[i] >= key[i]);
        nk += midkey[i] - key[i];
    }

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
        if (last[i + 256])
            model.addConstr(s[i] == 1);
        else
            model.addConstr(s[i] == 0);


    vector<GRBVar> midstate;
    for (int i = 0; i < 128; i++)
        midstate.emplace_back(midkey[i]);
    for (int i = 0; i < 128; i++)
        midstate.emplace_back(midiv[i]);
    for (int i = 0; i < 288; i++)
        midstate.emplace_back(mids[i]);

    mycallback cb = mycallback(midstate, sols, rounds, midround, last);
    model.setCallback(&cb);

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
            return SOLUTION;
        }
        else
            return NOSOLUTION;

    }


}
