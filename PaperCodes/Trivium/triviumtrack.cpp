#include<map>
#include<cmath>
#include<vector>
#include<bitset>
#include<algorithm>
#include<string>
#include"gurobi_c++.h" 
#include"node.h"
#include"log.h"

using namespace std;

const int MAX = 200000000; // the maximum value of PoolSearchMode, P625


// 288 iv key
// track non-zero variables
void triviumTrackCore(bitset<288>& trackRep, int i1, int i2, int i3, int i4, int i5)
{
	bool tmp = trackRep[i3] && trackRep[i4];
	trackRep[i5] = trackRep[i1] || trackRep[i2] || tmp || trackRep[i5];
}

void triviumTrack(bitset<80>& cube, vector<bitset<288>>& trackReps)
{
	trackReps.clear();

	// start state
	bitset<288> start;
	for (int i = 0; i < 80; i++)
		start[i] = 1;
	for (int i = 0; i < 80; i++)
		if (cube[i])
			start[93 + i] = 1;
	for (int i = 285; i < 288; i++)
		start[i] = 1;

	// iterate state til all 1
	bitset<288> state = start;
	while (state.count() != 288)
	{
		trackReps.emplace_back(state);

		triviumTrackCore(state, 65, 170, 90, 91, 92);
		triviumTrackCore(state, 161, 263, 174, 175, 176);
		triviumTrackCore(state, 242, 68, 285, 286, 287);

		bitset<288> tmp(state);
		for (int j = 0; j < 288; j++)
			state[(j + 1) % 288] = tmp[j];
	}
}

// track variables
#define IsVars1(x) (trackRep[x])
#define IsVars2(x,y) ( (trackRep[x] || trackRep[y] ) && (nonZero[x] && nonZero[y]) )
void triviumVarsTrackCore(bitset<288>& trackRep, int i1, int i2, int i3, int i4, int i5, bitset<288>& nonZero)
{
	bool tmp = IsVars2(i3, i4);
	trackRep[i5] = IsVars1(i1) || IsVars1(i2) || tmp || IsVars1(i5);
}

void triviumVarsTrack(bitset<80>& cube, vector<bitset<288>>& trackReps, vector<bitset<288>>& nonZeros)
{
	trackReps.clear();

	// start state
	bitset<288> start;
	for (int i = 0; i < 80; i++)
		if (cube[i])
			start[i + 93] = 1;


	// iterate state til all 1
	bitset<288> state = start;
	int r = 0;
	while (state.count() != 288)
	{
		trackReps.emplace_back(state);

		if (r < nonZeros.size())
		{
			triviumVarsTrackCore(state, 65, 170, 90, 91, 92, nonZeros[r]);
			triviumVarsTrackCore(state, 161, 263, 174, 175, 176, nonZeros[r]);
			triviumVarsTrackCore(state, 242, 68, 285, 286, 287, nonZeros[r]);
		}
		else
		{
			bitset<288> tmp;
			for (int i = 0; i < 288; i++)
				tmp[i] = 1;
			triviumVarsTrackCore(state, 65, 170, 90, 91, 92, tmp);
			triviumVarsTrackCore(state, 161, 263, 174, 175, 176, tmp);
			triviumVarsTrackCore(state, 242, 68, 285, 286, 287, tmp);
		}


		bitset<288> tmp(state);
		for (int j = 0; j < 288; j++)
			state[(j + 1) % 288] = tmp[j];

		r++;
	}
}