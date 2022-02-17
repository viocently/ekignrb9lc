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
void kreyviumTrackCore(bitset<288 + 256>& trackRep, int i1, int i2, int i3, int i4, int i5)
{
	bool tmp = trackRep[i3] && trackRep[i4];
	trackRep[i5] = trackRep[i1] || trackRep[i2] || tmp || trackRep[i5];
}

void kreyviumTrack(bitset<128>& cube, vector<bitset<288 + 256>>& trackReps)
{
	trackReps.clear();

	// start state
	bitset<288 + 256> start;
	for (int i = 0; i < 93; i++)
		start[i] = 1;
	for (int i = 0; i < 128; i++)
		if (cube[i])
		{
			start[93 + i] = 1;
			start[288 + 127 - i] = 1;
		}
	for (int i = 93 + 128; i < 287; i++)
		start[i] = 1;

	for (int i = 0; i < 128; i++)
		start[288 + 128 + i] = 1;

	// iterate state til all 1
	bitset<288 + 256> state = start;
	while (state.count() != (288 + 128 + cube.count()) )
	{
		trackReps.emplace_back(state);

		kreyviumTrackCore(state, 65, 170, 90, 91, 92);
		kreyviumTrackCore(state, 161, 263, 174, 175, 176);
		kreyviumTrackCore(state, 242, 68, 285, 286, 287);

		state[287] = 1;
		state[92] = state[92] || state[288];

		bitset<288 + 256> tmp(state);
		for (int j = 0; j < 288; j++)
			state[(j + 1) % 288] = tmp[j];

		// LFSR
		for (int j = 0; j < 128; j++)
		{
			state[288 + j] = tmp[288 + ((j + 1) % 128)];
			state[288 + 128 + j] = tmp[288 + 128 + ((j + 1) % 128)];
		}
	}
}

// track variables
#define IsVars1(x) (trackRep[x])
#define IsVars2(x,y) ( (trackRep[x] || trackRep[y] ) && (nonZero[x] && nonZero[y]) )
void kreyviumVarsTrackCore(bitset<288 + 256>& trackRep, int i1, int i2, int i3, int i4, int i5,bitset<288+256> & nonZero)
{
	bool tmp = IsVars2(i3, i4);
	trackRep[i5] = IsVars1(i1) || IsVars1(i2) || tmp || IsVars1(i5);
}

void kreyviumVarsTrack(bitset<128>& cube, vector<bitset<288 + 256>>& trackReps, vector<bitset<288 + 256>>& nonZeros)
{
	trackReps.clear();

	// start state
	bitset<288 + 256> start;
	for (int i = 0; i < 128; i++)
		if (cube[i])
		{
			start[93 + i] = 1;
			start[288 + 127-i] = 1;
		}


	// iterate state til all 1
	bitset<288 + 256> state = start;
	int r = 0;
	while (state.count() != (288 + cube.count()) )
	{
		trackReps.emplace_back(state);

		if (r < nonZeros.size())
		{
			kreyviumVarsTrackCore(state, 65, 170, 90, 91, 92,nonZeros[r]);
			kreyviumVarsTrackCore(state, 161, 263, 174, 175, 176,nonZeros[r]);
			kreyviumVarsTrackCore(state, 242, 68, 285, 286, 287,nonZeros[r]);
		}
		else
		{
			bitset<288 + 256> tmp;
			for (int i = 0; i < 288; i++)
				tmp[i] = 1;
			for (int i = 0; i < 128; i++)
				if (cube[127 - ((i + r) % 128)])
					tmp[288 + i] = 1;
			kreyviumVarsTrackCore(state, 65, 170, 90, 91, 92, tmp);
			kreyviumVarsTrackCore(state, 161, 263, 174, 175, 176, tmp);
			kreyviumVarsTrackCore(state, 242, 68, 285, 286, 287, tmp);
		}

		state[92] = state[92] || state[288];

		bitset<288 + 256> tmp(state);
		for (int j = 0; j < 288; j++)
			state[(j + 1) % 288] = tmp[j];

		// LFSR
		for (int j = 0; j < 128; j++)
		{
			state[288 + j] = tmp[288 + ((j + 1) % 128)];
			state[288 + 128 + j] = tmp[288 + 128 + ((j + 1) % 128)];
		}

		r++;
	}
}