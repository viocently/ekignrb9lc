#include<iostream>
#include<bitset>
#include<vector>

using namespace std;


bool grainTrackfuncG(bitset<256>& trackRep)
{
	if (trackRep[0] || trackRep[26] || trackRep[56] || trackRep[91] || trackRep[96]
		|| (trackRep[3] && trackRep[67]) || (trackRep[11] && trackRep[13]) || (trackRep[17] && trackRep[18])
		|| (trackRep[27] && trackRep[59]) || (trackRep[40] && trackRep[48]) || (trackRep[61] && trackRep[65]) 
		|| (trackRep[68] && trackRep[84]) || (trackRep[88] && trackRep[92] && trackRep[93] && trackRep[95])
		|| (trackRep[22] && trackRep[24] && trackRep[25]) || (trackRep[70] && trackRep[78] && trackRep[82]))
		return true;
	else
		return false;
}

bool grainTrackfuncF(bitset<256>& trackRep)
{
	if (trackRep[128+0] || trackRep[128+7] || trackRep[128+38] || trackRep[128+70] || trackRep[128+81] || trackRep[128+96])
		return true;
	else
		return false;
}

bool grainTrackfuncH(bitset<256>& trackRep)
{
	if ((trackRep[12] && trackRep[8 + 128]) || (trackRep[13 + 128] && trackRep[20 + 128]) || (trackRep[95] && trackRep[42 + 128])
		|| (trackRep[60 + 128] && trackRep[79 + 128]) || (trackRep[12] && trackRep[95] && trackRep[94 + 128]))
		return true;
	else
		return false;
}

bool grainTrackfuncZ(bitset<256>& trackRep)
{
	bool tmpH = grainTrackfuncH(trackRep);
	if (tmpH || trackRep[93 + 128] || trackRep[2] || trackRep[15] || trackRep[36] || trackRep[45]
		|| trackRep[64] || trackRep[73] || trackRep[89])
		return true;
	else
		return false;
}


void grainTrack(bitset<96>& cube, vector<bitset<256>>& trackReps)
{
	trackReps.clear();

	// start state
	bitset<256> start;
	for (int i = 0; i < 128; i++)
		start[i] = 1;
	for (int i = 0; i < 96; i++)
		if (cube[i])
			start[i + 128] = 1;

	for (int i = 96; i < 127; i++)
		start[i+128] = 1;


	// iterate state til all 1
	bitset<256> state = start;
	while (state.count() != 256)
	{
		trackReps.emplace_back(state);
		bool g = grainTrackfuncG(state);
		bool z = grainTrackfuncZ(state);
		bool f = grainTrackfuncF(state);
		bool s0 = state[0+128];

		for (int i = 0; i < 127; i++)
		{
			state[i] = state[i + 1];
			state[i + 128] = state[i + 129];
		}

		state[127] = g || s0 || z;
		state[255] = f || z;
	}
}


#define IsVars1(x) (trackRep[x])
#define IsVars2(x,y) ( (trackRep[x] || trackRep[y] ) && (nonZero[x] && nonZero[y]) )
#define IsVars3(x,y,z) ( (trackRep[x] || trackRep[y] || trackRep[z]) && (nonZero[x] && nonZero[y] && nonZero[z]) )
#define IsVars4(x,y,z,t) ( (trackRep[x] || trackRep[y] || trackRep[z] || trackRep[t]) && (nonZero[x] && nonZero[y] && nonZero[z] && nonZero[t]) )

bool grainVarsTrackfuncG(bitset<256>& trackRep,bitset<256>& nonZero)
{
	if (IsVars1(0) || IsVars1(26) || IsVars1(56) || IsVars1(91) || IsVars1(96) 
		|| IsVars2(3,67) || IsVars2(11,13) || IsVars2(17,18) || IsVars2(27,59)
		|| IsVars2(40,48) || IsVars2(61,65) || IsVars2(68,84) || IsVars4(88,92,93,95)
		|| IsVars3(22,24,25) || IsVars3(70,78,82))
		return true;
	else
		return false;
}

bool grainVarsTrackfuncF(bitset<256>& trackRep, bitset<256>& nonZero)
{
	if (IsVars1(128+0) || IsVars1(128+7) || IsVars1(128+38) || IsVars1(128+70) 
		|| IsVars1(128+81) || IsVars1(128+96) )
		return true;
	else
		return false;
}

bool grainVarsTrackfuncH(bitset<256>& trackRep, bitset<256>& nonZero)
{
	if (IsVars2(12,128+8) || IsVars2(128+13,128+20) || IsVars2(95,128+42) 
		|| IsVars2(128+60,128+79) || IsVars3(12,95,128+94))
		return true;
	else
		return false;
}

bool grainVarsTrackfuncZ(bitset<256>& trackRep, bitset<256>& nonZero)
{
	bool tmpH = grainVarsTrackfuncH(trackRep,nonZero);
	if (tmpH || IsVars1(93+128) || IsVars1(2) || IsVars1(15) || IsVars1(36)
		|| IsVars1(45) || IsVars1(64) || IsVars1(73) || IsVars1(89) )
		return true;
	else
		return false;
}


void grainVarsTrack(bitset<96>& cube, vector<bitset<256>>& trackReps,vector<bitset<256>> & nonZeros)
{
	trackReps.clear();

	// start state
	bitset<256> start;
	for (int i = 0; i < 96; i++)
		if (cube[i])
			start[i + 128] = 1;

	// iterate state til all 1
	bitset<256> state = start;
	int r = 0;
	while (state.count() != 256)
	{
		trackReps.emplace_back(state);
		bool g, z, f;
		if (r < nonZeros.size())
		{
			g = grainVarsTrackfuncG(state,nonZeros[r]);
			z = grainVarsTrackfuncZ(state,nonZeros[r]);
			f = grainVarsTrackfuncF(state, nonZeros[r]);
		}
		else
		{
			bitset<256> tmp;
			tmp.set();
			g = grainVarsTrackfuncG(state, tmp);
			z = grainVarsTrackfuncZ(state, tmp);
			f = grainVarsTrackfuncF(state, tmp);
		}
		bool s0 = state[0 + 128];

		for (int i = 0; i < 127; i++)
		{
			state[i] = state[i + 1];
			state[i + 128] = state[i + 129];
		}

		state[127] = g || s0 || z;
		state[255] = f || z;
		
		r++;
	}
}