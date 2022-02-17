#include<iostream>
#include<bitset>
#include<vector>
#include"log.h"

using namespace std;

// non-zero flags
void xorFTrack(bitset<293>& nonZeroTrack, int i, int j, int k)
{
	nonZeroTrack[i] = nonZeroTrack[i] || nonZeroTrack[j] || nonZeroTrack[k];
}

void updateNonZeroTrack(bitset<293> & nonZeroTrack)
{
	for (int i = 0; i < 292; i++)
		nonZeroTrack[i] = nonZeroTrack[i + 1];
	nonZeroTrack[292] = 1;
}

int acornNonZeroTrack(vector<bitset<293>>& nonZeroTrackReps)
{
	bitset<293> nonZeroTrack;

	int r = 0; // start from 256
	while (nonZeroTrack.count() != 293)
	{
		nonZeroTrackReps.emplace_back(nonZeroTrack);
		updateNonZeroTrack(nonZeroTrack);
		r++;
	}
	return r;
}

// variable-containing flags
#define IsVars2(x,y) ( (varsTrack[x] || varsTrack[y] ) && (nonZeroTrack[x] && nonZeroTrack[y]) )
bool majVarsTrack(bitset<293>& varsTrack, int i, int j, int k,bitset<293> & nonZeroTrack)
{
	bool t1 = IsVars2(i, j);
	bool t2 = IsVars2(i, k);
	bool t3 = IsVars2(j, k);
	return t1 || t2 || t3;
}

bool chVarsTrack(bitset<293>& varsTrack, int i, int j, int k, bitset<293>& nonZeroTrack)
{
	bool t1 = IsVars2(i, j);
	bool t2 = IsVars2(i, k);
	return t1 || t2 ;
}

void xorFVarsTrack(bitset<293>& varsTrack, int i, int j, int k)
{
	varsTrack[i] = varsTrack[i] || varsTrack[j] || varsTrack[k];
}

bool ksgVarsTrack(bitset<293>& varsTrack,bitset<293>& nonZeroTrack)
{
	bool t1 = majVarsTrack(varsTrack, 235,61,193,  nonZeroTrack);
	bool t2 = chVarsTrack(varsTrack, 230, 111, 66,  nonZeroTrack);
	return t1 || t2 || varsTrack[12] || varsTrack[154];
}

bool fbkVarsTrack(bitset<293>& varsTrack, bitset<293>& nonZeroTrack)
{
	bool t1 = majVarsTrack(varsTrack, 244, 23, 160, nonZeroTrack);
	return t1 || varsTrack[0] || varsTrack[107] || varsTrack[196];
}

bool genMVarsTrack(bitset<128>& cube, int r)
{
	if (r < 128 || r >= 256)
	{
		logger(__func__ + string(" : parameter rounds is out of range."));
		exit(-1);
	}

	if (cube[r - 128] == 1)
		return true;
	else
		return false;
}

void updateVarsTrack(bitset<293>& varsTrack, bitset<128>& cube, int r, bitset<293>& nonZeroTrack)
{
	xorFVarsTrack(varsTrack, 289, 235, 230);
	xorFVarsTrack(varsTrack, 230, 196, 193);
	xorFVarsTrack(varsTrack, 193, 160, 154);
	xorFVarsTrack(varsTrack, 154, 111, 107);
	xorFVarsTrack(varsTrack, 107, 66, 61);
	xorFVarsTrack(varsTrack, 61, 23, 0);

	bool t1 = ksgVarsTrack(varsTrack, nonZeroTrack);
	bool t2 = fbkVarsTrack(varsTrack, nonZeroTrack);
	bool t3 = false;
	if (r >= 128 && r < 256 && (cube[r - 128] == 1))
		t3 = genMVarsTrack(cube, r);

	for (int i = 0; i < 292; i++)
		varsTrack[i] = varsTrack[i + 1];
	varsTrack[292] = t1 || t2 || t3;
}

// start from round 0
int acornVarsTrack0(vector<bitset<293>> & varsTrackReps, bitset<128> & cube,vector<bitset<293> > & nonZeroTrackReps)
{
	bitset<293> varsTrack;
	int r = 0;
	while (varsTrack.count() != 293)
	{
		bitset<293> nonZeroTrack;
		if (r < nonZeroTrackReps.size())
			nonZeroTrack = nonZeroTrackReps[r];
		else
			nonZeroTrack.set();

		varsTrackReps.emplace_back(varsTrack);
		updateVarsTrack(varsTrack, cube, r, nonZeroTrack);
		r++;
	}
	return r;
}

int acornVarsTrack256(vector<bitset<293>>& varsTrackReps, bitset<293> & startcube,vector<bitset<293> >& nonZeroTrackReps)
{
	bitset<293> varsTrack = startcube;
	bitset<128> cube;
	int r = 0;
	while (r < 256)
	{
		bitset<293> paddingflags;
		varsTrackReps.emplace_back(paddingflags);
		r++;
	}
	while (varsTrack.count() != 293)
	{
		bitset<293> nonZeroTrack;
		if (r < nonZeroTrackReps.size())
			nonZeroTrack = nonZeroTrackReps[r];
		else
			nonZeroTrack.set();

		varsTrackReps.emplace_back(varsTrack);
		updateVarsTrack(varsTrack, cube, r,nonZeroTrack);
		r++;
	}
	return r;
}

