#include<iostream>
#include<fstream>
#include<vector>
#include<bitset>
#include<algorithm>
#include<thread>
#include<mutex>
#include<future>
#include<sstream>
#include"log.h"
#include"deg.h"
#include"thread_pool.h"
#include"kreyvium.h"
#include"kreyviumtrack.h"
#include"newkreyvium.h"
#include"kreyviumcallback.h"
#include"kreyviumframework.h"

using namespace std;
using namespace thread_pool;
mutex stream_mutex;
mutex map_mutex;
mutex term_mutex;
int THREAD_NUM = 64;

void filterMap( map<bitset<544>, int, CMPS<544>> & mp )
{ 
	int size0 = mp.size();
	map< bitset<544>, int, CMPS<544>> tmp( mp );
    mp.clear();

    for ( auto it : tmp )
    {
		if ( it.second % 2 == 1 )
			mp[it.first] = it.second;
    }
	int size1 = mp.size();

	logger( __func__ + string(": ") + to_string( size0 ) + string( "\t" ) +
			to_string( size1 ) );
}

void filterMap(map<bitset<128>, int, CMPS<128>>& mp)
{
	int size0 = mp.size();
	map< bitset<128>, int, CMPS<128>> tmp(mp);
	mp.clear();

	for (auto& it : tmp)
	{
		if (it.second % 2 == 1)
			mp[it.first] = it.second;
	}
	int size1 = mp.size();

	logger(__func__ + string(": ") + to_string(size0) + string("\t") +
		to_string(size1));
}

void printSol( int rounds, const bitset<544> & vector, const map<bitset<128>,
		int, CMPS<128>> & solutions )  
{
	logger( __func__ + string(": " ) + to_string( rounds ) + string("\t") +
			vector.to_string()  );

	string path = string ( "TERM/" ) + to_string( rounds ) + string ( ".txt" ); 
	ofstream os;
	os.open( path, ios::out | ios::app );
	os << rounds << "\t" << vector << endl;
	for ( auto it : solutions )
		os << it.first << "\t" << it.second << "\n";
	os << endl;
	os.close();
}


void filterTerms( const bitset<128> & cube, int rounds, vector<bitset<544>> & terms )
{ 
	int size0 = terms.size();
    vector<bitset<544>> tmp( terms.begin(), terms.end() );
    terms.clear();
    for ( auto & it : tmp )
    {
        auto d = computeDegree( cube, rounds, it );
        if ( d >= cube.count() )
			terms.push_back( it );
    }
	int size1 = terms.size();

	logger( __func__ + string(": ") + to_string( size0 ) + string( "\t" ) +
			to_string( size1 ) );
}

// this function is only called from SolutionFunc to get solutions from counterMap
void SecondSolutionFunc(int midround, bitset<128>& cube, map<bitset<256+288>, vector<bitset<256+288> >, CMPS<256+288> >& counterMap, map<bitset<128>, int, CMPS<128> >& solutions, vector<bitset<256+288>>& expandTerms, float time, int singlethread)
{
	solutions.clear();

	// first sort out midsolutions
	vector<bitset<256+288> > midsolutions;
	for (auto& x : counterMap)
	{
		if (x.second.size() % 2 == 1)
			midsolutions.emplace_back(x.first);
		else
		{
			bitset<256+288> midstate = x.second[0];
			int midsols = 0;
			auto status = newMidSolutionCounter(midround, cube, midstate, midsols, time, singlethread);
			if (status == EXPAND || status == NOSOLUTION)
			{
				logger(__func__ + string(" : midsolution count error."));
				exit(-1);
			}

			int sols = x.second.size();
			if ((sols / midsols) % 2 == 1)
				midsolutions.emplace_back(x.first);
		}
	}

	// get final solutions from midsolutions
	for (auto& x : midsolutions)
	{
		auto status = MidSolutionCounter(midround, cube, x, solutions, 240, 2);
		if (status == NOSOLUTION)
		{
			logger(__func__ + string(" : get final solutions error."));
			exit(-1);
		}
		else if (status == EXPAND)
		{
			lock_guard<mutex> guard(stream_mutex);
			expandTerms.emplace_back(x);
		}
	}
}


STATUS SolutionFunc(int rounds, int midround, bitset<128>& cube, bitset<288+256>& last, map<bitset<128>, int, CMPS<128>>& solutions, vector<bitset<256+288>>& expandTerms, float time, int singlethread)
{
	if (midround > rounds || rounds < 0 || midround < 0)
	{
		logger(__func__ + string(" : parameter error."));
		exit(-1);
	}

	map<bitset<256+288>, vector<bitset<256+288>>, CMPS<256+288> > counterMap;
	auto status = newMidSolutionCounter(rounds, cube, last, counterMap, time, singlethread, midround);

	if (status == SOLUTION)
		SecondSolutionFunc(midround, cube, counterMap, solutions, expandTerms, time, singlethread);

	return status;
}

// below is threadworker
void SolutionSearcherWorker(bitset<256+288>& vec, int ROUND, bitset<128>& cube,
	vector<bitset<256+288>>& layerTerms, float time, int singlethread, int midround, map<bitset<128>, int, CMPS<128>>& totalSolutions, vector<bitset<256+288>>& expandTerms)
{
	map<bitset<128>, int, CMPS<128>> solutions;

	auto status = SolutionFunc(ROUND, midround, cube, vec, solutions, expandTerms, time, singlethread);

	if (status == EXPAND)
	{
		lock_guard<mutex> guard(map_mutex);
		layerTerms.push_back(vec);
	}
	else if (status == SOLUTION)
	{
		lock_guard<mutex> guard(stream_mutex);
		for (auto& x : solutions)
			totalSolutions[x.first] = (totalSolutions[x.first] + x.second) % 2;
	}
}

// below is for testing intergrity
template<int N>
bool vecCmp(vector<bitset<N> > _Left, vector<bitset<N> > _Right)
{
	if (_Left.size() == _Right.size())
	{
		for (int i = 0; i < _Left.size(); i++)
		{
			if (_Left[i] == _Right[i])
				continue;
			else
				return false;
		}
		return true;
	}
	return false;
}



int main()
{

	/* change cubeIndex to test other cube */
	vector<int> cubeIndex{ 0, 1, 2, 3, 4, 5,6,7,8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
		18, 19, 20, 21, 22,/**/  23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
		37, 38, /**/ 39, 40, 41, 42, 43, /**/ 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
		57, 58, 59, 60, 61, /**/ 62, 63, 64, 65, 67, 68, 69, 70, 71, 74, 75, 76, 77, 79,
		80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
		98, 99, 100, 102, 103, 104, 105, 107, 108, 111, 112, 113, 114, 115, 116,
		117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127 };

	bitset<128> cube{ 0 };
	for ( auto & it : cubeIndex )
		cube[it] = 1;

	int midround = 120;
	ThreadPool thread_pool{};


	int ROUND = 895;
	int step = 300;
	int expandStep = 20;
	float expandTime = 14400;

	float time = 0;

	vector<bitset<256+288>> terms;
	vector<bitset<256+288>> layerTerms;

	BackExpandPolynomial(step, terms);

	ROUND = ROUND - step;
	filterTerms(cube, ROUND, terms);

	logger(__func__ + string("First expand"));

	int singlethread = 2;

	map<bitset<128>, int, CMPS<128>> totalSolutions;
	vector<bitset<256+288> > expandTerms;
	

	if (ROUND < CALLBACKBOUND)
		newExpandFunc(cube, terms, ROUND, expandStep, thread_pool, expandTime, singlethread);
	else
		expandFunc(cube, terms, ROUND, expandTime, singlethread);

	while (true)
	{
		if (ROUND > 600)
			time = 60;
		else if (ROUND > 500)
			time = 120;
		else if (ROUND > 400)
			time = 180;
		else if (ROUND > 300)
			time = 360;
		else if (ROUND > 200)
			time = 720;
		else
			time = 0;

		

		layerTerms.clear();
		vector<std::future<void>> futures;
		for (auto& it : terms)
			futures.emplace_back(thread_pool.Submit(SolutionSearcherWorker, it,
				ROUND, ref(cube), ref(layerTerms), time, singlethread, midround, ref(totalSolutions), ref(expandTerms)));
		for (auto& it : futures)
			it.get();

		filterMap(totalSolutions);


		logger(string("layerTermsSize ") + to_string(layerTerms.size()));

		logger(string("expandTermsSize ") + to_string(expandTerms.size()));

		if (layerTerms.size() == 0 && expandTerms.size() == 0)
			break;
		else if (layerTerms.size() == 0)
		{
			logger("Start to handle expandTerms.");
			ROUND = midround;
			midround = ROUND / 2;
			layerTerms = expandTerms;
			filterTerms(cube, ROUND, layerTerms);
			expandTerms.clear();
		}

		if (ROUND < 10)
		{
			logger("ROUND is to small.");
			exit(-1);
		}

		terms = layerTerms;


		if (ROUND < FASTBOUND)
			expandStep = 20;

		if (ROUND - expandStep < midround)
			expandStep = ROUND - midround;

		// second Expand
		if (ROUND < CALLBACKBOUND)
			newExpandFunc(cube, terms, ROUND, expandStep, thread_pool, expandTime, singlethread);
		else
			expandFunc(cube, terms, ROUND, expandTime, singlethread);

		// if ROUND == midround,combine expandterms with terms since they are all in midround
		if (ROUND == midround)
		{
			logger("current ROUND : " + to_string(ROUND));
			logger("start combine expandTerms with terms.");
			map<bitset<256+288>, int, CMPS<256+288>> termCounter;
			for (auto& x : terms)
				termCounter[x] ++;
			for (auto& x : expandTerms)
				termCounter[x] ++;

			expandTerms.clear();
			terms.clear();
			for (auto& x : termCounter)
				if (x.second % 2 == 1)
					terms.emplace_back(x.first);

			midround = ROUND / 2;
		}

	}
	cout << "Success!" << endl;


	cout << "Output to file." << endl;

	// output result to file
	string path = string("TERM/") + string("res.txt");
	ofstream os;
	os.open(path, ios::out | ios::app);
	for (auto& it : totalSolutions)
	{
		bitset<128> term = it.first;
		if (term.count() == 0)
		{
			os << "1" << endl;
		}
		else
		{
			for (int i = 0; i < 128; i++)
				if (term[i])
					os << "k" << i;
			os << endl;
		}
	}
	os << endl;
	os.close();


	cout << "Output to file finished." << endl;
}
