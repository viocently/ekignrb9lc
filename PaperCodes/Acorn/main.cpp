#include<iostream>
#include<fstream>
#include<vector>
#include<bitset>
#include<algorithm>
#include<thread>
#include<mutex>
#include<future>
#include<sstream>
#include<algorithm>
#include<thread>
#include<string>
#include<malloc.h>
#include"log.h"
#include"thread_pool.h"
#include"Acorn.h"
#include"newAcorn.h"
#include"uAcorn.h"
#include"Acorncallback.h"
#include"Acornframework.h"


using namespace std;
using namespace thread_pool;
const int THREAD_NUM = 64;

mutex stream_mutex;
mutex map_mutex;
mutex term_mutex;

// this function is only called from SolutionFunc to get solutions from counterMap
void SecondSolutionFunc(int midround, bitset<293> & startcube, map<bitset<549>, int, CMPS<549> >& counterMap, 
	map<bitset<549>, int, CMPS<549> >& solutions, vector<bitset<549>>& expandTerms, float time, int singlethread)
{
	solutions.clear();

	// first sort out midsolutions
	vector<bitset<549> > midsolutions;
	for (auto& x : counterMap)
	{
		if (x.second % 2 == 1)
			midsolutions.emplace_back(x.first);
		else
		{
			bitset<293> midstate;
			for (int i = 0; i < 293; i++)
				midstate[i] = x.first[i + 256];

			int midsols = 0;
			auto status = newMidSolutionCounter(midround, startcube, midstate, midsols, time, singlethread);
			if (status == EXPAND || status == NOSOLUTION)
			{
				logger(__func__ + string(" : midsolution count error."));
				exit(-1);
			}

			int sols = x.second;
			if ((sols / midsols) % 2 == 1)
				midsolutions.emplace_back(x.first);
		}
	}

	// get final solutions from midsolutions
	for (auto& x : midsolutions)
	{
		auto status = MidSolutionCounter(256, midround, startcube, x, solutions, 120, 2);
		if (status == NOSOLUTION)
		{
			logger(__func__ + string(" : get final solutions error."));
			exit(-1);
		}
		else if (status == EXPAND)
		{
			logger(__func__ + string("expandTerm added : ") + x.to_string());

			lock_guard<mutex> guard(term_mutex);
			expandTerms.emplace_back(x);
		}
	}
}

STATUS SolutionFunc(int rounds, int midround, bitset<293>& startcube, bitset<549>& last, map<bitset<549>, int, CMPS<549>>& solutions, vector<bitset<549>>& expandTerms, float time, int singlethread)
{
	if (midround > rounds || rounds < 0 || midround < 0)
	{
		logger(__func__ + string(" : parameter error."));
		exit(-1);
	}

	map<bitset<549>, int, CMPS<549> > counterMap;
	auto status = newMidSolutionCounter(rounds, startcube, last, counterMap, time, singlethread, midround);

	if (status == SOLUTION)
		SecondSolutionFunc(midround, startcube, counterMap, solutions, expandTerms, time, singlethread);

	return status;
}

// below is threadworker
void SolutionSearcherWorker(bitset<549>& vec, int ROUND, bitset<293>& startcube,
	vector<bitset<549>>& layerTerms, float time, int singlethread, int midround, map<bitset<549>, int, CMPS<549>>& totalSolutions, vector<bitset<549>>& expandTerms)
{
	map<bitset<549>, int, CMPS<549>> solutions;

	auto status = SolutionFunc(ROUND, midround, startcube, vec, solutions, expandTerms, time, singlethread);

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

// this threadworker is used to solve terms in round 256
void SecondSolutionSearcherWorker(const bitset<549>& vec, int ROUND, const
	bitset<128>& cube, float time, int singlethread,map<bitset<128>,int,CMPS<128>> & kSolutions)
{
	//logger(vec.to_string());
	vector<bitset<128>> res;
	res.emplace_back(bitset<128>());
	
	
	for(int i = 0 ;i < 549;i++)
	{
		if(vec[i])
		{
			bitset<549> bitVec;
			bitVec[i] = 1;

			map<bitset<128>, int, CMPS<128>> solutions;

			auto status = MidSolutionCounter2(ROUND, cube, bitVec,
				solutions, time, singlethread);
	
			if (status == EXPAND)
			{	
				logger(__func__ + string(" : EXPAND when solving terms in round 256."));
				exit(-1);
			}	
			else if (status == SOLUTION)
			{
				vector<bitset<128>> tmpres;
				for(auto & x : solutions)
					if(x.second % 2)
						tmpres.emplace_back(x.first);

				res = vecMul<128>(res,tmpres);
			}
			else
			{
				logger(__func__ + string(" : NOSOLUTION when solving terms in round 256."));
				exit(-1);
			}
		}
	}

	lock_guard<mutex> guard(stream_mutex);
	for(auto & x : res)
		kSolutions[x]++;
}


int main()
{
	// change cubeindex to set the cube
	vector<int> cubeIndex;
	for (int i = 0; i < 128; i++)
		if( i != 1 && i != 28)
			cubeIndex.push_back(i);


	bitset<128> cube;
	for (int i = 0; i < 128; i++)
		if (count(cubeIndex.begin(), cubeIndex.end(), i))
			cube[i] = 1;
		else
			cube[i] = 0;



	vector < bitset<293>> terms256;

	terms256.clear();
	newAcornForwardExpand(256, cube, terms256, 120, 2);
	if (terms256.size() != 1)
	{
		logger("Error : Size of startcubes are not unique.");
		exit(-1);
	}

	bitset<293> startcube = terms256[0];
	

	int midround = 350;
	ThreadPool thread_pool{};

	int ROUND = 776;
	int step = 100;
	int expandStep = 5;
	float expandTime = 7200;

	float time = 0;

	vector<bitset<549>> terms;
	vector<bitset<549>> layerTerms;

	BackExpandPolynomial(ROUND, step, terms);

	ROUND = ROUND - step;

	logger(__func__ + string("First expand"));

	int singlethread = 2;

	map<bitset<549>, int, CMPS<549>> midSolutions;
	vector<bitset<549> > expandTerms;

	
	while (true)
	{
		if (ROUND > 600)
			time = 600;
		else if (ROUND > 500)
			time = 1800;
		else if (ROUND > 400)
			time = 3600;
		else if (ROUND > 300)
			time = 7200;
		else if (ROUND > 200)
			time = 0;

		layerTerms.clear();
		vector<std::future<void>> futures;
		for (auto& it : terms)
			futures.emplace_back(thread_pool.Submit(SolutionSearcherWorker, ref(it),
				ROUND, ref(startcube), ref(layerTerms), time, singlethread, midround, ref(midSolutions), ref(expandTerms)));
		for (auto& it : futures)
			it.get();

		filterMap<549>(midSolutions);
		filterExpandTerms<549>(expandTerms);

#ifdef _WIN32
#else
		showProcessMemUsage();
		malloc_trim(0);
		showProcessMemUsage();
#endif

		logger(string("layerTermsSize ") + to_string(layerTerms.size()));
		logger(string("expandTermsSize ") + to_string(expandTerms.size()));


		if (layerTerms.size() == 0 && expandTerms.size() == 0)
			break;
		else if (layerTerms.size() == 0)
		{
			logger("Start to handle expandTerms.");
			ROUND = midround;
			midround = 256 + (ROUND - 256) / 2;
			layerTerms = expandTerms;
			expandTerms.clear();
		}

		if (ROUND < 10)
		{
			logger("ROUND is to small.");
			exit(-1);
		}

		terms = layerTerms;


		if (ROUND < FASTBOUND)
			expandStep = 5;

		if (ROUND - expandStep < midround)
			expandStep = ROUND - midround;

		// second Expand
		if (ROUND < CALLBACKBOUND)
			BackExpandCallbackFunc(startcube, terms, ROUND, expandStep, thread_pool, expandTime, singlethread);
		else
			BackExpandDivideFunc(terms, ROUND, expandTime, singlethread);


		// if ROUND == midround,combine expandterms with terms since they are all in midround
		if (ROUND == midround)
		{
			logger("current ROUND : " + to_string(ROUND));
			logger("start combine expandTerms with terms.");
			map<bitset<549>, int, CMPS<549>> termCounter;
			for (auto& x : terms)
				termCounter[x] ++;
			for (auto& x : expandTerms)
				termCounter[x] ++;

			expandTerms.clear();
			terms.clear();
			for (auto& x : termCounter)
				if (x.second % 2 == 1)
					terms.emplace_back(x.first);

			midround = 256 + (ROUND - 256) / 2;
		}

		// back up
		for (auto& x : midSolutions)
			loggerRes(x.first.to_string(), ROUND);
		for (auto& x : expandTerms)
			loggerExpandTerms(x.to_string(), ROUND);

	}
	cout << "Success!" << endl;
	
	// handle midSolutions
	logger("Start to handle midSolutions.");
	logger("midSolutions Size : " + to_string(midSolutions.size()));
		// back up
	for (auto& x : midSolutions)
		loggerState(x.first.to_string(), 256);
	
	


		// extract midcoefficients
	map<bitset<549>, int, CMPS<549>> midCoes;
	for (auto& x : midSolutions)
	{
		bitset<549> midCoe = x.first;
		for (int i = 0; i < 293; i++)
			if(startcube[i] == 1)
				midCoe[i + 256] = 0;
		midCoes[midCoe] ++;
	}
	cout<<midCoes.size()<<endl;

		// compute midcoefficients in K
	map<bitset<128>, int, CMPS<128>> kMidCoes;
	vector<std::future<void>> futures2;
	for (auto& x : midCoes)
		futures2.emplace_back(thread_pool.Submit(SecondSolutionSearcherWorker, ref(x.first), 256, ref(cube), 1200,2,ref(kMidCoes) ));
	for (auto& it : futures2)
		it.get();
	filterMap<128>(kMidCoes);
	cout<<"Compute midcoefficients in K finished."<<endl;
		// compute startcube in K
	map<bitset<128>, int, CMPS<128>> kMidCube;
	bitset<549> MidCube;
	for (int i = 0; i < 293; i++)
		MidCube[i + 256] = startcube[i];
	MidSolutionCounter(256, cube, MidCube, kMidCube, 120, 2);
	filterMap<128>(kMidCube);

		// combine 
	map<bitset<128>, int, CMPS<128>> kSolutions;
	for(auto & x : kMidCoes)
		for (auto& y : kMidCube)
		{
			kSolutions[x.first | y.first] ++;
		}
	filterMap<128>(kSolutions);
	cout<<kSolutions.size()<<endl;


	logger("Solving midSolutions finished.");

	cout << "Output to file." << endl;

	
	string path = string("TERM/") + string("res.txt");
	ofstream os;
	os.open(path, ios::out | ios::app);
	for (auto& it : kSolutions)
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
