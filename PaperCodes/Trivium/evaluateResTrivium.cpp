#include <iostream>
#include <string>
#include<bitset>
#include<fstream>
#include<vector>
#include<map>
#include<random>
#include<filesystem>
#include<set>
#include<malloc.h>
#include<ctime>
#include<unistd.h>
#include"thread_pool.h"
#include"dynamic_bitset.hpp"

#define KEYLEN 80
using namespace std;
using namespace boost;
using namespace thread_pool;

mutex stream_mutex;


template<size_t N>
bool CMP0(const bitset<N>& a, const bitset<N>& b)
{
    for (int i = 0; i < N; i++)
        if (a[i] < b[i]) return true;
        else if (a[i] > b[i]) return false;
    return false; // equal
}

template<size_t N>
bool CMP1(const vector<bitset<N>>& a, const vector<bitset<N> >& b)
{
    if (a.size() < b.size())
        return true;
    else if (a.size() > b.size())
        return false;

    auto cmp = CMP0<N>;
    for (int i = 0; i < a.size(); i++)
        if (cmp(a[i], b[i]) == true)
            return true;
        else if (cmp(b[i], a[i]) == true)
            return false;
        else
            continue;
    return false;
}

template<size_t N>
bool CMP2(const vector<vector<bitset<N>>>& a, const vector<vector<bitset<N> > >& b)
{
    if (a.size() < b.size())
        return true;
    else if (a.size() > b.size())
        return false;

    auto cmp = CMP1<N>;
    for (int i = 0; i < a.size(); i++)
        if (cmp(a[i], b[i]) == true)
            return true;
        else if (cmp(b[i], a[i]) == true)
            return false;
        else
            continue;
    return false;
}



template<size_t N>
struct cmp
{
    bool operator()(const bitset<N>& a, const bitset<N>& b) const
    {
        return CMP0<N>(a, b);
    }
};

template<size_t N>
struct vcmp
{
    bool operator()(const vector<bitset<N>>& a, const vector<bitset<N> >& b) const
    {
        return CMP1<N>(a, b);
    }
};

template<size_t N>
struct vvcmp
{
    bool operator()(const vector<vector<bitset<N>>>& a, const vector<vector<bitset<N> > >& b) const
    {
        return CMP2<N>(a, b);
    }
};


template<int n>
vector<bitset<n>> vecMul(const vector<bitset<n>>& _left, const vector<bitset<n>>& _right)
{
    map<bitset<n>, int, cmp<n>> resMap;
    vector<bitset<n>> res;

    for (auto& x : _left)
        for (auto& y : _right)
            resMap[x | y] ++;

    for (auto& x : resMap)
        if (x.second % 2)
            res.emplace_back(x.first);

    return res;
}

template<int n>
vector<bitset<n>> vecAdd(vector<bitset<n>>& _left, vector<bitset<n>>& _right)
{
    map<bitset<n>, int, cmp<n>> resMap;
    vector<bitset<n>> res;
    for (auto& x : _left)
        resMap[x]++;
    for (auto& y : _right)
        resMap[y] ++;

    for (auto& x : resMap)
        if (x.second % 2)
            res.emplace_back(x.first);

    return res;
}

void getJustCurrentFile(const string& path, vector<string>& files)
{
    for (const auto& entry : filesystem::directory_iterator(path))
        files.emplace_back(entry.path().string());

}



template<int n>
void readTerm(const string& inStr, bitset<n>& _Out)
{
    _Out.reset();

    if (inStr == "1")
    {
        return;
    }

    int iScan = 0;
    while (iScan != inStr.length())
    {
        string bitStr;
        bitStr += inStr[iScan++];
        while (inStr[iScan] != 'k' && iScan != inStr.length())
            bitStr += inStr[iScan++];

        _Out[stoi(bitStr.substr(1))] = 1;
    }
}

template<int n>
void readTerms(const string& inStr, vector<bitset<n>>& _Out)
{
    int iScan = 0;
    string termStr;
    bitset<n> term;
    while (iScan != inStr.length())
    {
        if (inStr[iScan] != '+')
            termStr += inStr[iScan];
        else
        {
            readTerm<n>(termStr, term);
            _Out.emplace_back(term);
            termStr.clear();
        }
        iScan++;
    }

    readTerm<n>(termStr, term);
    _Out.emplace_back(term);
}

template<int n>
void readLine(const string& inStr, vector<vector<bitset<n>>>& _Out)
{
    int iScan = 0;
    string termsStr;
    while (iScan != inStr.length())
    {
        if (inStr[iScan] == '(' || inStr[iScan] == '*')
            ;
        else if (inStr[iScan] == ')')
        {
            vector<bitset<n>> terms;
            readTerms<n>(termsStr, terms);
            _Out.emplace_back(terms);
            termsStr.clear();
        }
        else
            termsStr += inStr[iScan];

        iScan++;
    }
}


void readFileThread(const string filename, map<vector<vector<bitset<KEYLEN>>>, int, vvcmp<KEYLEN>>& lineMap)
{

    fstream fs;
    fs.open(filename, ios::in);
    string line;
    map<vector<vector<bitset<KEYLEN>>>, int, vvcmp<KEYLEN>> threadLineMap;
    while (getline(fs, line))
    {
        if (line.empty())
            continue;
        vector<vector<bitset<KEYLEN>>> lineExp;
        readLine<KEYLEN>(line, lineExp);
        threadLineMap[lineExp]++;

    }
    fs.close();

    lock_guard<mutex> guard(stream_mutex);
    cout << "Read " << filename << " completed." << endl;
    for (auto& x : threadLineMap)
        if (x.second % 2)
            lineMap[x.first]++;


}


void clearMem()
{
    while (1)
    {
        sleep(120);
        malloc_trim(0);
    }
}



template<int n>
string termToStr(bitset<n>& term)
{
    if (term.count() == 0)
        return "1";

    string res;

    for (int i = 0; i < n; i++)
        if (term[i])
            res += 'k' + to_string(i);

    return res;
}


void testBalancedNessThread(const vector<vector<bitset<KEYLEN>>> & polys,const vector<dynamic_bitset<> > & lineBitExps, 
    const bitset<KEYLEN> k, int& cnt_1)
{
    int polyNum = polys.size();

    dynamic_bitset<> lineBitVal(polyNum);


    for (int i = 0; i < polyNum;i++)
    {
        int tmpcnt = 0;
        for (auto& y : polys[i])
        {
            if ((y & k) == y)
                tmpcnt++;
        }

        if (tmpcnt % 2)
            lineBitVal[i] = 1;
    }

    int tmpcnt2 = 0;
    for (auto& x : lineBitExps)
        if ((x & lineBitVal) == x)
            tmpcnt2++;

    if (tmpcnt2 % 2)
        cnt_1++;

}




int main(int argc, char* argv[])
{
    vector<string> files;
    for (int i = 1; i < argc; i++)
    {
        string dirpath(argv[i]);
        getJustCurrentFile(dirpath, files);
    }






    map<vector<vector<bitset<KEYLEN>>>, int, vvcmp<KEYLEN>> lineMap;

    ThreadPool thread_pool{};

    thread t(clearMem);

    vector<std::future<void>> futures;

    for (auto& file : files)
        if (file.find("term") != string::npos)
            futures.emplace_back(thread_pool.Submit(readFileThread, file, ref(lineMap)));

    for (auto& x : futures)
        x.get();

    t.detach();
    malloc_trim(0);

    map<vector<bitset<KEYLEN>>, int, vcmp<KEYLEN>> polyMap;
    for (auto& x : lineMap)
        if (x.second % 2)
        {
            for (auto& y : x.first)
                polyMap[y]++;
        }

    map<vector<bitset<KEYLEN>>, int, vcmp<KEYLEN>> polyIndexMap;
    int index = 0;
    for (auto& x : polyMap)
            polyIndexMap[x.first] = index++;

    int polyNum = polyIndexMap.size();
    cout << "Number of polys : " << polyNum << endl;

    vector<dynamic_bitset<> > lineBitExps;

    for(auto & x : lineMap)
        if (x.second % 2)
        {
            dynamic_bitset<> lineBitExp(polyNum);
            for (auto& y : x.first)
                lineBitExp[polyIndexMap[y]] = 1;
            lineBitExps.emplace_back(lineBitExp);
        }

    vector<vector<bitset<KEYLEN>>> polys(polyNum);
    for (auto& x : polyIndexMap)
        polys[x.second] = x.first;

    // test balancedness
    cout << "Start testing balancedness." << endl;


    default_random_engine e;
    e.seed(time(0));
    uniform_int_distribution<unsigned long long> u(1, UINT64_MAX);
    int uN = KEYLEN / 64 + 1;
    int uL = KEYLEN % 64;

    int threadnum = thread_pool.num_threads();
    vector<int> cnt_1s(threadnum, 0);
    futures.clear();

    for (unsigned int i = 0; i < (1 << 15); i++)
    {
        // generate random key values
        string kStr;
        for (int j = 0; j < uN - 1; j++)
        {
            bitset<64> kj(u(e));
            kStr += kj.to_string();
        }

        bitset<64> kl(u(e));
        kStr += kl.to_string().substr(0, uL);

        bitset<KEYLEN> k(kStr);


        // count the value of superpoly
        futures.emplace_back(thread_pool.Submit(testBalancedNessThread, ref(polys), ref(lineBitExps),k, ref(cnt_1s[i % threadnum])));
    }

    for (auto& x : futures)
        x.get();


    int cnt_1 = 0;
    for (auto& x : cnt_1s)
        cnt_1 += x;



    cout << "Number of output 1: " << cnt_1 << endl;
    cout << "Balancedness : " << ((double)cnt_1) / (1 << 15) << endl;












}

