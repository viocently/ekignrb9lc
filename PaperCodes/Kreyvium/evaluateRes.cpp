

#include <iostream>
#include <string>
#include<bitset>
#include<fstream>
#include<vector>
#include<map>
#include<random>

using namespace std;

template<size_t N>
bool CMP0(const bitset<N>& a, const bitset<N>& b)
{
    for (int i = 0; i < N; i++)
        if (a[i] < b[i]) return true;
        else if (a[i] > b[i]) return false;
    return false; // equal
}

template<size_t N>
struct cmp
{
    bool operator()(const bitset<N>& a, const bitset<N>& b) const
    {
        return CMP0<N>(a, b);
    }
};



template<int n> 
void readTerm(const string& inStr, bitset<n>& _Out)
{
    if (inStr == "1")
    {
        _Out.reset();
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

int main(int argc, char *argv[])
{
    // first read the superpoly into a vector of bitset
    fstream fs;
    string filepath(argv[1]);
    fs.open(filepath);
    
    string termStr;
    vector<bitset<128>> terms;
    map<bitset<128>,string,cmp<128>> termMap;
    while (getline(fs, termStr))
    {
        if (!termStr.empty())
        {
            bitset<128> term;
            vector<int> inList;
            readTerm<128>(termStr, term);
            terms.emplace_back(term);
            termMap[term] = termStr;
        }
    }
    
    fs.close();

    // analyze
    bitset<128> inBits; // involved key bits
    bitset<128> maxTerm; // term with maximum degree
    vector<vector<bitset<128>>> inMap(128); // classify terms according to involved bits

    for (auto& x : terms)
    {
        inBits |= x;
        if (maxTerm.count() < x.count())
            maxTerm = x;
        for (int i = 0; i < 128; i++)
            if (x[i])
                inMap[i].emplace_back(x);
    }

    // output
    cout << "Involved key bits: " << inBits.count() << endl;
    cout << "Degree: " << maxTerm.count() << endl;
    for(int i = 0; i < 128;i++)
        if (inMap[i].size() < 2)
        {
            cout << "Bit " << i << " involves terms: " << inMap[i].size() << endl;
            if (inMap[i].size() == 1)
            {
                        cout << termMap[inMap[i][0]] << endl;
            }
        }

    // test balancedness
    default_random_engine e;
    e.seed(time(0));
    uniform_int_distribution<unsigned long long> u(1, UINT64_MAX);
    int cnt_1 = 0;
    for (unsigned int i = 0; i < (1 << 15); i++)
    {
        bitset<64> kLow(u(e));
        bitset<64> kHigh(u(e));
        bitset<128> k(kHigh.to_string() + kLow.to_string());

        int tmp = 0;
        for (auto& x : terms)
            if ((x & k) == x)
                tmp++;
        if (tmp % 2)
            cnt_1++;
    }

    cout << "Number of output 1: " << cnt_1 << endl;
    cout << "Balancedness : " << ((double)cnt_1) / (1 << 15) << endl;



}

