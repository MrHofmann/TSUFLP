#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iterator>
#include <vector>
#include <map>

using namespace std;


typedef pair<unsigned, unsigned> key;


class InputData
{
public:
    unsigned _I;
    unsigned _J;
    unsigned _K;

    vector<unsigned> _D;
    vector<unsigned> _f;
    vector<unsigned> _g;

    vector<vector<double > > _c;
    vector<vector<double > > _d;

    InputData(){}
    InputData(unsigned I, unsigned J, unsigned K,
              vector<unsigned> D, vector<unsigned> f, vector<unsigned> g,
              vector<vector<double > > c, vector<vector<double > > d)
        :_I(I), _J(J), _K(K), _D(D), _f(f), _g(g), _c(c), _d(d)
    {}
};

class Solution
{
private:
    vector<map<key, bool> > _x;
    vector<bool> _y;
    vector<bool> _z;

public:
    Solution()
    {}
    Solution(vector<map<key, bool> > x, vector<bool> y, vector<bool> z)
        :_x(x), _y(y), _z(z)
    {}

    vector<map<key, bool> > GetX() const;
    vector<bool> GetY() const;
    vector<bool> GetZ() const;

    void PrintSolution() const;
};


unsigned temp = 50000;
double alpha = 0.75;
Solution global;

//------------------------------------------------------------------------------------//

vector<map<key, bool> > Solution::GetX() const
{
    return _x;
}

vector<bool> Solution::GetY() const
{
    return _y;
}

vector<bool> Solution::GetZ() const
{
    return _z;
}


void Solution::PrintSolution() const
{
    for(unsigned i=0; i<_x.size(); i++)
    {
        cout << "Consumer " << i << ": ";
        map<key, bool>::const_iterator it;
        for(it=_x[i].begin(); it!=_x[i].end(); it++)
            cout << "(" << (it->first).first << ", " << (it->first).second << ") = "
                 << it->second << " ";
        cout << endl;
    }
    cout << endl;

    for(unsigned j=0; j<_y.size(); j++)
        cout << _y[j] << " ";
    cout << endl;

    for(unsigned k=0; k<_z.size(); k++)
        cout << _z[k] << " ";
    cout << endl << endl;
}


//------------------------------------------------------------------------------------//

void SA_Load(const string &path, InputData &D)
{
    ifstream file(path);
    if(!file.is_open())
    {
        cerr << "Error open file" << endl;
        exit(EXIT_FAILURE);
    }

    vector<string> data;
    string line;


    getline(file, line);
    {
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(data));

        D._I = stoul(data[0].c_str());
        D._J = stoul(data[1].c_str());
        D._K = stoul(data[2].c_str());
    }


    getline(file, line);
    getline(file, line);


    getline(file, line);
    {
        data.clear();
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(data));
    }
    for(unsigned i=0; i<data.size(); i++)
        D._D.push_back(stoul(data[i]));


    getline(file, line);
    {
        data.clear();
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(data));
    }
    for(unsigned i=0; i<data.size(); i++)
        D._f.push_back(stoul(data[i]));


    getline(file, line);
    {
        data.clear();
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(data));
    }
    for(unsigned i=0; i<data.size(); i++)
        D._g.push_back(stoul(data[i]));


    getline(file, line);
    getline(file, line);


    for(unsigned i=0; i<D._I; i++)
    {
        getline(file, line);

        data.clear();
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(data));

        vector<double> v;
        for(unsigned j=0; j<D._J; j++)
             v.push_back(stod(data[j]));

        D._c.push_back(v);
    }


    getline(file, line);
    getline(file, line);


    for(unsigned j=0; j<D._J; j++)
    {
        getline(file, line);

        data.clear();
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(data));

        vector<double> v;
        for(unsigned k=0; k<D._K; k++)
            v.push_back(stod(data[k]));

        D._d.push_back(v);
    }

}

void SA_InitSolution(Solution &s,unsigned I, unsigned J, unsigned K)
{
    //RANDOMIZE Y
    vector<bool> y;
    vector<unsigned> y_ind;
    while(y_ind.empty())
    {
        y.clear();
        y_ind.clear();
        for(unsigned j=0; j<J; j++)
        {
            unsigned s = rand()%2;
            y.push_back(s);

            if(s == 1)
                y_ind.push_back(j);
        }
    }


    //RANDOMIZE Z
    vector<bool> z;
    vector<unsigned> z_ind;
    while(z_ind.empty())
    {
        z.clear();
        z_ind.clear();
        for(unsigned k=0; k<K; k++)
        {
            unsigned s = rand()%2;
            z.push_back(s);

            if(s == 1)
                z_ind.push_back(k);
        }
    }


    //RANDOMIZE X
    vector< pair<unsigned, unsigned> > pairs;
    for(unsigned j=0; j<y_ind.size(); j++)
        for(unsigned k=0; k<z_ind.size(); k++)
            pairs.push_back(make_pair(y_ind[j], z_ind[k]));

    vector<map<key, bool> > x;
    for(unsigned i=0; i<I; i++)
    {
        map<key, bool> m;
        unsigned k = rand()%pairs.size();
        m[pairs[k]] = true;

        x.push_back(m);
    }

    s = Solution(x, y, z);
}

//MOZE MALO BOLJI IZBOR OKOLINE
Solution SA_GetRandomNeighbor1(const Solution &s)
{
    vector<map<key, bool> > x = s.GetX();
    vector<bool> y = s.GetY();
    vector<bool> z = s.GetZ();

    bool b = true;
    while(b)
    {
        y = s.GetY();

        unsigned r = rand()%y.size();
        y[r] = !y[r];

        for(unsigned j=0; j<y.size(); j++)
            if(y[j])
                b = false;
    }

    b = true;
    while(b)
    {
        z = s.GetZ();

        unsigned r = rand()%z.size();
        z[r] = !z[r];

        for(unsigned k=0; k<z.size(); k++)
            if(z[k])
                b = false;
    }

    vector<key> keys;
    for(unsigned j=0; j<y.size(); j++)
        for(unsigned k=0; k<z.size(); k++)
            if(y[j] && z[k])
                keys.push_back(key(j,k));

    for(unsigned i=0; i<x.size(); i++)
    {
        unsigned j = ((x[i].begin())->first).first;
        unsigned k = ((x[i].begin())->first).second;
        if(!y[j] || !z[k])
        {
            x[i].clear();
            unsigned r3 = rand()%keys.size();
            x[i][keys[r3]] = true;
        }
    }
    return Solution(x,y,z);
}

Solution SA_GetRandomNeighbor2(const Solution &s)
{
    vector<map<key, bool> > x = s.GetX();
    vector<bool> y = s.GetY();
    vector<bool> z = s.GetZ();

    vector<key> keys;
    for(unsigned j=0; j<y.size(); j++)
        for(unsigned k=0; k<z.size(); k++)
            if(y[j] && z[k])
                keys.push_back(key(j,k));

    unsigned r = rand()%x.size();
    unsigned k = rand()%keys.size();

    x[r].clear();
    x[r][keys[k]] = true;

    return Solution(x, y, z);
}

double SA_Evaluate(const InputData &data, const Solution &s)
{
    vector<map<key, bool> > x = s.GetX();
    vector<bool> y = s.GetY();
    vector<bool> z = s.GetZ();

    int sum1 = 0;
    for(unsigned j=0; j<data._J; j++)
        sum1 += data._f[j]*y[j];

    int sum2 = 0;
    for(unsigned k=0; k<data._K; k++)
        sum2 += data._g[k]*z[k];

    double sum3 = 0;
    for(unsigned i=0; i<data._I; i++)
    {
        map<key, bool> m = x[i];
        for(map<key, bool>::const_iterator it=m.begin(); it!=m.end(); it++)
        {
            unsigned j = (it->first).first;
            unsigned k = (it->first).second;

            sum3 += data._D[i]*(it->second)*(data._c[i][j]+data._d[j][k]);
        }
     }

    return sum1+sum2+sum3;
}

void SA_ApplyGeometricCooling(double &t)
{
    t *= alpha;
}

//---------------------------------------------------------------------//


int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        cerr << "Enter instance path!" << endl;
        exit(EXIT_FAILURE);
    }

    srand(time(0));

    InputData data;
    Solution s;

    SA_Load(argv[1], data);
    SA_InitSolution(s, data._I, data._J, data._K);

    unsigned iterations = data._I+data._J+data._K;
    double t = temp;

    global=s;
    while(t > 0.0001)
    {
        for(unsigned i=0; i<iterations/2; i++)
        {
            Solution sp;
            if(i%10 == 0)
                sp = SA_GetRandomNeighbor1(s);
            else
                sp = SA_GetRandomNeighbor2(s);

            if(SA_Evaluate(data,sp) < SA_Evaluate(data, global))
                global = sp;
            if(SA_Evaluate(data, sp) < SA_Evaluate(data, s))
                s = sp;
            else
            {
                //BOLJE UNIFORMNA RASPODELA
                double delta_f = SA_Evaluate(data, sp) - SA_Evaluate(data, s);
                double p = (double)rand()/RAND_MAX;

                if(p > 1/exp(delta_f/t))
                    s = sp;
            }
        }

        SA_ApplyGeometricCooling(t);
    }

    cout << SA_Evaluate(data, global) << endl;

    return 0;
}
