#include "sa.h"

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
