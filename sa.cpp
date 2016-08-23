#include <limits>
#include "sa.h"

void SA_Load(const string &s, InputData &D)
{
    string str;
    str = "instances/input_" + s + ".txt";

    ifstream file(str);
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

void SA_Save(const string s, Solution &g, double &v, int t1, int t2)
{
    string str;
    str = "results/output_" + s + ".txt";

    ofstream file(str, ios::app);
    if(!file.is_open())
    {
        cerr << "Output failed!" << endl;
        exit(EXIT_FAILURE);
    }

    stringstream ss;
    time_t t = time(0);
    struct tm * now = localtime( & t );
    ss <<  now->tm_mday << '-'
        << (now->tm_mon + 1) << '-'
        << (now->tm_year + 1900) << " "
        << now->tm_hour << ":"
        << now->tm_min << ":"
        << now->tm_sec << endl;


    vector<bool> y = g.GetY();
    unsigned depots = 0;
    for(unsigned j=0; j<y.size(); j++)
        if(y[j])
            depots++;

    vector<bool> z = g.GetZ();
    unsigned plants = 0;
    for(unsigned k=0; k<z.size(); k++)
        if(z[k])
            plants++;


    file << "INSTANCE " << s << ":\t" << ss.str() << endl
         << "Plants:\t\t" << plants << endl
         << "Depots:\t\t" << depots << endl
         << "Function Value:\t" << v << endl
         << "Time:\t\t" << double(t2-t1)/CLOCKS_PER_SEC << endl << endl << endl;

}

void SA_InitSolution(Solution &s,unsigned J, unsigned K)
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

    s = Solution(y, z);
}


Solution SA_GetRandomNeighbor1(const Solution &s)
{
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

    return Solution(y,z);
}

Solution SA_GetRandomNeighbor2(const Solution &s)
{
    vector<bool> y = s.GetY();
    vector<bool> z = s.GetZ();

    bool b = true;
    while(b)
    {
        z = s.GetZ();

        unsigned r = rand()%z.size();
        z[r] = !z[r];

        for(unsigned k=0; k<z.size(); k++)
            if(z[k])
                b = false;
    }

    return Solution(y,z);
}


double SA_Evaluate(const InputData &data, const Solution &s)
{
    vector<bool> y = s.GetY();
    vector<bool> z = s.GetZ();

    int sum1 = 0;
    for(unsigned j=0; j<data._J; j++)
        sum1 += data._f[j]*y[j];

    if(sum1 == 0)
        return numeric_limits<double>::max();

    int sum2 = 0;
    for(unsigned k=0; k<data._K; k++)
        sum2 += data._g[k]*z[k];

    if(sum2 == 0)
        return numeric_limits<double>::max();

    double sum3 = 0;
    for(unsigned i=0; i<data._I; i++)
    {
        double tmp1 = numeric_limits<double>::max();

        for(unsigned j=0; j<y.size(); j++)
            for(unsigned k=0; k<z.size(); k++)
                if(y[j] && z[k])
                {
                    double tmp2 = data._D[i]*(data._c[i][j] + data._d[j][k]);
                    if(tmp2 < tmp1)
                        tmp1 = tmp2;
                }

        sum3 += tmp1;
    }

    return sum1+sum2+sum3;
}

void SA_ApplyGeometricCooling(double &t)
{
    t *= alpha;
}

