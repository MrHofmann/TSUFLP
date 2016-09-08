#include "ps.h"
//ZAOKRUZIVANJE PRILIKOM UCITAVANJA
void PS_Load(const string &s, InputData &D)
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

    file.close();
}

void PS_Save(const string s, const vector<Solution> &g, const vector<double> &v, int t1, int t2, int t3)
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

    double tmp = numeric_limits<double>::max();
    unsigned l = 0;
    for(unsigned i=0; i<v.size(); i++)
        if(v[i] < tmp)
        {
            tmp = v[i];
            l = i;
        }


    vector<bool> y = g[l].GetY();
    unsigned depots = 0;
    for(unsigned j=0; j<y.size(); j++)
        if(y[j])
            depots++;

    vector<bool> z = g[l].GetZ();
    unsigned plants = 0;
    for(unsigned k=0; k<z.size(); k++)
        if(z[k])
            plants++;


    file << "INSTANCE " << s << ":\t" << ss.str() << endl
         << "Plants:\t\t" << plants << endl
         << "Depots:\t\t" << depots << endl
         << "Function Value:\t" << gvalues[l] << endl
         << "Time:\t\t" << double(t2-t1)/CLOCKS_PER_SEC << endl
         << "Last Update:\t" << double(t3-t1)/CLOCKS_PER_SEC << endl << endl << endl;

}

void PS_InitParticles(vector<Particle> &swarm, unsigned n, unsigned J, unsigned K)
{
    for(unsigned p=0; p<n; p++)
    {
        //RANDOMIZE Y
        vector<bool> y;
        for(unsigned j=0; j<J; j++)
        {
            unsigned s = rand()%2;
            y.push_back(s);
        }


        //RANDOMIZE Z
        vector<bool> z;
        for(unsigned k=0; k<K; k++)
        {
            unsigned s = rand()%2;
            z.push_back(s);
        }



        vector<double> v_y;
        for(unsigned j=0; j<J; j++)
            v_y.push_back(PS_GetRandomUniform1(-vmax, vmax));

        vector<double> v_z;
        for(unsigned k=0; k<K; k++)
            v_z.push_back(PS_GetRandomUniform1(-vmax, vmax));


        Particle particle(p, Solution(y, z), v_y, v_z);
        swarm.push_back(particle);

        if(p%hsize == 0)
        {
            globals.push_back(Solution(y, z));
            gvalues.push_back(numeric_limits<double>::max());
        }
    }
}

void PS_SortMatrix(const vector<vector<double> > &c, const vector<vector<double> > &d,
                         vector<multimap<double, key> > &vm)
{
    for(unsigned i=0; i<c.size(); i++)
    {
        multimap<double, key> m;
        for(unsigned j=0; j<d.size(); j++)
        {
            for(unsigned k=0; k<d[j].size(); k++)
            {
                double s = c[i][j] + d[j][k];
                m.insert(make_pair(s, key(j, k)));
            }
        }

        vm.push_back(m);
    }
}

double PS_Evaluate(const InputData &data, const Solution &s)
{
    vector<bool> y = s.GetY();
    vector<bool> z = s.GetZ();

    double sum1 = 0;
    for(unsigned j=0; j<data._J; j++)
        sum1 += data._f[j]*y[j];

    if(sum1 == 0)
        return numeric_limits<double>::max();

    double sum2 = 0;
    for(unsigned k=0; k<data._K; k++)
        sum2 += data._g[k]*z[k];

    if(sum2 == 0)
        return numeric_limits<double>::max();


    double sum3 = 0;
    for(unsigned i=0; i<data._vm.size(); i++)
    {
        key k;
        double s = 0;
        for(multimap<double, key>::const_iterator it=data._vm[i].begin();
            it!=data._vm[i].end(); it++)
        {
            k = it->second;
            if(y[k.first] && z[k.second])
            {
                s = it->first;
                break;
            }
        }

        sum3 += data._D[i]*s;
        //sum3 += data._D[i]*(data._c[i][k.first] + data._d[k.first][k.second]);
    }

    return sum1+sum2+sum3;
}

void PS_Compute(const InputData &data, Particle &p)
{
    Solution s(p.GetCurrentPosition().GetY(),p.GetCurrentPosition().GetZ());
    unsigned h = p.GetID()/hsize;

    double tmp = PS_Evaluate(data, s);
    if(tmp < p.GetLocal())
    {
        p.SetLocal(tmp);
        p.SetLocalBest(s);
    }
    if(tmp < gvalues[h])
    {
        globals[h] = s;
        gvalues[h] = tmp;
    }
}

void PS_Debug(int n)
{
    cout << "Here " << n << endl;
}

void PS_Debug(const Particle &p)
{
    vector<bool> y = p.GetCurrentPosition().GetY();
    vector<bool> z = p.GetCurrentPosition().GetZ();

    if(y.empty())
        cerr << "Y prazan" << endl;
    else if(z.empty())
        cerr << "Z prazan" << endl;
}

double PS_GetRandomUniform1(double left, double right)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(left, right);

    return distribution(generator);
}

double PS_GetRandomUniform2(double left, double right)
{
    double rnd = (double)rand()/RAND_MAX;
    return rnd*(right - left) + left;
}
