#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
//#include <algorithm>
#include <ctime>
#include "particle.h"

using namespace std;


double vmax = 8.2;
double omega = 0.9;
double fi_p = 3.0;
double fi_g = 0.2;

unsigned parts = 20;
unsigned hsize = 20;
unsigned hoods = ceil(parts/hsize);

vector<Solution> globals;
vector<double> gvalues;


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
    vector<multimap<double, key> > _vm;

    InputData(){}
    InputData(unsigned I, unsigned J, unsigned K,
              vector<unsigned> D, vector<unsigned> f, vector<unsigned> g,
              vector<vector<double > > c, vector<vector<double > > d)
        :_I(I), _J(J), _K(K), _D(D), _f(f), _g(g), _c(c), _d(d)
    {}
};
//-------------------------------------------------------------------------------------------------

void H_Load(const string &s, InputData &D);

void H_Save(const string s, const vector<Solution> &g, const vector<double> &v, int t1, int t2, int t3);

void H_SortMatrix(const vector<vector<double> > &c, const vector<vector<double> > &d,
                    vector<multimap<double, key> > &vm);

bool H_ApplySA(const InputData &data, Solution &global, double &gvalue, double temp, double alpha);

double H_Evaluate1(const InputData &data, const Solution &s);

double H_Evaluate2(const InputData &data, const Solution &s);

double H_Compute(const InputData &data, Particle &p);

double H_GetRandomUniform1(double left, double right);

double H_GetRandomUniform2(double left, double right);


void PS_InitParticles(vector<Particle> &swarm, unsigned n, unsigned J, unsigned K);

void PS_Debug(int n);

void PS_Debug(const Particle &p);


Solution SA_GetRandomNeighbor1(const Solution &s);

Solution SA_GetRandomNeighbor2(const Solution &s);

void SA_ApplyGeometricCooling(double &t, double alpha);

//-------------------HYBRID------------------------------------------------------------------------

//ZAOKRUZIVANJE PRILIKOM UCITAVANJA
void H_Load(const string &s, InputData &D)
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

void H_Save(const string s, const vector<Solution> &g, const vector<double> &v, int t1, int t2, int t3)
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

void H_SortMatrix(const vector<vector<double> > &c, const vector<vector<double> > &d,
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

double H_Evaluate1(const InputData &data, const Solution &s)
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

double H_Evaluate2(const InputData &data, const Solution &s)
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
        for(multimap<double, key>::const_iterator it=data._vm[i].begin();
            it!=data._vm[i].end(); it++)
        {
            k = it->second;
            if(y[k.first] && z[k.second])
                break;
        }

        sum3 += data._D[i]*(data._c[i][k.first] + data._d[k.first][k.second]);
    }

    return sum1+sum2+sum3;
}

double H_Compute(const InputData &data, Particle &p)
{
    Solution s = p.GetCurrentPosition();

    double tmp = H_Evaluate2(data, s);
    if(tmp < p.GetLocal())
    {
        p.SetLocal(tmp);
        p.SetLocalBest(s);
    }

    return tmp;
}

bool H_ApplySA(const InputData &data, Solution &global, double &gvalue, double temp, double alpha)
{
    bool ret = false;

    Solution s = global;
    while(temp > 10000)
    {
        for(unsigned it=0; it<30; it++)
        {
            Solution sp;
            double rnd = H_GetRandomUniform1(0.0, 1.0);

            if(rnd < 0.75)
                sp = SA_GetRandomNeighbor1(s);
            else
                sp = SA_GetRandomNeighbor2(s);

            double sp_value = H_Evaluate2(data, sp);

            if(sp_value < H_Evaluate2(data, global))
            {
                global = sp;
                gvalue = sp_value;
                ret = true;
            }
            if(sp_value < H_Evaluate2(data, s))
                s = sp;
            else
            {
                double delta_f = sp_value - H_Evaluate2(data, s);
                double p = H_GetRandomUniform1(0.0, 1.0);

                if(p > 1/exp(delta_f/temp))
                    s = sp;
            }
        }

        SA_ApplyGeometricCooling(temp, alpha);
    }

    return ret;
}

double H_GetRandomUniform1(double left, double right)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(left, right);

    return distribution(generator);
}

double H_GetRandomUniform2(double left, double right)
{
    double rnd = (double)rand()/RAND_MAX;
    return rnd*(right - left) + left;
}


//------------------PARTICLE-SWARM-------------------------------------------------------------

//ZAOKRUZIVANJE PRILIKOM UCITAVANJA
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
            v_y.push_back(H_GetRandomUniform1(-vmax, vmax));

        vector<double> v_z;
        for(unsigned k=0; k<K; k++)
            v_z.push_back(H_GetRandomUniform1(-vmax, vmax));


        Particle particle(p, Solution(y, z), v_y, v_z);
        swarm.push_back(particle);

        if(p%hsize == 0)
        {
            globals.push_back(Solution(y, z));
            gvalues.push_back(numeric_limits<double>::max());
        }
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


//-----------------SIMULATED-ANNEALING---------------------------------------------------------

//MOZE MALO BOLJI IZBOR OKOLINE(PROMENA SAMO JEDNOG INDEKSA U X[i] A NE OBA)
//MOZDA IPAK BOLJE OVAKO(VECA DIVERSIFIKACIJA)

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

void SA_ApplyGeometricCooling(double &t, double alpha)
{
    t *= alpha;
}


//----------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    srand(time(0));

    if(argc < 2)
    {
        cerr << "Enter instance code!" << endl;
        exit(EXIT_FAILURE);
    }

    InputData data;
    H_Load(argv[1], data);

    int t1 = clock();
    int t2, t3;
    vector<Particle> swarm;

    PS_InitParticles(swarm, parts, data._J, data._K);
    H_SortMatrix(data._c, data._d, data._vm);


    //MOZDA BOLJE OVAKO
    //q < argv[2]
    for(unsigned i=0, q=0; q<30; i++, q++)
    {
        vector<double> v(hoods, numeric_limits<double>::max());
        vector<Solution> g(hoods, Solution());

        if(i%100 == 0)
        {
            cout << "ITERATION " << i << ":" << endl;
            cout << "------------------------------------------------------------" << endl;
        }

        for(unsigned j=0; j<swarm.size(); j++)
        {
            swarm[j].UpdateVelocity();

            double rnd = H_GetRandomUniform1(0.0, 1.0);
            if(rnd < 0.8)
            {
                swarm[j].UpdatePosition1();
//                swarm[j].UpdateVelocity1();
            }
            else
            {
                swarm[j].UpdatePosition2();
//                swarm[j].UpdateVelocity2();
            }

            double tmp = H_Compute(data, swarm[j]);
            if(tmp < v[j/hsize])
            {
                g[j/hsize] = swarm[j].GetCurrentPosition();
                v[j/hsize] = tmp;
            }
        }

        bool b;
        for(unsigned l=0; l<v.size(); l++)
        {
            b = false;
            if(v[l] < gvalues[l])
            {
                cout << "New global best " << l << ": " << v[l] << "(PSO)" << endl;
                g[l].PrintSolution();

                b = true;
                t3 = clock();
            }
            if(H_ApplySA(data, g[l], v[l], 50000, 0.9))
            {
                if(v[l] < gvalues[l])
                {
                    cout << "New global best " << l << ": " << v[l] << "(SA)" << endl;
                    g[l].PrintSolution();

                    b = true;
                    t3 = clock();
                }
            }
            if(b)
            {
                globals[l] = g[l];
                gvalues[l] = v[l];
                q = 0;
            }
        }
        if(b)
        {
            cout << "-----------------------------------------------------" << endl;
        }
    }

    t2 = clock();
    H_Save(argv[1], globals, gvalues, t1, t2, t3);
    cout << (double)(t3-t1)/CLOCKS_PER_SEC << endl;

    return 0;
}

