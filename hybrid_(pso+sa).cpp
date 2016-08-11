#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "particle.h"

using namespace std;


double vmax = 10.0;
double omega = 0.9;
double fi_p = 2.7;
double fi_g = 2.0;

unsigned parts = 60;
unsigned hsize = 60;
unsigned hoods = ceil(parts/hsize);

vector<Solution> globals;
vector<double> gvalues;


unsigned temp = 50000;
double alpha = 0.9;
Solution global;


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


//--------------------------------------------------------------------------------------------

//ZAOKRUZIVANJE PRILIKOM UCITAVANJA
void PS_Load(const string &path, InputData &D)
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

void PS_InitParticles(vector<Particle> &swarm, unsigned n, unsigned I, unsigned J, unsigned K)
{
    for(unsigned p=0; p<n; p++)
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

/*
        for(unsigned j=0; j<J; j++)
            if(j%2)
            {
                y[j] = true;
            }
            else
                y[j] = false;

*/
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
/*
        for(unsigned k=0; k<K; k++)
            if(k==1)
                z[k] = true;
            else
                z[k] = false;
*/
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


        //RANDOMIZE VELOCITY (bolje uniformnom raspodelom)
        vector< map<key, double> > v_x;
        for(unsigned i=0; i<I; i++)
        {
            map<key, double> m;
            map<key, bool>::iterator it = x[i].begin();
            double r = (double)rand()/RAND_MAX;

            m[it->first] = r*2*vmax - vmax;
            v_x.push_back(m);
        }

        vector<double> v_y;
        for(unsigned j=0; j<J; j++)
        {
            double r = (double)rand()/RAND_MAX;
            v_y.push_back(r*2*vmax - vmax);
        }

        vector<double> v_z;
        for(unsigned k=0; k<K; k++)
        {
            double r = (double)rand()/RAND_MAX;
            v_z.push_back(r*2*vmax - vmax);
        }


        Particle particle(p, Solution(x, y, z), v_x, v_y, v_z);
        swarm.push_back(particle);

        if(p%hsize == 0)
        {
            globals.push_back(Solution(x, y, z));
            gvalues.push_back(numeric_limits<double>::max());
        }
    }
}

double PS_Evaluate(const InputData &data, const Solution &s)
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

void PS_Compute(const InputData &data, Particle &p)
{
    Solution s = p.GetCurrentPosition();
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
/*
        cout << "New global best " << h << ": " << tmp << endl;
        vector<bool> y = s.GetY();
        for(unsigned j=0; j<y.size(); j++)
            cout << y[j] << " ";
        cout << endl;

        vector<bool> z = s.GetZ();
        for(unsigned k=0; k<z.size(); k++)
            cout << z[k] << " ";
        cout << endl;
*/    }
}

void PS_Debug(int n)
{
    cout << "Here " << n << endl;
}


//----------------------------------------------------------------------------------------------
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


//----------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    srand(time(0));

    if(argc < 2)
    {
        cerr << "Enter instance path!" << endl;
        exit(EXIT_FAILURE);
    }

    InputData data;
    PS_Load(argv[1], data);

    vector<Particle> swarm;
    PS_InitParticles(swarm, parts, data._I, data._J, data._K);


    for(unsigned i=0; i<1000; i++)
    {
        vector<double> v = gvalues;

        if(i%100 == 0)
        {
            cout << "ITERATION " << i << ":" << endl;
            cout << "------------------------------------------------------------" << endl;
        }

        for(unsigned j=0; j<swarm.size(); j++)
        {
            swarm[j].UpdateVelocity();
            swarm[j].UpdatePosition0();
            PS_Compute(data, swarm[j]);

            swarm[j].UpdateVelocity();
            swarm[j].UpdatePosition1();
            PS_Compute(data, swarm[j]);

            swarm[j].UpdateVelocity();
            swarm[j].UpdatePosition2();
            PS_Compute(data, swarm[j]);
        }

        bool a = false;
        for(unsigned l=0; l<v.size(); l++)
            if(gvalues[l] < v[l])
            {
                cout << "New global best " << l << ": " << gvalues[l] << "(PSO)" << endl;
                vector<bool> y = globals[l].GetY();
                for(unsigned j=0; j<y.size(); j++)
                    cout << y[j] << " ";
                cout << endl;

                vector<bool> z = globals[l].GetZ();
                for(unsigned k=0; k<z.size(); k++)
                    cout << z[k] << " ";
                cout << endl;

//-------------------------------------------------------------------------------
                unsigned iterations = data._I+data._J+data._K;
                double t = temp;
                bool b = false;

                global=globals[l];
                while(t > 0.0001)
                {
                    for(unsigned it=0; it<iterations/2; it++)
                    {
                        Solution sp;
                        if(it%10 == 0)
                            sp = SA_GetRandomNeighbor1(global);
                        else
                            sp = SA_GetRandomNeighbor2(global);


                        double sp_value = SA_Evaluate(data, sp);

                        if(sp_value < SA_Evaluate(data, globals[l]))
                        {
                            globals[l] = sp;
                            gvalues[l] = sp_value;

                            b = true;
                        }
                        if(sp_value < SA_Evaluate(data, global))
                            global = sp;
                        else
                        {
                            //BOLJE UNIFORMNA RASPODELA
                            double delta_f = SA_Evaluate(data, sp) - SA_Evaluate(data, global);
                            double p = (double)rand()/RAND_MAX;

                            if(p > 1/exp(delta_f/t))
                                global = sp;
                        }
                    }

                    SA_ApplyGeometricCooling(t);
                }
                if(b)
                {
                    cout << "New global best " << l << ": " << gvalues[l] << "(SA)" << endl;
                    vector<bool> y = globals[l].GetY();
                    for(unsigned j=0; j<y.size(); j++)
                        cout << y[j] << " ";
                    cout << endl;

                    vector<bool> z = globals[l].GetZ();
                    for(unsigned k=0; k<z.size(); k++)
                        cout << z[k] << " ";
                    cout << endl;
                }
//-------------------------------------------------------------------------------------------------

                a = true;
            }
        if(a)
            cout << "-----------------------------------------------------" << endl;
    }

    return 0;
}
