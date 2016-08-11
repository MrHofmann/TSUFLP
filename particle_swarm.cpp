#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <tuple>
#include <ctime>
#include <numeric>
#include <limits>
#include <cmath>
#include <set>
#include <map>
#include <iterator>

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

class Particle
{
private:
    unsigned                    _id;
    double                      _local;

    Solution                    _current;
    Solution                    _localp;
    Solution                    _global;

    vector< map<key, double> >  _x_velocity;
    vector<double>              _y_velocity;
    vector<double>              _z_velocity;

public:
    Particle()
    {}
    Particle(unsigned id,
             const Solution &curr, vector< map<key, double> > x_vel,
             vector<double> y_vel, vector<double> z_vel)
        : _id(id), _current(curr), _localp(curr),
          _x_velocity(x_vel), _y_velocity(y_vel), _z_velocity(z_vel)
    {
        _local = numeric_limits<double>::max();
    }

    unsigned GetID() const;
    double GetLocal() const;
    Solution GetCurrentPosition() const;
    Solution GetLocalBest() const;
    Solution GetGlobalBest() const;

    void SetLocal(double local);
    void SetGlobalBest(const Solution &s);
    void SetLocalBest(const Solution &s);

    void PrintParticle() const;
    void PrintCurrentPosition() const;
    void PrintGlobalBest() const;
//    void PrintLocalBest() const;
    void PrintVelocity() const;

    double Evaluate(unsigned I, unsigned J, unsigned K,
                  vector<unsigned> &D, vector<unsigned> &f, vector<unsigned> &g,
                  vector< vector<double> > &c, vector< vector<double> > &d);
    void UpdateVelocity();
    void UpdatePosition();
};

double vmax = 10.0;
double omega = 0.9;
double fi_p = 2.7;
double fi_g = 2.0;

unsigned parts = 60;
unsigned hsize = 60;
unsigned hoods = ceil(parts/hsize);

vector<Solution> globals;
vector<double> gvalues;

int jota = 0;
//------------------------------------------------------------------------------------//


void PS_Load(const string &path, InputData &data);

void PS_InitParticles(vector<Particle> &swarm, unsigned n, unsigned I, unsigned J, unsigned K);

void PS_ComputeGlobal(vector<Particle> &swarm, const InputData &data);

void PS_Debug(int n);

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

unsigned Particle::GetID() const
{
    return _id;
}

double Particle::GetLocal() const
{
    return _local;
}

Solution Particle::GetCurrentPosition() const
{
    return _current;
}

Solution Particle::GetLocalBest() const
{
    return _localp;
}

Solution Particle::GetGlobalBest() const
{
    return _global;
}


void Particle::SetLocalBest(const Solution &s)
{
    _localp = s;
}

void Particle::SetGlobalBest(const Solution &s)
{
    _global = s;
}

void Particle::SetLocal(double local)
{
    _local = local;
}


void Particle::PrintParticle() const
{
    cout << "Particle " << _id << ":" << endl << endl;

    PrintCurrentPosition();
    PrintVelocity();

    cout << "----------------------------------------------------------------" << endl;

}

void Particle::PrintCurrentPosition() const
{
    _current.PrintSolution();
}

void Particle::PrintGlobalBest() const
{
    _global.PrintSolution();
}

void Particle::PrintVelocity() const
{
    for(unsigned i=0; i<_x_velocity.size(); i++)
    {
        cout << "Consumer " << i << ": ";
        map<key, double>::const_iterator it;
        for(it=_x_velocity[i].begin(); it!=_x_velocity[i].end(); it++)
            cout << "(" << (it->first).first << ", " << (it->first).second << ") = "
                 << it->second << " ";
        cout << endl;
    }
    cout << endl;

    for(unsigned j=0; j<_y_velocity.size(); j++)
        cout << _y_velocity[j] << " ";
    cout << endl;
    for(unsigned k=0; k<_z_velocity.size(); k++)
        cout << _z_velocity[k] << " ";
    cout << endl;
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


//BOLJE JE UNIFORMNOM RASPODELOM, UPDATE ZA X NIJE GOTOV
void Particle::UpdateVelocity()
{
    double fp = (double)rand()/RAND_MAX;
    double fg = (double)rand()/RAND_MAX;
    fp = fp*(1.0 - 0.5) + 0.5;
    fg = fg*(1.0 - 0.5) + 0.5;

//    cout << "PARTICLE(velocity) " << GetID() << ":" << endl << endl;
    for(unsigned i=0; i<_x_velocity.size(); i++)
    {
        map<key, bool> xc = GetCurrentPosition().GetX()[i];
        map<key, bool> xl = GetLocalBest().GetX()[i];
        map<key, bool> xg = globals[GetID()/hsize].GetX()[i];
        map<key, double> vi = _x_velocity[i];

        set<key> keys;

        keys.insert(xc.begin()->first);
        keys.insert(xl.begin()->first);
        keys.insert(xg.begin()->first);
        for(map<key, double>::const_iterator it=vi.begin(); it!=vi.end(); it++)
            keys.insert(it->first);


        for(set<key>::const_iterator it=keys.begin(); it!=keys.end(); it++)
        {
            key k = *it;
            double v = omega*vi[k] + fi_p*fp*(xl[k] - xc[k])+ fi_g*fg*(xg[k] - xc[k]);

            v = (v > vmax)? vmax : v;
            v = (v < -vmax)? -vmax : v;

            _x_velocity[i][k] = v;
        }
    }

    for(unsigned j=0; j<_y_velocity.size(); j++)
    {
        vector<bool> yc = _current.GetY();
        vector<bool> yl = _localp.GetY();
        vector<bool> yg = _global.GetY();

        _y_velocity[j] = omega*_y_velocity[j] +
                fi_p*fp*(yl[j] - yc[j]) +
                fi_g*fg*(yg[j] - yc[j]);

        if(_y_velocity[j] > vmax)
            _y_velocity[j] = vmax;
        else if(_y_velocity[j] < -vmax)
            _y_velocity[j] = -vmax;
    }


    for(unsigned k=0; k<_z_velocity.size(); k++)
    {
        vector<bool> zc = _current.GetZ();
        vector<bool> zl = _localp.GetZ();
        vector<bool> zg = _global.GetZ();

        _z_velocity[k] = omega*_z_velocity[k] +
                fi_p*fp*(zl[k] - zc[k]) +
                fi_g*fg*(zg[k] - zc[k]);

        if(_z_velocity[k] > vmax)
            _z_velocity[k] = vmax;
        else if(_z_velocity[k] < -vmax)
            _z_velocity[k] = -vmax;
    }

/*
    PrintVelocity();
    cout << "---------------------------------------------------" << endl << endl;
*/
}


//BOLJE JE UNIFORMNOM RASPODELOM, UPOREDI SA LOKALNIM RESENJEM OVDE, UPDATE X TREBA DA SE POPRAVI
void Particle::UpdatePosition()
{
    double u = (double)rand()/RAND_MAX;

//    cout << "PARTICLE(position) " << GetID() << ":" << endl;
    //UPDATE Y

    int t1 = clock();
    vector<bool> y = _current.GetY();
    vector<unsigned> y_ind;
    while(y_ind.empty())
    {
        u = (double)rand()/RAND_MAX;
        for(unsigned j=0; j<y.size(); j++)
        {
            double sig = 1/(1+exp(-_y_velocity[j]));
            if(u<sig)
            {
                y[j] = true;
                y_ind.push_back(j);
            }
            else
                y[j] = false;

/*
            y[j] = true;
            y_ind.push_back(j);
*/

/*
            if(j%2)
            {
                y[j] = true;
                y_ind.push_back(j);
            }
            else
                y[j] = false;
*/

        }
    }

    //UPDATE Z
    vector<bool> z = _current.GetZ();
    vector<unsigned> z_ind;
    for(unsigned k=0; k<z.size(); k++)
    {
        double sig = 1/(1+exp(-_z_velocity[k]));
        if(u<sig)
        {
            z[k] = true;
            z_ind.push_back(k);
        }
        else
            z[k] = false;


//        z[k] = true;
//        z_ind.push_back(k);

/*
       if(k==1)
       {
           z[k] = true;
           z_ind.push_back(k);
       }
       else
           z[k] = false;
*/
    }

    if(z_ind.empty())
    {
        unsigned k = rand()%z.size();
        z[k] = true;
        z_ind.push_back(k);
    }


//    PrintParticle();

    //UPDATE X, mozda bi trebalo da inicijalizujem celu matricu xvelocity
    vector< pair<unsigned, unsigned> > pairs;
    for(unsigned j=0; j<y_ind.size(); j++)
        for(unsigned k=0; k<z_ind.size(); k++)
            pairs.push_back(make_pair(y_ind[j], z_ind[k]));


    vector<map<key, bool> > x = _current.GetX();
    for(unsigned i=0; i<x.size(); i++)
    {
        vector< pair<unsigned, unsigned> > ppairs = pairs;
        map<key, double> m = _x_velocity[i];
        x[i].clear();

        //PRVI NACIN, MOGUCA BESKONACNA PETLJA(nije se skoro javljala)
/*        bool b = true;
        while(b)
        {
            u = (double)rand()/RAND_MAX;

            unsigned j = 0;
            while(j < ppairs.size())
            {
                double sig = 1/(1+exp(-m[ppairs[j]]));

                if(ppairs.size() == 1)
                {
                    b = false;
                    break;
                }
                else if(u>=sig)
                {
                    ppairs[j] = ppairs.back();
                    ppairs.pop_back();
                }
                else
                    j++;
            }
        }

        _x_current[i][ppairs[0]] = true;
*/

        //DRUGI NACIN
        do{
            ppairs = pairs;
            u = (double)rand()/RAND_MAX;

            unsigned j = 0;
            while(j < ppairs.size())
            {
                double sig = 1/(1+exp(-m[ppairs[j]]));

                if(u>=sig)
                {
                    ppairs[j] = ppairs.back();
                    ppairs.pop_back();
                }
                else
                    j++;
            }
        }while(ppairs.empty());

        unsigned k = rand()%ppairs.size();
        x[i][ppairs[k]] = true;
    }

    _current = Solution(x, y, z);
    int t2 = clock();
    if(jota%30 == 0)
        cout << t2 - t1 << endl;

//    PrintCurrentPosition();
//    cout << "-----------------------------------------" << endl;
//    cout << "-----------------------------------------" << endl << endl;

}

//------------------------------------------------------------------------------------

//NE VALJA ISTRINGSTREAM ISS, ZAOKRUZIVANJE PRILIKOM UCITAVANJA
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

//BOLJA IMPLEMENTACIJA AZURIRANJA GLOBALNOG RESENJA
void PS_ComputeGlobal(vector<Particle> &swarm, const InputData &data)
{
    /*IMPLEMENTACIJA SA NUMBER OF NAEIGHBORHOODS
    //OVO MORA GLOBALNO I INICIJALIZOVANO
    vector<Solution> glob;

    int k = 4;
    int n = swarm.size();

    unsigned q1 = (n/k +1);
    unsigned q2 = (n%k)*(n/k + 1);
    unsigned r = n%k;
    for(unsigned h=0; h<r; h++)
    {
        for(unsigned i=h*q1; i<(h+1)*q1; i++)
        {
            Solution s = swarm[i].GetCurrentPosition();
            double tmp = PS_Evaluate(data, s);

            if(tmp < swarm[i].GetLocal())
            {
                swarm[i].SetLocal(tmp);
                swarm[i].SetLocalBest(swarm[i].GetCurrentPosition());
            }

            if(tmp < PS_Evaluate(data, glob[h]))
                glob[h] = swarm[i].GetCurrentPosition();

        }
    }
    for(unsigned h=0; h<k-r; h++)
    {
        for(unsigned i=q2+h*(n/k); i<q2+(h+1)*(n/k); i++)
        {
            for(unsigned i=h*q1; i<(h+1)*q1; i++)
            {
                Solution s = swarm[i].GetCurrentPosition();
                double tmp = PS_Evaluate(data, s);

                if(tmp < swarm[i].GetLocal())
                {
                    swarm[i].SetLocal(tmp);
                    swarm[i].SetLocalBest(swarm[i].GetCurrentPosition());
                }

                if(tmp < PS_Evaluate(data, glob[h]))
                    glob[h] = swarm[i].GetCurrentPosition();

            }
        }
    }
*/

    //IMPLEMENTACIJA SA HOOD SIZE
    unsigned h=0;
    for(unsigned i=0; i<parts; i++)
    {
        Solution s = swarm[i].GetCurrentPosition();
        double tmp = PS_Evaluate(data, s);

        if(tmp < swarm[i].GetLocal())
        {
            swarm[i].SetLocal(tmp);
            swarm[i].SetLocalBest(swarm[i].GetCurrentPosition());
        }
        if(tmp < gvalues[h])
        {
            globals[h] = swarm[i].GetCurrentPosition();
            gvalues[h] = tmp;

            cout << "New global best " << h << ": " << tmp << endl;
            vector<bool> y = s.GetY();
            for(unsigned j=0; j<y.size(); j++)
                cout << y[j] << " ";
            cout << endl;

            vector<bool> z = s.GetZ();
            for(unsigned k=0; k<z.size(); k++)
                cout << z[k] << " ";
            cout << endl;
        }

        if((i+1)%hsize)
            h++;
    }
}

void PS_Debug(int n)
{
    cout << "Here " << n << endl;
}

//------------------------------------------------------------------------------------

int main(int argc, char **argv)
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
        if(i%100 == 0)
        {
            cout << "ITERATION " << i << ":" << endl;
            cout << "------------------------------------------------------------" << endl;
        }


        PS_ComputeGlobal(swarm, data);


        for(unsigned j=0; j<swarm.size(); j++)
        {
            swarm[j].UpdateVelocity();
            swarm[j].UpdatePosition();

            jota++;
        }

    }

    return 0;
}
