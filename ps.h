#ifndef PS_H
#define PS_H

#include <fstream>
#include <sstream>
#include <iterator>
#include "particle.h"

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

void PS_Load(const string &path, InputData &data);

void PS_InitParticles(vector<Particle> &swarm, unsigned n, unsigned I, unsigned J, unsigned K);

void PS_Compute(const InputData &data, Particle &p);

void PS_Debug(int n);

#endif // PS_H
