#ifndef SA_H
#define SA_H

#include <fstream>
#include <sstream>
#include <iterator>
#include <cmath>
#include "solution.h"

extern unsigned temp;
extern double alpha;
extern Solution global;

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


void SA_Load(const string &s, InputData &D);

void SA_Save(const string s, Solution &g, double &v, int t1, int t2);

void SA_InitSolution(Solution &s, unsigned J, unsigned K);

//MOZE MALO BOLJI IZBOR OKOLINE
Solution SA_GetRandomNeighbor0(const Solution &s);

Solution SA_GetRandomNeighbor1(const Solution &s);

Solution SA_GetRandomNeighbor2(const Solution &s);

double SA_Evaluate(const InputData &data, const Solution &s);

void SA_ApplyGeometricCooling(double &t);

#endif // SA_H
