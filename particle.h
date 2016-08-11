#ifndef PARTICLE_H
#define PARTICLE_H

#include <limits>
#include <set>
#include <cmath>
#include "solution.h"

extern double vmax;
extern double omega;
extern double fi_p;
extern double fi_g;

extern unsigned parts;
extern unsigned hsize;
extern unsigned hoods;

extern vector<Solution> globals;
extern vector<double> gvalues;

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

    void UpdateX();
    void UpdateY();
    void UpdateZ();

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
    void PrintVelocity() const;

    void UpdateVelocity();
    void UpdatePosition0();
    void UpdatePosition1();
    void UpdatePosition2();
};


#endif // PARTICLE_H
