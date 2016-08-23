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

    vector<double>              _y_velocity;
    vector<double>              _z_velocity;

public:
    Particle()
    {}
    Particle(unsigned id,
             const Solution &curr, vector<double> y_vel, vector<double> z_vel)
        : _id(id), _current(curr), _localp(curr), _y_velocity(y_vel), _z_velocity(z_vel)
    {
        _local = numeric_limits<double>::max();
    }

/*    Particle(const Particle &p)
        :_id(p.GetID()), _current(p.GetCurrentPosition()), _localp(p.GetCurrentPosition())
    {
        _y_velocity = vector<double>(_current.GetY().size(), 0);
        _z_velocity = vector<double>(_current.GetZ().size(), 0);
        _local = numeric_limits<double>::max();
    }

    Particle &operator = (const Particle& p)
    {
        _id = p.GetID();
        _current = p.GetCurrentPosition();
        _localp = p.GetLocalBest();
        _local = p.GetLocal();

        _y_velocity = vector<double>(_current.GetY().size(), 0);
        _z_velocity = vector<double>(_current.GetZ().size(), 0);

        return *this;
    }

    ~Particle(){}
*/

    unsigned GetID() const;
    double GetLocal() const;
    Solution GetCurrentPosition() const;
    Solution GetLocalBest() const;
    Solution GetGlobalBest() const;

    void SetLocal(double local);
    void SetCurrentPosition(const Solution &s);
    void SetGlobalBest(const Solution &s);
    void SetLocalBest(const Solution &s);

    void PrintParticle() const;
    void PrintCurrentPosition() const;
    void PrintGlobalBest() const;
    void PrintVelocity() const;

    void UpdateVelocity();
    void UpdateVelocity1();
    void UpdateVelocity2();

    void UpdatePosition1();
    void UpdatePosition2();
};


#endif // PARTICLE_H
