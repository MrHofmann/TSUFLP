#include "particle.h"



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

void Particle::SetCurrentPosition(const Solution &s)
{
    _current = s;
}

void Particle::SetLocalBest(const Solution &s)
{
    _localp = s;
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

void Particle::PrintVelocity() const
{

    for(unsigned j=0; j<_y_velocity.size(); j++)
        cout << _y_velocity[j] << " ";
    cout << endl;

    for(unsigned k=0; k<_z_velocity.size(); k++)
        cout << _z_velocity[k] << " ";
    cout << endl;
}


//BOLJE JE UNIFORMNOM RASPODELOM
void Particle::UpdateVelocity()
{
    for(unsigned j=0; j<_y_velocity.size(); j++)
    {
        double fp = GetRandomUniform1(2.0, 2.5);
        double fg = GetRandomUniform1(2.0, 2.5);

        vector<bool> yc = _current.GetY();
        vector<bool> yl = _localp.GetY();
        vector<bool> yg = globals[GetID()/hsize].GetY();

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
        double fp = GetRandomUniform1(2.0, 2.5);
        double fg = GetRandomUniform1(2.0, 2.5);

        vector<bool> zc = _current.GetZ();
        vector<bool> zl = _localp.GetZ();
        vector<bool> zg = globals[GetID()/hsize].GetZ();

        _z_velocity[k] = omega*_z_velocity[k] +
                fi_p*fp*(zl[k] - zc[k]) +
                fi_g*fg*(zg[k] - zc[k]);

        if(_z_velocity[k] > vmax)
            _z_velocity[k] = vmax;
        else if(_z_velocity[k] < -vmax)
            _z_velocity[k] = -vmax;
    }
}

void Particle::UpdateVelocity1()
{
    for(unsigned j=0; j<_y_velocity.size(); j++)
    {
        double fp = GetRandomUniform1(2.0, 2.5);
        double fg = GetRandomUniform1(2.0, 2.5);

        vector<bool> yc = _current.GetY();
        vector<bool> yl = _localp.GetY();
        vector<bool> yg = globals[GetID()/hsize].GetY();

        _y_velocity[j] = omega*_y_velocity[j] +
                fi_p*fp*(yl[j] - yc[j]) +
                fi_g*fg*(yg[j] - yc[j]);

        if(_y_velocity[j] > vmax)
            _y_velocity[j] = vmax;
        else if(_y_velocity[j] < -vmax)
            _y_velocity[j] = -vmax;
    }
}

void Particle::UpdateVelocity2()
{
    for(unsigned k=0; k<_z_velocity.size(); k++)
    {
        double fp = GetRandomUniform1(2.0, 2.5);
        double fg = GetRandomUniform1(2.0, 2.5);

        vector<bool> zc = _current.GetZ();
        vector<bool> zl = _localp.GetZ();
        vector<bool> zg = globals[GetID()/hsize].GetZ();

        _z_velocity[k] = omega*_z_velocity[k] +
                fi_p*fp*(zl[k] - zc[k]) +
                fi_g*fg*(zg[k] - zc[k]);

        if(_z_velocity[k] > vmax)
            _z_velocity[k] = vmax;
        else if(_z_velocity[k] < -vmax)
            _z_velocity[k] = -vmax;
    }
}


void Particle::UpdatePosition1()
{
    vector<bool> y = _current.GetY();
    double u = GetRandomUniform1(0.0, 1.0);
    for(unsigned j=0; j<y.size(); j++)
    {
        double sig = 1/(1+exp(-_y_velocity[j]));
        if(u<sig)
            y[j] = true;
        else
            y[j] = false;
    }

    _current.SetY(y);
}

void Particle::UpdatePosition2()
{
    double u = GetRandomUniform1(0.0, 1.0);

    vector<bool> z = _current.GetZ();
    for(unsigned k=0; k<z.size(); k++)
    {
        double sig = 1/(1+exp(-_z_velocity[k]));
        if(u<sig)
            z[k] = true;
        else
            z[k] = false;
    }

    _current.SetZ(z);
}
