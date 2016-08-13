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

/*
    PrintVelocity();
    cout << "---------------------------------------------------" << endl << endl;
*/
}


//BOLJE JE UNIFORMNOM RASPODELOM
void Particle::UpdatePosition0()
{
    UpdateX();
}

//BILO JE BOLJE ONAKO
void Particle::UpdatePosition1()
{
//    cout << "PARTICLE(position) " << GetID() << ":" << endl;
    vector<map<key, bool> > x = _current.GetX();
    vector<bool> z = _current.GetZ();

    UpdateY();
    vector<bool> y = _current.GetY();
    vector<unsigned> y_ind;
    for(unsigned j=0; j<y.size(); j++)
        if(y[j])
            y_ind.push_back(j);

    for(unsigned i=0; i<x.size(); i++)
    {
        unsigned j = (x[i].begin()->first).first;
        unsigned k = (x[i].begin()->first).second;

        if(y[j] == 0)
        {
            unsigned r = rand()%y_ind.size();

            x[i].clear();
            x[i][key(y_ind[r], k)] = true;
        }
    }

    _current = Solution(x, y, z);


//    PrintCurrentPosition();
//    cout << "-----------------------------------------" << endl;
//    cout << "-----------------------------------------" << endl << endl;
}

void Particle::UpdatePosition2()
{
//    cout << "PARTICLE(position) " << GetID() << ":" << endl;
    vector<map<key, bool> > x = _current.GetX();
    vector<bool> y = _current.GetY();

    UpdateZ();
    vector<bool> z = _current.GetZ();
    vector<unsigned> z_ind;
    for(unsigned k=0; k<z.size(); k++)
        if(z[k])
            z_ind.push_back(k);

    for(unsigned i=0; i<x.size(); i++)
    {
        unsigned j = (x[i].begin()->first).first;
        unsigned k = (x[i].begin()->first).second;

        if(z[k] == 0)
        {
            unsigned r = rand()%z_ind.size();

            x[i].clear();
            x[i][key(j, z_ind[r])] = true;
        }
    }

    _current = Solution(x, y, z);


//    PrintCurrentPosition();
//    cout << "-----------------------------------------" << endl;
//    cout << "-----------------------------------------" << endl << endl;
}


void Particle::UpdateX()
{
    vector<map<key, bool> > x = _current.GetX();
    vector<bool> y = _current.GetY();
    vector<bool> z = _current.GetZ();

    vector< pair<unsigned, unsigned> > pairs;
    for(unsigned j=0; j<y.size(); j++)
        for(unsigned k=0; k<z.size(); k++)
            if(y[j] && z[k])
                pairs.push_back(make_pair(j, k));


    for(unsigned i=0; i<x.size(); i++)
    {
        vector<key> keys = pairs;
        map<key, double> m = _x_velocity[i];
        x[i].clear();

        do{
            keys = pairs;
            double u = (double)rand()/RAND_MAX;

            unsigned k = 0;
            while(k < keys.size())
            {
                double sig = 1/(1+exp(-m[keys[k]]));

                if(u>=sig)
                {
                    keys[k] = keys.back();
                    keys.pop_back();
                }
                else
                    k++;
            }
        }while(keys.empty());

        unsigned k = rand()%keys.size();
        x[i][keys[k]] = true;
    }

    _current.SetX(x);
}

void Particle::UpdateY()
{
    vector<bool> y = _current.GetY();
    vector<unsigned> y_ind;
    unsigned q = 0;
    while(y_ind.empty())
    {
        q++;
        double u = (double)rand()/RAND_MAX;
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


            if(j%2)
            {
                y[j] = true;
                y_ind.push_back(j);
            }
            else
                y[j] = false;
*/
        }
        if(q%100 == 0)
            cout << q << endl;
    }

    _current.SetY(y);
}

void Particle::UpdateZ()
{
    double u = (double)rand()/RAND_MAX;

    vector<bool> z = _current.GetZ();
    vector<unsigned> z_ind;
    unsigned q = 0;
    while(z_ind.empty() && q < 15)
    {
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

            /*
        //        z[k] = true;
        //        z_ind.push_back(k);


           if(k==1)
           {
               z[k] = true;
               z_ind.push_back(k);
           }
           else
               z[k] = false;
        */
        }
        q++;
    }
    while(z_ind.empty())
        for(unsigned k=0; k<z.size(); k++)
        {
            z[k] = rand()%2;
            if(z[k])
                z_ind.push_back(k);
        }

    _current.SetZ(z);
}
