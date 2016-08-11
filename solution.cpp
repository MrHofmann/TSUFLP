#include "solution.h"


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


void Solution::SetX(const vector<map<key, bool> > x)
{
    _x = x;
}

void Solution::SetY(const vector<bool> y)
{
    _y = y;
}

void Solution::SetZ(const vector<bool> z)
{
    _z = z;
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

