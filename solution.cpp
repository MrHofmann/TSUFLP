#include "solution.h"


vector<bool> Solution::GetY() const
{
    return _y;
}

vector<bool> Solution::GetZ() const
{
    return _z;
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
    for(unsigned j=0; j<_y.size(); j++)
        cout << _y[j] << " ";
    cout << endl;

    for(unsigned k=0; k<_z.size(); k++)
        cout << _z[k] << " ";
    cout << endl << endl;
}

