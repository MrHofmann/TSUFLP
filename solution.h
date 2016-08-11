#ifndef SOLUTION_H
#define SOLUTION_H

#include <iostream>
#include <vector>
#include <map>

using namespace std;


typedef pair<unsigned, unsigned> key;

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

    void SetX(const vector<map<key, bool> > x);
    void SetY(const vector<bool> y);
    void SetZ(const vector<bool> z);

    void PrintSolution() const;
};

#endif // SOLUTION_H
