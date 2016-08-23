#ifndef SOLUTION_H
#define SOLUTION_H

#include <iostream>
#include <vector>

using namespace std;


typedef pair<unsigned, unsigned> key;

class Solution
{
private:
    vector<bool> _y;
    vector<bool> _z;

public:
    Solution()
    {}
    Solution(vector<bool> y, vector<bool> z)
        :_y(y), _z(z)
    {}

/*
    Solution(const Solution &s)
        :_y(s.GetY()), _z(s.GetZ())
    {

    }

    Solution &operator = (const Solution &s)
    {
        _y = s.GetY();
        _z = s.GetZ();

        return *this;
    }

    ~Solution()
    {

    }
*/

    vector<bool> GetY() const;
    vector<bool> GetZ() const;

    void SetY(const vector<bool> y);
    void SetZ(const vector<bool> z);

    void PrintSolution() const;
};

#endif // SOLUTION_H
