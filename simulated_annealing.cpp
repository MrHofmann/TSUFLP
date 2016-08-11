#include <iostream>
#include "sa.h"


unsigned temp = 50000;
double alpha = 0.75;
Solution global;


//------------------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        cerr << "Enter instance path!" << endl;
        exit(EXIT_FAILURE);
    }

    srand(time(0));

    InputData data;
    Solution s;

    SA_Load(argv[1], data);
    SA_InitSolution(s, data._I, data._J, data._K);

    unsigned iterations = data._I+data._J+data._K;
    double t = temp;

    global=s;
    while(t > 0.0001)
    {
        for(unsigned i=0; i<iterations/2; i++)
        {
            Solution sp;
            if(i%10 == 0)
                sp = SA_GetRandomNeighbor1(s);
            else
                sp = SA_GetRandomNeighbor2(s);

            if(SA_Evaluate(data,sp) < SA_Evaluate(data, global))
                global = sp;
            if(SA_Evaluate(data, sp) < SA_Evaluate(data, s))
                s = sp;
            else
            {
                //BOLJE UNIFORMNA RASPODELA
                double delta_f = SA_Evaluate(data, sp) - SA_Evaluate(data, s);
                double p = (double)rand()/RAND_MAX;

                if(p > 1/exp(delta_f/t))
                    s = sp;
            }
        }

        SA_ApplyGeometricCooling(t);
    }

    cout << SA_Evaluate(data, global) << endl;

    return 0;
}
