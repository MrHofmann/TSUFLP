#include <iostream>
#include "sa.h"


unsigned temp = 5000;
double alpha = 0.9;
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
    SA_Load(argv[1], data);

    int t1 = clock();
    int t2, t3;
    SA_InitSolution(global, data._J, data._K);
    SA_SortMatrix(data._c, data._d, data._vm);

    unsigned iterations = data._J + data._K;
    double t = temp;
    for(unsigned k=0; t>0.1 && k<6; k++)
    {
        double tmp = SA_Evaluate(data, global);
        Solution s;
        SA_InitSolution(s, data._J, data._K);

        for(unsigned i=0; i<iterations*10; i++)
        {
            Solution sp;
            double rnd = SA_GetRandomUniform1(0.0, 1.0);

            if(rnd < 0.75)
                sp = SA_GetRandomNeighbor1(s);
            else
            {
                sp = SA_GetRandomNeighbor1(s);
                SA_Compute(data, t, s, sp, global);

                sp = SA_GetRandomNeighbor2(s);
            }

            SA_Compute(data, t, s, sp, global);
        }

        SA_ApplyGeometricCooling(t);

        if(SA_Evaluate(data, global) < tmp)
        {
            t3 = clock();
            k=0;
        }
    }

    t2 = clock();
    double tmp = SA_Evaluate(data, global);
    SA_Save(argv[1], global, tmp, t1, t2);
    cout << (double)(t3-t1)/CLOCKS_PER_SEC;

    return 0;
}
