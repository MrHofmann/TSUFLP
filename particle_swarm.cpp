#include <iostream>
#include <ctime>
#include "ps.h"


double vmax = 8.2;
double omega = 0.99;
double fi_p = 3.0;
double fi_g = 0.2;

unsigned parts = 20;
unsigned hsize = 20;
unsigned hoods = ceil(parts/hsize);

vector<Solution> globals;
vector<double> gvalues;


int main(int argc, char **argv)
{
    srand(time(0));

    if(argc < 2)
    {
        cerr << "Enter instance path!" << endl;
        exit(EXIT_FAILURE);
    }

    InputData data;
    PS_Load(argv[1], data);
    PS_SortMatrix(data._c, data._d, data._vm);

    int t1 = clock();
    int t2, t3;
    vector<Particle> swarm;
    PS_InitParticles(swarm, parts, data._J, data._K);

    unsigned q = 0;
    for(unsigned i=0; q<100; i++)
    {
        vector<double> v = gvalues;

        for(unsigned j=0; j<swarm.size(); j++)
        {
            swarm[j].UpdateVelocity();

            double rnd = PS_GetRandomUniform1(0.0, 1.0);
            if(rnd < 0.8)
            {
                swarm[j].UpdatePosition1();
//                swarm[j].UpdateVelocity1();
            }
            else
            {
                swarm[j].UpdatePosition2();
//                swarm[j].UpdateVelocity2();
            }

            PS_Compute(data, swarm[j]);
        }

        bool b = false;
        for(unsigned l=0; l<v.size(); l++)
            if(gvalues[l] < v[l])
            {
                cout << "New global best " << l << ": " << gvalues[l] << endl;
                vector<bool> y = globals[l].GetY();
                for(unsigned j=0; j<y.size(); j++)
                    cout << y[j] << " ";
                cout << endl;

                vector<bool> z = globals[l].GetZ();
                for(unsigned k=0; k<z.size(); k++)
                    cout << z[k] << " ";
                cout << endl;

                b = true;
                q = 0;
                t3 = clock();
            }
        if(b)
        {
            cout << "-----------------------------------------------------" << endl;
        }
        q++;
    }

    t2 = clock();
    PS_Save(argv[1], globals, gvalues, t1, t2);
    cout << (double)(t3-t1)/CLOCKS_PER_SEC;

    return 0;
}
