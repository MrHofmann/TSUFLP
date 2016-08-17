#include <iostream>
#include <ctime>
#include "ps.h"


double vmax = 10.0;
double omega = 0.9;
double fi_p = 2.7;
double fi_g = 1.0;

unsigned parts = 60;
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

    vector<Particle> swarm;
    PS_InitParticles(swarm, parts, data._I, data._J, data._K);


    int t1 = clock();
    int t2, t3;
    unsigned q = 0;
    for(unsigned i=0; i<1000 && q<100; i++)
    {
        vector<double> v = gvalues;

        if(i%100 == 0)
        {
            cout << "ITERATION " << i << ":" << endl;
            cout << "------------------------------------------------------------" << endl;
        }

        for(unsigned j=0; j<swarm.size(); j++)
        {
            swarm[j].UpdatePosition0();
            swarm[j].UpdateVelocity0();
            PS_Compute(data, swarm[j]);

            swarm[j].UpdatePosition1();
            swarm[j].UpdateVelocity1();
            PS_Compute(data, swarm[j]);

            swarm[j].UpdatePosition0();
            swarm[j].UpdateVelocity0();
            PS_Compute(data, swarm[j]);

            swarm[j].UpdatePosition2();
            swarm[j].UpdateVelocity2();
            PS_Compute(data, swarm[j]);

            swarm[j].UpdatePosition0();
            swarm[j].UpdateVelocity0();
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
            }
        if(b)
        {
            cout << "-----------------------------------------------------" << endl;

            t2 = clock();
        }
        q++;
    }

    t3 = clock();
    cout << (t2 - t1)/CLOCKS_PER_SEC << " " << (t3 - t1)/CLOCKS_PER_SEC << endl;

    return 0;
}
