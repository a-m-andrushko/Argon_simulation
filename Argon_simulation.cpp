#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>

using namespace std;

int n, m, e, S_o, S_d, S_out, S_xyz;
double R, f, L, a, T_0, tau;
int N;
double k = 8.31 * 0.001; // BOLTZMANN CONSTANT

struct atom
{
    double r[3];
    double E[3];
    double p[3];
    double F[3];
};

double calc_VF(atom* atoms, double* V_P, double* V_S, double* P)
{
    for(int i = 0; i < N; i++)
    {
        atoms[i].F[0] = 0;
        atoms[i].F[1] = 0;
        atoms[i].F[2] = 0;
    }
    *V_S = 0;
    *V_P = 0;
    *P = 0;
    double distance = 0;
    double r_ij = 0;

    for(int i = 0; i < N; i++)
    {
        distance = sqrt(atoms[i].r[0] * atoms[i].r[0] + atoms[i].r[1] * atoms[i].r[1] + atoms[i].r[2] * atoms[i].r[2]);

        if(distance >= L)
        {
            *V_S += 1 / 2. * f * (distance - L) * (distance - L);
            atoms[i].F[0] += f * (L - distance) * atoms[i].r[0] / distance;
            atoms[i].F[1] += f * (L - distance) * atoms[i].r[1] / distance;
            atoms[i].F[2] += f * (L - distance) * atoms[i].r[1] / distance;
        }

        *P += 1 / (4 *M_PI * L * L) * sqrt(atoms[i].F[0] * atoms[i].F[0] + atoms[i].F[1] * atoms[i].F[1] + atoms[i].F[2] * atoms[i].F[2]);

        for(int j = 0; j < i; j++)
        {
            r_ij = 0;
            r_ij = sqrt((atoms[i].r[0] - atoms[j].r[0]) * (atoms[i].r[0] - atoms[j].r[0]) + (atoms[i].r[1] - atoms[j].r[1]) * (atoms[i].r[1] - atoms[j].r[1]) + (atoms[i].r[2] - atoms[j].r[2]) * (atoms[i].r[2] - atoms[j].r[2]));
            *V_P += pow(R / r_ij, 12) - 2 * pow(R / r_ij, 6);
            atoms[i].F[0] += 12 * (pow(R / r_ij, 12) - pow(R / r_ij, 6)) / (r_ij * r_ij) * (atoms[i].r[0] - atoms[j].r[0]);
            atoms[i].F[1] += 12 * (pow(R / r_ij, 12) - pow(R / r_ij, 6)) / (r_ij * r_ij) * (atoms[i].r[1] - atoms[j].r[1]);
            atoms[i].F[2] += 12 * (pow(R / r_ij, 12) - pow(R / r_ij, 6)) / (r_ij * r_ij) * (atoms[i].r[2] - atoms[j].r[2]);
            atoms[j].F[0] += -12 * (pow(R / r_ij, 12) - pow(R / r_ij, 6)) / (r_ij * r_ij) * (atoms[i].r[0] - atoms[j].r[0]);
            atoms[j].F[1] += -12 * (pow(R / r_ij, 12) - pow(R / r_ij, 6)) / (r_ij * r_ij) * (atoms[i].r[1] - atoms[j].r[1]);
            atoms[j].F[2] += -12 * (pow(R / r_ij, 12) - pow(R / r_ij, 6)) / (r_ij * r_ij) * (atoms[i].r[2] - atoms[j].r[2]);
        }
    }

    return 0;
}

double calc_T(atom* atoms, double* H, double* V_P, double* V_S, double* T, double* E_k)
{
    *E_k = 0;
    
    for(int i = 0; i < N; i++)
        *E_k += (atoms[i].p[0] * atoms[i].p[0] + atoms[i].p[1] * atoms[i].p[1] + atoms[i].p[2] * atoms[i].p[2]) / (2 * m);

    *H = *E_k + *V_S + *V_P;

    *T = 2 / (3 * N * k) * (*E_k);

    return 0;
}

double time_step(atom* atoms, double* V_P, double* V_S, double* P)
{
    for(int i = 0; i < N; i++)
    {
        atoms[i].p[0] += 0.5 * atoms[i].F[0] * tau;
        atoms[i].p[1] += 0.5 * atoms[i].F[1] * tau;
        atoms[i].p[2] += 0.5 * atoms[i].F[2] * tau;

        atoms[i].r[0] += 1. / m * atoms[i].p[0] * tau;
        atoms[i].r[1] += 1. / m * atoms[i].p[1] * tau;
        atoms[i].r[2] += 1. / m * atoms[i].p[2] * tau; 
    }

    calc_VF(atoms, V_P, V_S, P);

    for(int i = 0; i < N; i++)
    {
        atoms[i].p[0] += 0.5 * atoms[i].F[0] * tau;
        atoms[i].p[1] += 0.5 * atoms[i].F[1] * tau;
        atoms[i].p[2] += 0.5 * atoms[i].F[2] * tau;
    }

    return 0;
}

int main(int argc, char* argv[])
{
    // OPEN INPUT FILE WITH PARAMETERS
    ifstream infile("parameters.txt");
    if(!infile)
    {
        cerr << "Unable to open file!";
        return 1;
    }

    // READ AND ASSIGN VALUES
    infile >> n >> m >> e >> R >> f >> L >> a >> T_0 >> tau >> S_o >> S_d >> S_out >> S_xyz;

    // CLOSE INPUT FILE
    infile.close();
    

    N = n * n * n;

    atom atoms[N];


    // BASE VECTORS
    double b_0[3] = {a, 0, 0};
    double b_1[3] = {a / 2., a * sqrt(3.) / 2., 0};
    double b_2[3] = {a / 2., a * sqrt(3.) / 6., a * sqrt(2. / 3.)};    

    // CALCULATE POSITION (r) FOR EACH ATOM
    for(int i_0 = 0; i_0 < n; i_0++)
    {
        for(int i_1 = 0; i_1 < n; i_1++)
        {
            for(int i_2 = 0; i_2 < n; i_2++)
            {
                atoms[i_0 + i_1 * n + i_2 * n * n].r[0] = (i_0 - ((n - 1) / 2.)) * b_0[0] + (i_1 - ((n - 1) / 2.)) * b_1[0] + (i_2 - ((n - 1) / 2.)) * b_2[0];

                atoms[i_0 + i_1 * n + i_2 * n * n].r[1] = (i_0 - ((n - 1) / 2.)) * b_0[1] + (i_1 - ((n - 1) / 2.)) * b_1[1] + (i_2 - ((n - 1) / 2.)) * b_2[1];

                atoms[i_0 + i_1 * n + i_2 * n * n].r[2] = (i_0 - ((n - 1) / 2.)) * b_0[2] + (i_1 - ((n - 1) / 2.)) * b_1[2] + (i_2 - ((n - 1) / 2.)) * b_2[2];
            }
        }
    }


    // GENERATE KINETIC ENERGY
    srand48(time(0));
    double sum_E = 0;

    random_device rd;  // SEED
    mt19937 gen(rd()); // MERSENNE TWISTER ENGINE
    uniform_real_distribution<> dis(0.0, 1.0);

    double random_number_x = 0, random_number_y = 0, random_number_z = 0;

    for(int i = 0; i < N; i++)
    {
        random_number_x = dis(gen);
        random_number_y = dis(gen);
        random_number_z = dis(gen);

        atoms[i].E[0] = -1 * 0.5 * k * T_0 * log(random_number_x);
        atoms[i].E[1] = -1 * 0.5 * k * T_0 * log(random_number_y);
        atoms[i].E[2] = -1 * 0.5 * k * T_0 * log(random_number_z);
    }


    // CALCULATE MOMENTUM (p) FOR EACH ATOM
    uniform_int_distribution<> distr(1, 100);

    int random_number_px = 0, random_number_py = 0, random_number_pz = 0;

    for(int i = 0; i < N; i++)
    {
        random_number_px = distr(gen);
        random_number_py = distr(gen);
        random_number_pz = distr(gen);

        if(random_number_px % 2 == 0)    
            atoms[i].p[0] = sqrt(2 * m * atoms[i].E[0]);
        else
            atoms[i].p[0] = -sqrt(2 * m * atoms[i].E[0]);

        if(random_number_py % 2 == 0)    
            atoms[i].p[1] = sqrt(2 * m * atoms[i].E[1]);
        else
            atoms[i].p[1] = -sqrt(2 * m * atoms[i].E[1]);

        if(random_number_pz % 2 == 0)    
            atoms[i].p[2] = sqrt(2 * m * atoms[i].E[2]);
        else
            atoms[i].p[2] = -sqrt(2 * m * atoms[i].E[2]);       
    }


    // CALCULATE INITIAL FORCES AND POTENTIALS
    double V_P = 0;
    double V_S = 0;
    double P = 0;    
    
    calc_VF(atoms,&V_P, &V_S, &P);


    // OPEN OUTPUT FILES FOR WRITING
    ofstream outfile1("state.txt");
    if (!outfile1) {
        cerr << "Unable to open state.txt!";
        return 1;
    }

    ofstream outfile2("xyz.txt");
    if (!outfile2) {
        cerr << "Unable to open xyz.txt!";
        return 1;
    }

    // SET PRECISION
    outfile1 << fixed << setprecision(5);

    outfile2 << fixed << setprecision(5);

    
    double T = 0;
    P = 0;
    double H = 0;
    double E_k = 0;

    double T_sum = 0, P_sum = 0;
    
    // SIMULATION LOOP
    for(int t = 0; t < S_o + S_d; t++)
    {
        time_step(atoms, &V_P, &V_S, &P);

        calc_T(atoms, &H, &V_P, &V_S, &T, &E_k);
        
        if(t >= S_o)
        {
            T_sum += T;
            P_sum += P;

            if(t % S_out == 0)
                outfile1 << setw(12) << t << setw(12) << H << setw(12) << V_P + V_S << setw(12) << T << setw(12) << P << "\n";
        }

        if(t % S_xyz == 0)
        {
            for(int i = 0; i < N; i++)
                outfile2 << setw(12) << atoms[i].r[0] << setw(12) << atoms[i].r[1] << setw(12) << atoms[i].r[2] << "\n";

            outfile2 << "\n" << "\n";
        }
    }

    cout << "T_avg: " << T_sum / S_d << endl;
    cout << "P_avg: " << P_sum / S_d << endl;

    // CLOSE OUTPUT FILES
    outfile1.close();
    outfile2.close();

    return 0;
}