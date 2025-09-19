#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>  

#define PI 3.141592653589793
#define sigma 1
#define iterations 2
//#define k 1
#define epsilon 0.000001

using namespace std;

class Vcycle
{
    public:
    double n; // Number of grid points
    double f[520];
    double F[520];
    double u[520];
    double v[520];
    double u_exact[520];
    double v_corr[520];
    double v1[520];
    double res[520];
    double err[520];
    double err_inter[520];
    double A[200][200];
    double dx;

    void exactsolution(int k, int n, double u_exact[520])
    {
        double dx = 1/((double)n-1);
        for(int i = 0; i< n; i++)
        {
            u_exact[i] = sin(k*PI*i*dx);
        }
    }
    void createvectorf(int k, int n, double f[520])
    {
        double dx = 1/((double)n-1);
        for(int i = 1; i < n; i++)
        {
            f[i] = (pow(PI*k,2)+sigma)*sin(k*PI*i*dx);
        }
    }
    void gaussseidel(double u[520], double v[520], double f[520])
    {
        double dx = 1/((double)n-1);
        double a = -pow(dx,-2);
        double b = 2*pow(dx,-2) + sigma;
        for(int i = 0; i < n; i++)
        {
            v[i] = u[i];
        }
        for(int  iter = 0; iter < iterations; iter++)
        {
            for(int i = 1; i < n-1; i++)
            {
                v[i] = (f[i] - a*v[i+1] - a*v[i-1])/b;
            }
        }
    }
    void compute_residual(int n, double res[520], double f[520], double v[520])
    {
        double dx = 1/((double)n-1);
        double a = -pow(dx,-2);
        double b = 2*pow(dx,-2) + sigma;
        for(int i = 1; i< n-1; i++)
        {
            res[i] = f[i] - a*v[i-1] - b*v[i] - a*v[i+1];
        }
    }
    void prolongation(int n,double res_f[520], double res_c[520])
    {
        for(int i = 1; i< n-1; i++)
        {
            res_c[i] = 0.25*(res_f[2*i+1] + 2*res_f[2*i] + res_f[2*i-1]);
        }
    }
    void solve_error_eqn_coarsest_grid(double err[520], double f[520])
    {
        double dx = 0.5;
        double b = 2*pow(dx,-2) + sigma;
        err[1] = f[1]/b;
        err[0] = 0.0;
        err[2] = 0.0;
    }
    void interpolate_error(int n1, int n2, double err_c[520], double err_f[520])
    {
        for(int i = 0; i< n1; i++)
        {
            err_f[2*i] = err_c[i];
        }
        for(int i =1; i<n2-1; i=i+2)
        {
            err_f[i] = 0.5*(err_f[i+1] + err_f[i-1]);
        }
    }
    void correction(int n, double v_corr[520], double v[520], double err_inter[520])
    {
        for(int i = 0; i < n; i++)
        {
            v_corr[i] = v[i] + err_inter[i];
        }
    }
    void infinitynorm(int n, double m[520], double p[520], double& value)
    {
        value = 0.0;
        for (int i = 0; i < n; i++)
        {
            if(fabs(m[i]-p[i])>value)
                value = fabs(m[i]-p[i]);
        }
    }
    double twonorm(int n, double res[520])
    {
    	double value = 0.0;
    	for(int i =0; i<n; i++)
    	{
    		value = value + pow(res[i],2);
    	}
    	return sqrt(value);
    }

};

// V-cycle function
void V_cycle(Vcycle* level, int a)
{
    for(int i = a; i > 1; i--)
    {
        level[i].gaussseidel(level[i].u, level[i].v, level[i].f);
        level[i].compute_residual(level[i].n, level[i].res, level[i].f, level[i].v);
        level[i-1].prolongation(level[i-1].n, level[i].res, level[i-1].f);
    }
    // Solving error equation at the coarsest
    level[1].solve_error_eqn_coarsest_grid(level[1].v1, level[1].f);
    // Going upwards from coarsest grid to the finest grid
    for(int i = 1; i < a; i++)
    {
        level[i].interpolate_error(level[i].n, level[i+1].n, level[i].v1, level[i+1].err_inter);
        level[i+1].correction(level[i+1].n, level[i+1].v_corr, level[i+1].v, level[i+1].err_inter);
        level[i+1].gaussseidel(level[i+1].v_corr, level[i+1].v1, level[i+1].f);
    }
}
void FMG(int k, Vcycle* level, int a)
{
    level[1].createvectorf(k,level[1].n, level[1].f);
    level[1].solve_error_eqn_coarsest_grid(level[1].u, level[1].F);
    for(int i = 1; i < a; i++)
    {
        level[i].interpolate_error(level[i].n, level[i+1].n, level[i].u, level[i+1].u);
        level[i+1].createvectorf(k,level[i+1].n, level[i+1].f);
        V_cycle(level, i+1);
    }
    
}

int main()
{
    Vcycle level[20];
    cout << "Enter the no. of grid points in Finest grid(eg: 512,256,64 etc)" << endl;
    int g,k;
    cin >> g;
    cout << "Enter the wave number:" << endl;
    cin >> k;
    int N = log(g)/log(2);
    //int N = 9;
    double norm = 0.0;
    int Vcycle_iter = 0;
    // Total grid points in a level
    for(int i = 1; i<10; i++)
    {
        level[i].n = pow(2,i)+1;
    }
    // Initial conditions in each level
    for(int i = 2; i < 10; i++)
    {
        for(int j = 0; j<level[i].n; j++)
        {
            level[i].u[j] = 0.0;
        }
    }
    level[N].exactsolution(k,level[N].n, level[N].u_exact);
    level[N].createvectorf(k,level[N].n,level[N].f);
    
    
    ofstream outfile("FMG_Results.txt");
    outfile << "Residual Norm:\t\tIterations:" << endl;
    FMG(k,level, N);
        
    do
    {
        for(int i = 0; i < level[N].n; i++)
        {
            level[N].u[i] = level[N].v1[i];
        }

        V_cycle(level, N);
        norm = level[N].twonorm(level[N].n,level[N].res);
        for(int i = 2; i < N-1; i++)
        {
            for(int j = 0; j < level[i].n; j++)
            {
                level[i].u[j] = 0.0;
            }
        }
        
        Vcycle_iter++;
        outfile << norm << "\t\t\t" <<  Vcycle_iter  << endl;
    } while (norm > epsilon);
    
    cout << "Results printed in FMG_Results.txt file." << endl;
    outfile << "\n\nExact solution:\t\t" << "Numerical solution:" << endl;
    for(int i = 0; i< level[N].n; i++)
    {
        outfile << level[N].u_exact[i] << "\t\t" << level[N].v1[i] << endl;
    }

    outfile << "Residual two Norm of the velocities is: " << endl;
    outfile << norm << endl;

    outfile << "NO.of V cycle iterations in FMG for convergence is: " << endl;
    outfile << Vcycle_iter << endl;

    outfile << "Your code is running Fine. Congrats!" << endl;

 
    outfile.close();

    return 0;
}
