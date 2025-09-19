#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>  

#define PI 3.141592653589793
#define sigma 1
//#define k 1
#define epsilon 0.000001
using namespace std;

double infinity_norm(double u[520], double v[520])
{
    double sum = 0;
    for(int i = 0; i<513; i++)
    {
        if(fabs(u[i] - v[i]) > sum)
        sum = fabs(u[i] - v[i]);
    }
    return sum;
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

void create_vectorf(int k, double f[520])
{
    int n =512;
    double dx = 1/(double)n;
    for(int i = 0; i< 513; i++)
    {
        f[i] = (pow(PI*k,2) + sigma)*sin(PI*k*i*dx);
    }
}
void compute_residual(int n, double res[520], double f[520], double v[520])
    {
        double dx = 1/((double)n);
        double a = -pow(dx,-2);
        double b = 2*pow(dx,-2) + sigma;
        for(int i = 1; i< n-1; i++)
        {
            res[i] = f[i] - a*v[i-1] - b*v[i] - a*v[i+1];
        }
    }

int main()
{
    double u_old[520], u_new[520], v[520], f[520], res[520];
    int n;
    cout << "Enter the number of grid points on the domain" << endl;
    cin >> n;
    double dx = 1/(double)n;
    double norm;
    int iter = 0;
    double a = -1/pow(dx,2);
    double b = (2/pow(dx,2)) + sigma;
    double w;
    int k;
    cout << "\nEnter the wave number(k):" << endl;
    cin >> k;
    cout << "\nEnter the weighting factor(w):" << endl;
    cin >> w;
    create_vectorf(k,f);
    ofstream outfile("Weighted_Jacobi_Results.txt");
    //Initial guess values
    for(int i = 0; i < n+1; i++)
    {
        u_new[i] = 0;
    }
    do
    {
       for(int i = 0; i< n+1; i++)
       {
        u_old[i] = u_new[i];
       }
       for(int i = 1; i<n; i++)
       {
        u_new[i] = (1-w)*u_old[i] + w*(f[i] - a*u_old[i+1] - a*u_old[i-1])/b;
       }
       iter++;
       compute_residual(n, res,f, u_new);
       norm  = twonorm(n, res);

    } while (norm > epsilon);
    cout << "Results are printed in WEighted_Jacobi_Results.text" << endl;
    outfile << "\nResidual Norm:\t" << norm << endl;
    outfile << "Iterations:\t" << iter << endl;
    outfile << "Your code is running Fine. Congrats!" << endl;
    return 0;
}
