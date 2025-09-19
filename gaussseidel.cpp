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
void create_vectorf(int k, double f[520])
{
    int n =512;
    double dx = 1/(double)n;
    for(int i = 0; i< 513; i++)
    {
        f[i] = (pow(PI*k,2) + sigma)*sin(PI*k*i*dx);
    }
}
int main()
{
    double u[520], v[520], f[520], res[520];
    int n,k;
    cout << "Enter the number of grid points on the domain( eg: 512,64,32 etc):" << endl;
    cin >> n;
    cout << "Enter the wave number: " << endl;
    cin >> k;
    double dx = 1/(double)n;
    double norm;
    int iter = 0;
    double a = -1/pow(dx,2);
    double b = (2/pow(dx,2)) + sigma;
    //Initializing velocities to zero
    for(int i = 0; i<n+1; i++)
    {
        v[i] = 0;
    }
    create_vectorf(k,f);
    
    ofstream outfile("Gauss_Seidel_Results.txt");
    
    do
    {
        for(int i = 0; i < n+1; i++)
        {
            u[i] = v[i];
        }
        for(int i = 1; i < n; i++)
        {
            v[i] = (f[i] - a*v[i+1] - a*v[i-1])/b;
        }
        compute_residual(n, res, f, v);
        norm = twonorm(n,res);
        iter++;
    } while (norm > epsilon);
       
    cout << "Results printrd in Gauss_Seidel_Results.txt." << endl;
    outfile <<"Iterations:\t" << iter << endl;
    outfile << "Residual norm:\t" << norm << endl;
    outfile << "Your code is running Fine. Congrats!" << endl;
    outfile.close();
    return 0;
    
}
