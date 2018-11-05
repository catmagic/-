#include <iostream>
#include<math.h>
using namespace std;

double convert_F_to_K(double T_f)
{
    return (T_f-32.0)*5.0/9.0+273.15;
}

double convert_psi_to_bar(double P_psi)
{
    return P_psi/14.504;
}

void   normalization_P(double* P,double * P_crit,double *P_r ,int n)
{
    for(int i=0; i<n ;i++)
    {
        P_r[i]=P[i]/P_crit[i];
    }
}

void   normalization_T(double* T,double * T_crit,double *T_r ,int n)
{
    for(int i=0; i<n ;i++)
    {
        T_r[i]=T[i]/T_crit[i];
    }
}
void Omega_a_calc(double Omega_a0,double* omega,double* Omega_a,double* T_r, int n)
{
    double temp_calc;
    for(int i=0; i<n;i++)
    {
        temp_calc=1+(0.37464+1.54226*omega[i]-0.26992*omega[i]*omega[i])*(1-sqrt(T_r[i]));
        Omega_a[i]=Omega_a0*temp_calc*temp_calc;
    }
}
void Omega_b_calc(double Omega_b0,double* Omega_a, int n)
{
    for(int i=0; i<n;i++)
    {
        Omega_a[i]=Omega_b0;
    }
}

int main()
{
    cout << "Hello world!" << endl;
    double p[1],p_crit[1],p_r[1];
    p[0]=1.0;
    p_crit[0]=4.0;
    normalization_P(p,p_crit,p_r,1);
    cout<<p_r[0];
    return 0;
}
