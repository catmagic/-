#include <iostream>

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
