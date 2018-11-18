#include <iostream>
#include<math.h>
#include <complex>
using namespace std;
#define  EPSILON 1e-10
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
void calc_A(double** A_i_j ,double *c,double &A,int n)
{
    A=0;
    for(int i=0; i<n;i++)
    {
        for(int j=0; j<n;j++)
        {
            A+=A_i_j[i][j]*c[i]*c[j];
        }
    }
}
void calc_B(double* B_i ,double *c,double &B,int n)
{
    B=0;
    for(int i=0; i<n;i++)
    {
        B+=B_i[i]*c[i];
    }
}
void calc_A_i(double* Omega_A,double* P_r,double* T_r,double *A_i ,int n)
{
    for(int i=0; i<n;i++)
    {
        A_i[i]=Omega_A[i]*P_r[i]/(T_r[i]*T_r[i]);
    }
}
void calc_B_i(double* Omega_B,double* P_r,double* T_r,double *B_i ,int n)
{
    for(int i=0; i<n;i++)
    {
        B_i[i]=Omega_B[i]*P_r[i]/T_r[i];
    }
}
void calc_A_i_j(double** Betta,double* A_i,double** A_i_j ,int n)
{
    for(int i=0; i<n;i++)
    {
        for(int j=0; j<n;j++)
        {
            A_i_j[i][j]=(1-Betta[i][j])*sqrt((A_i[i]*A_i[j]));
        }
    }
}
void cube_root(complex<double> a,complex<double> root[3])
{
 double angle=arg(a),radius=abs(a),rotation=1;
 if(angle<0.0){rotation=-1;}
 root[0]=polar(cbrt(radius),angle/3.0);
 root[1]=polar(cbrt(radius),(angle+rotation*2*M_PI)/3.0);
 root[2]=polar(cbrt(radius),(angle-rotation*2*M_PI)/3.0);
}
void cube_solver(double E2,double E1,double E0,double root[3],int count_double_root)
{
    double p,q;
    p=(3.0*E1-E2*E2)/3.0;                   //canonical form
    q=(2.0*E2*E2*E2-9.0*E2*E1+27.0*E0)/27.0;//cube equation
    double Q;
    Q=(p/3.0)*(p/3.0)*(p/3.0)+(q/2.0)*(q/2.0);
    complex<double> temp_root1[3],temp_root2[3],a;
    if(abs(Q)<EPSILON)//Q=0
    {
        //matching roots
        if((abs(p)<EPSILON)&&(abs(q)<EPSILON))
        {
            //3 matching roots
            count_double_root=3;
            root[0]=-E2/3.0;
            root[1]=-E2/3.0;
            root[2]=-E2/3.0;
        }
        else
        {
            //2 matching roots
            count_double_root=3;
            a.real()=q/2.0;
            a.imag()=0.0;
            cube_root(a,temp_root1);
            cube_root(a,temp_root2);
            root[0]=-2.0*cbrt((q/2.0))-E2/3.0;
            root[1]=-temp_root1[1].real()-temp_root1[2].real()-E2/3.0;
            root[2]=-temp_root1[1].real()-temp_root1[2].real()-E2/3.0;
        }
    }
    else
    {
        if(Q>0.0)
        {
            //1 real ,2 complex root
            count_double_root=1;
            root[0]=cbrt((-q/2.0)+sqrt(Q))+cbrt((-q/2.0)-sqrt(Q))-E2/3.0;
            root[1]=cbrt((-q/2.0)+sqrt(Q))+cbrt((-q/2.0)-sqrt(Q))-E2/3.0;
            root[2]=cbrt((-q/2.0)+sqrt(Q))+cbrt((-q/2.0)-sqrt(Q))-E2/3.0;
        }
        else
        {
            //3 real different root
            count_double_root=3;
            a.real()=-q/2.0;
            a.imag()=sqrt((-Q));
            cube_root(a,temp_root1);
            a.real()=-q/2.0;
            a.imag()=-sqrt((-Q));
            cube_root(a,temp_root2);
            root[0]=temp_root1[0].real()+temp_root2[0].real()-E2/3.0;
            root[1]=temp_root1[1].real()+temp_root2[1].real()-E2/3.0;
            root[2]=temp_root1[2].real()+temp_root2[2].real()-E2/3.0;
        }
    }
}
//void phi(double* c,int n)
int main()
{
    cout << "Hello world!" << endl;
    double p[1],p_crit[1],p_r[1];
    int i=0;
    p[0]=1.0;
    p_crit[0]=4.0;
    normalization_P(p,p_crit,p_r,1);
    complex<double> a,root[3];
    double root_d[3];
    a.real()=0.0;
    a.imag()=-1.0;
    cube_root(a,root);
     cube_solver(-8.0,21.0,-18.0,root_d,i);
      cout<<root_d[0]<<"   "<<'\n';
    cout<<root_d[1]<<"   "<<'\n';
    cout<<root_d[2]<<"   "<<'\n';
    cout<<root[0].real()<<"   "<<root[0].imag()<<'\n';
    cout<<root[1].real()<<"   "<<root[1].imag()<<'\n';
    cout<<root[2].real()<<"   "<<root[2].imag()<<'\n';
    cout<<p_r[0];
    return 0;
}
