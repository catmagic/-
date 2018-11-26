#include <iostream>
#include<stdlib.h>
#include <fstream>
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
void Omega_a_calc(double Omega_a0,double* omega,double * Omega_a,double* T_r, int n)
{
    double temp_calc;
    for(int i=0; i<n;i++)
    {
        temp_calc=1+(0.37464+1.54226*omega[i]-0.26992*omega[i]*omega[i])*(1-sqrt(T_r[i]));
        Omega_a[i]=Omega_a0*temp_calc*temp_calc;
    }
}
void Omega_b_calc(double Omega_b0,double* Omega_b, int n)
{
    for(int i=0; i<n;i++)
    {
        Omega_b[i]=Omega_b0;
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
        cout<<A_i[i]<<'\n';
    }
    cout<<'\n';
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
            //cout<<A_i_j[i][j]<<"    "<<A_i[i]<<"    "<<A_i[j]<<"     "<<Betta[i][j]<<'\n';
        }
    }
}
void calc_V(double *k_i,double *z,double &v,int n)
{
    double k_max=k_i[0],k_min=k_i[0],I_max,I_min,l,r;
    for(int i=0;i<n;i++)
    {
        if(k_max<k_i[i]){k_max=k_i[i];}
        if(k_min>k_i[i]){k_min=k_i[i];}
    }
    I_min=1/(1-k_max);//+infinity
    I_max=1/(1-k_min);//-infinity
    l=I_min;
    r=I_max;
     double sum;
    v=(l+r)/2.0;
    do{

       sum=0.0;
        for(int i=0;i<n;i++)
        {
            sum+=((k_i[i]-1)*z[i])/(1+(k_i[i]-1)*v);
        }
        if(abs(sum)<EPSILON){return;}
    if(sum>0.0){l=v;v=(l+r)/2.0;}
    else {r=v;v=(l+r)/2.0;}
    }while(abs(sum)>EPSILON);

}
void calc_x(double *k_i,double *z,double*x, double v,int n)
{
    for(int i=0;i<n;i++)
    {
        x[i]=z[i]/(1+(k_i[i]-1)*v);
    }
}
void calc_y(double *k_i,double *z,double*y, double v,int n)
{
    for(int i=0;i<n;i++)
    {
        y[i]=z[i]*k_i[i]/(1+(k_i[i]-1)*v);
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
void normalization_c(double *c,int n)
{
    double sum=0.0;
    for(int i=0;i<n;i++)
    {
        sum+=c[i];
    }
    for(int i=0;i<n;i++)
    {
        c[i]=c[i]/sum;
    }
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
void calc_phi(double C,double * S_i,double *B_i,double* Phi_i,double m1,double m2,double Z,double B,double A,int n)
{
    double M=log((Z+m2*B)/(Z+m1*B));
    double *Ln_phi_i;
    Ln_phi_i=(double *)malloc(n*sizeof(double));

    for(int i=0;i<n;i++)
    {
        Ln_phi_i[i]=log(C/(Z-B))+M*A/((m1-m2)*B)*(2.0*S_i[i]/A-B_i[i]/B)+B_i[i]*(Z-C)/B;
        Phi_i[i]=exp(Ln_phi_i[i]);
    }
}
void coeff(double E2,double E1,double E0,double m1,double m2,double A,double B,double C)
{
    E2=(m1+m2)*B-(B+C);
    E1=A+m1*m2*B*B-(m1+m2)*B*(B+C);
    E0=-A*B-m1*m2*B*B*(B+C);
}
void calc_c(double *c_i,double c,int n)
{
    c=0.0;
    for(int i=0;i<n;i++)
    {
        c+=c_i[i];
    }
}
void calc_si(double* si,double** A_ij,double *c,int n)
{
    for(int i=0;i<n; i++)
    {
        si[i]=0.0;
    }

    for(int i=0;i<n ;i++)
    {

        for(int j=0;j<n;j++)
        {
            si[i]+=A_ij[i][j]*c[j];
        }
    }
}
double delta(int i,int j)
{
    return double(i==j);
}
void practical_liq(double* z,double *K,double V,double** practical_x_practical_K,int n)
{
    double sum=0.0;
    double temp;
    for(int i=0;i<n;i++)
    {
        sum+=((K[i]-1)*(K[i]-1)*z[i])/((1+(K[i]-1)*V)*(1+(K[i]-1)*V));
    }
    for(int i=0;i<n;i++)
    {
        temp=z[i]/((1+(K[i]-1)*V)*(1+(K[i]-1)*V));
        for(int j=0;j<n;j++)
        {
            practical_x_practical_K[i][j]=temp*(-V*delta(i,j)-(K[i]-1)*temp/sum);
        }
    }
}
void practical_gas(double* z,double *K,double V,double** practical_y_practical_K,int n)
{
    double sum=0.0;
    double temp;
    for(int i=0;i<n;i++)
    {
        sum+=((K[i]-1)*(K[i]-1)*z[i])/((1+(K[i]-1)*V)*(1+(K[i]-1)*V));
    }
    for(int i=0;i<n;i++)
    {
        temp=z[i]/((1+(K[i]-1)*V)*(1+(K[i]-1)*V));
        for(int j=0;j<n;j++)
        {
            practical_y_practical_K[i][j]=temp*((1-V)*delta(i,j)-K[i]*(K[i]-1)*temp/sum);
        }
    }
}
void practical_Z(double Z,double M1,double M2,double E2,double E1,double A,double B,double C,double *B_i,double *S_i,double *practical_Z_practical_C,int n)
{
    for(int i=0;i<n ; i++)
    {
        double dE2_dc,dE1_dc,dE0_dc;
        dE2_dc=(M1+M2-1.0)*B_i[i]-1.0;
        dE1_dc=2.0*S_i[i]+2.0*M1*M2*B*B_i[i]-(M1+M2)*(B_i[i]*(B+C)+B*(B_i[i]+1.0));
        dE0_dc=-2.0*S_i[i]*B-A*B_i[i]-M1*M2*(2.0*B*B_i[i]*(B+C)+B*B*(B_i[i]+1.0));
        practical_Z_practical_C[i]=-(Z*Z*dE2_dc+Z*dE1_dc+dE0_dc)/(3.0*Z*Z+2.0*Z*E2+E1);
    }
}
void practical_ln_phi(double M1,double M2,double A,double B,double C,double Z,double *B_i,double* S_i,double *practical_Z_practical_C,double **practical_c_practical_K,double **practical_ln_phi_practical_K,double **A_ij, int n)
{
    double W1,W2,W3,W4,W5,W6,M;
    double **practical_ln_phi_practical_c;
    practical_ln_phi_practical_c        =   (double**)malloc(n*sizeof(double*));
    W1=A/((M1-M2)*B);
    M=log((Z+M2*B)/(Z+M1*B));
    for(int i=0;i<n;i++)
    {
        practical_ln_phi_practical_c[i]         =   (double*)malloc(n*sizeof(double));
        W5=2.0*S_i[i]/A-B_i[i]/B;
        for(int j=0;j<n;j++)
        {
            W2=1.0/C-(practical_Z_practical_C[j]-B_i[j])/(Z-B);
            W3=(practical_Z_practical_C[j]+M2*B_i[j])/(Z+M2*B)-(practical_Z_practical_C[j]+M1*B_i[j])/(Z+M1*B);
            W4=M*(2.0*S_i[j]*B-A*B_i[j])/((M1-M2)*B*B);
            W6=(B_i[i]*B_i[j])/(B*B);
            practical_ln_phi_practical_c[i][j]=W2+W1*W5*W3+W4*W5+M*W1*((2.0*A*A_ij[i][j]-4.0*S_i[i]*S_i[j])/(A*A)+W6)-(Z-C)*W6+(B_i[i]/B)*(practical_Z_practical_C[j]-1);
        }
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            practical_ln_phi_practical_K[i][j]=0.0;

            for(int k=0;k<n;k++)
            {
                practical_ln_phi_practical_K[i][j]+=practical_ln_phi_practical_c[i][k]*practical_c_practical_K[k][j];
            }
        }
    }
}
void A_linar_calc(double **A_linar_ij,double** practical_ln_phi_practical_K_liq,double **practical_ln_phi_practical_K_gas,double* K_i,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int  j=0;j<n;j++)
        {
            A_linar_ij[i][j]=delta(i,j)+K_i[j]*(-practical_ln_phi_practical_K_liq[i][j]+practical_ln_phi_practical_K_gas[i][j]);
        }
    }
}
void B_linar_calc(double*B_linar_i, double* phi_liq,double* phi_gas,double* K_i,int n)
{
    for(int i=0;i<n;i++)
    {
        B_linar_i[i]=log(phi_liq[i]/(K_i[i]*phi_gas[i]));
    }
}
void K_calk(double* K_i,double* tetta_i,int n)
{
    for(int i;i<n;i++)
    {
        K_i[i]*=exp(tetta_i[i]);
    }
}
void Gauss(double** A,double* B,double* x,int n)
{
    for(int i=0;i<n;i++)
    {

            bool flag=false;
            int k=i;
            for(;(k<n)&&!flag;k++)
            {
                if(abs(A[k][i])>EPSILON){flag=true;}
            }
            double swap_temp;
            double a,b;
            for(int l=i;l<n;l++)
            {
                swap_temp=A[k][l];
                A[k][l]=A[i][l];
                A[i][l]=swap_temp;
            }
            swap_temp=B[k];
            B[k]=B[i];
            B[i]=swap_temp;
            a=A[i][i];
            for(k=i+1;k<n;k++)
            {
                b=A[k][i];
                for(int l=i;l<n;l++)
                {
                    A[k][l]=a*A[k][l]-b*A[i][l];
                }
                B[k]=a*B[k]-b*B[i];
            }

    }


}
int main()
{
    ifstream data;
    data.open("data.txt");
    cout<<delta(1,0)<<" "<<delta(1,1);
    double *p,*p_crit,*p_r,*s_liq_i;
    double *t,*t_crit,*t_r,*s_gas_i;
    double *omega,*Omega_A,*Omega_B,*A_i,*B_i;
    double **A_ij,**Betta;
    double A_liq=0.0,B_liq=0.0,A_gas=0.0,B_gas=0.0,C_liq,C_gas,temp,T,P,v=0.0,Omega_A0=0.457235529,Omega_B0=0.077796074;
    double *x,*y,*z,*c,*k_i;
    double m1,m2;
    double root_liq[3],root_gas[3];
    m1=1.0+sqrt(2.0);
    m2=1.0-sqrt(2.0);
    double E2_liq,E1_liq,E0_liq,z_liq;
    double E2_gas,E1_gas,E0_gas,z_gas;
    int n,count_root_gas,count_root_liq;
    double *phi_liq_i,*phi_gas_i;
    data>>n>>T>>P;
    s_liq_i     =   (double*)malloc(n*sizeof(double));
    s_gas_i     =   (double*)malloc(n*sizeof(double));
    phi_liq_i   =   (double*)malloc(n*sizeof(double));
    phi_gas_i   =   (double*)malloc(n*sizeof(double));
    p           =   (double*)malloc(n*sizeof(double));
    p_crit      =   (double*)malloc(n*sizeof(double));
    p_r         =   (double*)malloc(n*sizeof(double));
    t           =   (double*)malloc(n*sizeof(double));
    t_crit      =   (double*)malloc(n*sizeof(double));
    t_r         =   (double*)malloc(n*sizeof(double));
    omega       =   (double*)malloc(n*sizeof(double));
    Omega_A     =   (double*)malloc(n*sizeof(double));
    Omega_B     =   (double*)malloc(n*sizeof(double));
    A_i         =   (double*)malloc(n*sizeof(double));
    B_i         =   (double*)malloc(n*sizeof(double));
    x           =   (double*)malloc(n*sizeof(double));
    y           =   (double*)malloc(n*sizeof(double));
    z           =   (double*)malloc(n*sizeof(double));
    c           =   (double*)malloc(n*sizeof(double));
    k_i         =   (double*)malloc(n*sizeof(double));

    A_ij        =   (double**)malloc(n*sizeof(double*));
    Betta       =   (double**)malloc(n*sizeof(double*));
    int i=0,j=0;
    for(i=0;i<n;i++)
    {
        data>>t_crit[i]>>p_crit[i]>>omega[i]>>temp>>c[i]>>z[i];
        t[i]=T;
        p[i]=P;
        k_i[i]=(i<n/2)?(0.5):(2.0);

        A_ij[i]         =   (double*)malloc(n*sizeof(double));
        Betta[i]        =   (double*)malloc(n*sizeof(double));
    }
    normalization_P(p,p_crit,p_r,n);
    normalization_T(t,t_crit,t_r,n);

    normalization_c(z,n);
    calc_V(k_i,z,v,n);
    calc_x(k_i,z,x,v,n);
    calc_y(k_i,z,y,v,n);

    Omega_a_calc(Omega_A0,omega,Omega_A,t_r,n);
    Omega_b_calc(Omega_B0,Omega_B,n);
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            data>>Betta[i][j];
        }
    }

    calc_B_i(Omega_B,p_r,t_r,B_i,n);
    calc_A_i(Omega_A,p_r,t_r,A_i,n);
    calc_A_i_j(Betta,A_i,A_ij,n);
    calc_A(A_ij,x,A_liq,n);
    calc_A(A_ij,y,A_gas,n);
    calc_B(B_i,x,B_liq,n);
    calc_B(B_i,y,B_gas,n);
    calc_c(x,C_liq,n);
    calc_c(y,C_gas,n);
    coeff(E2_liq,E1_liq,E0_liq,m1,m2,A_liq,B_liq,C_liq);
    coeff(E2_gas,E1_gas,E0_gas,m1,m2,A_gas,B_gas,C_gas);
    cube_solver(E2_liq,E1_liq,E0_liq,root_liq,count_root_liq);
    cube_solver(E2_gas,E1_gas,E0_gas,root_gas,count_root_gas);
    z_liq=max(max(root_liq[0],root_liq[1]),root_liq[2]);
    if((z_liq<root_liq[0])&&(root_liq[0]>B_liq)){z_liq=root_liq[0];}
    if((z_liq<root_liq[1])&&(root_liq[1]>B_liq)){z_liq=root_liq[1];}
    if((z_liq<root_liq[2])&&(root_liq[2]>B_liq)){z_liq=root_liq[2];}
    //z min positive root
    z_gas=max(max(root_gas[0],root_gas[1]),root_gas[2]);
    calc_si(s_liq_i,A_ij,x,n);
    calc_si(s_gas_i,A_ij,y,n);
    calc_phi(C_liq,s_liq_i,B_i,phi_liq_i,m1,m2,z_liq,B_liq,A_liq,n);
    calc_phi(C_gas,s_gas_i,B_i,phi_gas_i,m1,m2,z_gas,B_gas,A_gas,n);
    double kriterij=0.0;
    for(int k=0;k<n;k++)
    {
        k+=abs((x[i]*phi_liq_i[i])/(y[i]*phi_liq_i[i])-1.0);
    }
    if(kriterij>EPSILON)
    {
    }
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
