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

int main()
{
    cout << "Hello world!" << endl;
    return 0;
}
