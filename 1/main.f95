! A fortran95 program for G95
! By WQY
program main

    real T_f,T_k,P_psi,P_bar
    T_f     =   500.
    P_psi   =   1400.
    call covert_F_to_K(T_f,T_k)
    write(*, *)T_f
    write(*, *)T_k
    call covert_psi_to_bar(P_psi,P_bar)
    write(*, *)P_psi
    write(*, *)P_bar
end

subroutine covert_F_to_K(T_f,T_k)
    real T_f,T_k
    intent(in)  T_f
    intent(out) T_k
    T_k = (T_f - 32) * 5. / 9. + 273.15
end subroutine

subroutine covert_psi_to_bar(P_psi,P_bar)
    real P_psi,P_bar
    intent(in)  P_psi
    intent(out) P_bar
    P_bar = P_psi/14.504
end subroutine
subroutine omega_a_calc(T_ri,omega,n,Omega_A)
    integer, parameter :: n
    real ::T_ri(n),omega(n),Omega_A(n)
    intent(in)  n,T_ri,omega
    intent(out) P_bar
    P_bar = P_psi/14.504
end subroutine
