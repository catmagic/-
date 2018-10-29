! A fortran95 program for G95
! By WQY
program main
    integer,parameter::n=2
    real T_f,T_k,P_psi,P_bar
    real T(2),T_crit(2),T_ri(2)
    T_f     =   500.
    P_psi   =   1400.
    T=(/1.,1./)
    T_crit=(/2.,4./)
    call covert_F_to_K(T_f,T_k)
    write(*, *)T_f
    write(*, *)T_k
    call covert_psi_to_bar(P_psi,P_bar)
    write(*, *)P_psi
    write(*, *)P_bar
    call T_ri_calc(T,T_crit,T_ri,n)
    write(*, *)T
    write(*, *)T_crit
    write(*, *)T_ri
end

subroutine covert_F_to_K(T_f,T_k)
    real T_f,T_k

    T_k = (T_f - 32) * 5. / 9. + 273.15
end subroutine

subroutine covert_psi_to_bar(P_psi,P_bar)
    real P_psi,P_bar

    P_bar = P_psi/14.504
end subroutine

subroutine T_ri_calc(T,T_crit,T_ri,n)
    integer n
    real ::T_ri(n),T(n),T_crit(n)

     T_ri=T/T_crit
end subroutine

subroutine omega_a_calc(T_ri,omega,n,Omega_A)
    integer n
    real ::T_ri(n),omega(n),Omega_A(n)


end subroutine
