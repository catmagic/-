! A fortran95 program for G95
! By WQY
program main

    real T_f,T_k
    t_f=500.
    call covert_F_to_K(T_f,T_k)
    write(*, *)T_f
    write(*, *)T_k
end

subroutine covert_F_to_K(T_f,T_k)
    real T_f,T_k
    intent(in)  T_f
    intent(out) T_k
    T_k = (T_f - 32) * 5. / 9. + 273.15
end subroutine
