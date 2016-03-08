module globals
    double precision    :: del_min = 0.025D0
    double precision    :: conv_tol = 2e-10
    double precision    :: M_in
    double precision, parameter :: gam = 1.40d0
    double precision, parameter :: R = 1.0d0
    double precision, parameter :: Tt_in = 1.0d0
    double precision, parameter :: p_in = 1.0d0
    double precision :: T_in, a_in, U_in, a02

    contains

    subroutine eval_constants
        T_in = Tt_in / (1.0d0 + (gam - 1.0d0) / 2.0d0 * M_in * M_in)
        a_in = sqrt(gam * R * T_in)
        U_in = M_in * a_in
        a02  = a_in**2 + 0.50d0 * (gam - 1.0d0) * U_in * U_in
    end subroutine

end module globals
