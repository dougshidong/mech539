MODULE SOLVER
    INTEGER :: NI, NJ
    REAL    :: W, TOL
    CONTAINS
    SUBROUTINE SOLVELAPLACE(U, B, RESI, TIMES, ITER, ISOLV, W)
        IMPLICIT REAL (A-H,O-Z)
    
        REAL, DIMENSION(:) :: U, B, RESI, TIMES
        REAL, ALLOCATABLE  :: UOLD(:)
        IF(ISOLV == 1) ALLOCATE(UOLD(SIZE(U)))
        ERROR = 9999
        RES = 0
        ITER = 0
        CALL CPU_TIME(TSTART)
        DO WHILE(ERROR > TOL)
            ITER = ITER + 1
            IF(ISOLV == 1) UOLD = U
            RES = 0
            DO I = 2, NI - 1
                DO J = 2, NJ - 1
                    IJ = IGETIJ(I, J, NI)
                    IM = IGETIJ(I - 1, J, NI)
                    IP = IGETIJ(I + 1, J, NI)
                    JM = IGETIJ(I, J - 1, NI)
                    JP = IGETIJ(I, J + 1, NI)
                    UTEMP = U(IJ)
                    IF(ISOLV == 1) U(IJ) = ( B(IJ) - (UOLD(IM) + UOLD(IP) &
                                                    + UOLD(JM) + UOLD(JP)) ) / (-4.0)
                    IF(ISOLV == 2) U(IJ) = ( B(IJ) - (U(IM) + U(IP) &
                                                    + U(JM) + U(JP)) ) / (-4.0)
                    IF(ISOLV == 3) THEN
                       UTEMP = U(IJ)
                       UGS = ( B(IJ) - (U(IM) + U(IP) + U(JM) + U(JP)) ) / (-4.0)
                       U(IJ) = (1 - W) * U(IJ) + W * UGS
                    ENDIF
                    DU = ABS(UTEMP - U(IJ))
                    RES = MAX(RES, DU)
                ENDDO
            ENDDO
            ERROR = RES
            RESI(ITER) = RES
            CALL CPU_TIME(TIMES(ITER))
            IF(MOD(ITER, 10000) == 0) WRITE(*,*) ITER, RES
        ENDDO
        IF(ISOLV == 1) DEALLOCATE(UOLD)
        TIMES(:) = TIMES(:) - TSTART
        WRITE(*,*) ITER, RES, TIMES(ITER)
    END SUBROUTINE
    !***********************************************************************
    FUNCTION IGETIJ(I, J, NI)
        INTEGER :: IGETIJ, I, J, NI
        IGETIJ = (I - 1) * NI + J
    END FUNCTION
END MODULE
