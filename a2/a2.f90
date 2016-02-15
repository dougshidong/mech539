PROGRAM A2MAIN
    USE SOLVER
    
    IMPLICIT REAL (A-H,O-Z)
    
    REAL, ALLOCATABLE            :: U(:), B(:), Z(:), R(:)
    REAL, DIMENSION (10000000)   :: RESI
    CHARACTER(LEN = 11)          :: FNAME
    INTEGER                      :: ISIZE(3)

    W = 1.5
    ISIZE(1:3) = (/ 100, 200, 400 /)
    ICOND = 1 ! CHECK CONDITION NUMBER
    DO ISOLV = 3, 3
        DO INI = 1, 3
            NI = ISIZE(INI)
            NJ = NI
            NIJ = NI * NJ
            ALLOCATE(U(NIJ), B(NIJ))

            B = 0.0
            
            ! INITIALIZE
            DX = 1.0 / (NI - 1.0)
            DY = 1.0 / (NJ - 1.0)
 
            DO I = 1, NI
              DO J = 1, NJ
                  IJ = IGETIJ(I, J, NI)
                  U(IJ) = 0
                  IF(I == NI) U(IJ) = 1.0
              ENDDO
            ENDDO
            ! SOLVE
            TOL = 1E-15
            CALL SOLVELAPLACE(U, B, RESI, ITER, ISOLV, W)

            WRITE(FNAME,'("RESULTS", I1, I3 )' )  ISOLV, NI
            OPEN(UNIT = 11, FILE = FNAME, FORM='FORMATTED')  
            WRITE(11,*) NI
            WRITE(11,*) NJ 
            DO IJ = 1, NI * NJ
                WRITE(11,*) U(IJ)
            ENDDO

            WRITE(11,*) ITER
            DO I = 1, ITER
                WRITE(11,*) RESI(I)
            ENDDO
            CLOSE(11)

            IF(ICOND == 1) THEN
                ALLOCATE(R(NIJ), Z(NIJ))
                R = 0.0
                Z = 0.0
                JSOLV = 3
                TOL = 1E-15
                ! INITIALIZE Z AND R
                DO I = 2, NI - 1
                    DO J = 2, NJ - 1
                        IJ = IGETIJ(I, J, NI)
                        IM = IGETIJ(I - 1, J, NI)
                        IP = IGETIJ(I + 1, J, NI)
                        JM = IGETIJ(I, J - 1, NI)
                        JP = IGETIJ(I, J + 1, NI)
                        R(IJ) = B(IJ) - (U(IM) + U(IP) + U(JM) + U(JP) - 4 * U(IJ))
                        Z(IJ) = 1.0
                    ENDDO
                ENDDO
                CALL SOLVELAPLACE(Z, R, RESI, ITER, JSOLV, W)
                COND = NORM2(Z) / NORM2(U) / 3E-16
                WRITE(*,*) 'CONDITION NUMBER FOR N = ', NI, ' : ', COND
                DEALLOCATE(Z, R)
            ENDIF

            DEALLOCATE(U, B)
        ENDDO
    ENDDO

END PROGRAM
