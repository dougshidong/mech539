PROGRAM A2MAIN
    USE SOLVER
    
    IMPLICIT REAL (A-H,O-Z)
    
    REAL, ALLOCATABLE            :: U(:), B(:), Z(:), R(:), WEIGHTS(:)
    REAL, DIMENSION (10000000)   :: RESI, TIMES
    CHARACTER(LEN = 11)          :: FNAME
    INTEGER                      :: ISIZE(3)

    NW = 7
    ALLOCATE(WEIGHTS(NW))
    WEIGHTS(:) = (/0.75, 1.00, 1.25, 1.5, 1.75, 1.90, 1.99/)
    ISIZE(1:3) = (/ 100, 200, 400 /)
!    ISIZE(1:5) = (/ 50, 100, 200, 300, 400/)
    ICOND = 0 ! CHECK CONDITION NUMBER
    DO IW = 1, 1
    W = WEIGHTS(IW)
    W = 1.5
    DO ISOLV = 1, 3
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
            TOL = 1E-6
            CALL SOLVELAPLACE(U, B, RESI, TIMES, ITER, ISOLV, W)

            WRITE(FNAME,'("RESULTS", I1, I3 )' )  ISOLV, NI
            !WRITE(FNAME,'("RESULTS", I1, I3)' )  IW, NI
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
            DO I = 1, ITER
                WRITE(11,*) TIMES(I)
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
                CALL SOLVELAPLACE(Z, R, RESI, TIMES, ITER, JSOLV, W)
                COND = NORM2(Z) / NORM2(U) / 3E-16
                WRITE(*,*) 'CONDITION NUMBER FOR N = ', NI, ' : ', COND
                DEALLOCATE(Z, R)
            ENDIF

            DEALLOCATE(U, B)
        ENDDO
    ENDDO
    ENDDO

END PROGRAM
