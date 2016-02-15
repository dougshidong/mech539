      PROGRAM A2MAIN

      USE PREC
      IMPLICIT NONE

      INTEGER, PARAMETER        :: NI = 11, NJ = NI, NIJ = NI * NJ
      INTEGER                   :: I, J, IJ
      REAL*8, DIMENSION (NIJ) :: U, X, Y
      REAL*8                  :: DX, DY
      REAL*8                  :: TOL

      DX = 1.0 / (NI - 1)
      DY = 1.0 / (NJ - 1)
 
      DO I = 1, NI
        DO J = 1, NJ
            IJ = (I - 1) * NI + J
            X(IJ) = (I - 1) * DX
            Y(IJ) = (J - 1) * DY
            U(IJ) = 0
            
        ENDDO
      ENDDO

      WRITE(*,*) 'X'
      WRITE(*,*) X
      WRITE(*,*) 'Y'
      WRITE(*,*) Y
      WRITE(*,*) 'U'
      WRITE(*,*) U
      WRITE(*,*) 'DX'
      WRITE(*,*) DX
      WRITE(*,*) 'DY'
      WRITE(*,*) DY
      WRITE(*,*) 1.0 / 10.0

      END PROGRAM

!***********************************************************************
