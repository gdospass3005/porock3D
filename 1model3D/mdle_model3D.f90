MODULE mdle_model3D

CONTAINS

SUBROUTINE model3Dsub(Io,panel,d1,d2)

IMPLICIT NONE
  INTEGER, INTENT(IN) :: Io
  INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: panel
  INTEGER, INTENT(IN) :: d1,d2

CALL plan3D(panel,Io,d1,d2)

END SUBROUTINE model3Dsub


SUBROUTINE plan3D(A,Io,a1,a2)
IMPLICIT NONE
  INTEGER, INTENT(IN) :: Io
  INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: A
  INTEGER, INTENT(IN) :: a1,a2

A(:,:,:Io) = a1
A(:,:,(Io +1):) = a2

END SUBROUTINE plan3D

END MODULE mdle_model3D




