PROGRAM model3D

  USE mdle_model3D
  IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! VARIABLES USED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER :: i,j,Nx,Ny,Nz,Io
  INTEGER :: d1,d2

  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: panel

  CHARACTER(LEN=20) :: infile='data.dat'
  CHARACTER(LEN=20) :: outfile1='panel.ad'

  LOGICAL :: verbose = .TRUE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  OPEN(25,FILE=infile,STATUS='UNKNOWN')

  READ(25,*) Nx, Ny, Nz
  READ(25,*) Io
  READ(25,*) d1,d2



  PRINT*,'PROGRAMA GERADOR DE PAINEIS PARA O MODELAMENTO POROELASTICO'

  IF (verbose .EQV. .TRUE.) THEN
     PRINT*, Nx,Ny,Nz,' = grid size in x, y and z'
     PRINT*, Io,' = depth of the plane interface in grid units'
     PRINT*, d1,d2,' = layers'
  END IF

  ALLOCATE(panel(Nx,Ny,Nz))

  CALL model3Dsub(Io,panel,d1,d2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  OPEN(35,FILE=outfile1,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nx*4)
  DO J=1,Nx  
  DO I=1,Nx  
     WRITE(35, REC=i+(j-1)*Nx) panel(i,j,:) 
  END DO
  END DO

END PROGRAM model3D
