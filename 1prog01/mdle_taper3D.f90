MODULE mdle_taper3D

CONTAINS

  SUBROUTINE bt_apply_multiple3D(next_u1,next_w1,next_u2,next_w2,next_u3,next_w3,&
       curr_u1,curr_w1,curr_u2,curr_w2,curr_u3,curr_w3,Nx,Ny,Nz,nb,taper,&
       MaxCP,t,dx,lx,ly,lz)
    IMPLICIT NONE
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: next_u1,next_w1,next_u2,next_w2,next_u3,next_w3
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: curr_u1,curr_w1,curr_u2,curr_w2,curr_u3,curr_w3
    INTEGER, INTENT(IN) :: Nx,Ny,Nz,nb
    REAL, INTENT(IN) :: MaxCP,t,dx,lx,ly,lz
    REAL, DIMENSION(nb), INTENT(IN) :: taper
    INTEGER :: ireach, idist1, idist2, idist3

    ireach = CEILING(t*MaxCP/dx)
    idist1 = MIN((lx - nb),((Nx-nb) - lx))
    idist2 = MIN((ly - nb),((Ny-nb) - ly))
    idist3 = MIN((lz - nb),((Nz-nb) - lz))

    IF ((idist1<=ireach).OR.(idist2<=ireach).OR.(idist3<=ireach)) THEN

       CALL bt_apply3D(next_u1,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(next_w1,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(next_u2,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(next_w2,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(next_u3,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(next_w3,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(curr_u1,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(curr_w1,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(curr_u2,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(curr_w2,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(curr_u3,Nx,Ny,Nz,nb,taper)
       CALL bt_apply3D(curr_w3,Nx,Ny,Nz,nb,taper)

    END IF

  END SUBROUTINE bt_apply_multiple3D


SUBROUTINE bt_exp_create(taper,nb,F)
IMPLICIT NONE
  REAL :: F
  INTEGER :: nb
  REAL, DIMENSION(nb) :: taper
  INTEGER :: i

  DO i=1,nb
    taper(i) = exp( -(F*(REAL(nb) - REAL(i) ))**2 )
!    PRINT*, taper(i)
  END DO

END SUBROUTINE bt_exp_create

SUBROUTINE bt_create(taper,nb,cb)
IMPLICIT NONE
  REAL :: cb
  INTEGER :: nb
  REAL, DIMENSION(nb) :: taper
  INTEGER :: i

  DO i=1,nb
    taper(i) = 1 - (REAL(i)/REAL(nb))*(1.0 - cb)
!    PRINT*, taper(i)
  END DO

END SUBROUTINE bt_create

SUBROUTINE bt_apply3D(pp,Nx,Ny,Nz,nb,taper)
  IMPLICIT NONE
  INTEGER :: NX,NY,NZ,nb
  REAL, DIMENSION(nb) :: taper
  REAL, DIMENSION(NX,NY,NZ) :: pp
  INTEGER :: i,i2

  DO i2=1,NY
     DO i=1,NZ
        pp(1:nb,i2,i) = pp(1:nb,i2,i) * taper
        pp(NX:NX-nb+1:-1,i2,i) =  pp(NX:NX-nb+1:-1,i2,i) * taper
     END DO
  END DO

  DO i2=1,NY
     DO i=1,NX
        pp(i,i2,1:nb) = pp(i,i2,1:nb) * taper
        pp(i,i2,NZ:NZ-nb+1:-1) =  pp(i,i2,NZ:NZ-nb+1:-1) * taper  
     END DO
  END DO

  DO i2=1,NZ
     DO i=1,NX
        pp(i,1:nb,i2) = pp(i,1:nb,i2) * taper
        pp(i,NZ:NZ-nb+1:-1,i2) =  pp(i,NZ:NZ-nb+1:-1,i2) * taper  
     END DO
  END DO


END SUBROUTINE bt_apply3D

END MODULE mdle_taper3D



