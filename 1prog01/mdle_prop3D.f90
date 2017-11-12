MODULE mdle_prop3D

IMPLICIT NONE

CONTAINS

  SUBROUTINE prop_porous3D(next_u1,next_w1,next_u2,next_w2,next_u3,next_w3,&
       curr_u1,curr_w1,curr_u2,curr_w2,curr_u3,curr_w3,&
       a,b,c,d,e,f,g,h,panel,&
       dx,dt,Nx,Ny,Nz,MaxCP,t,lx,ly,lz,ir,cells) 
    IMPLICIT NONE
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: next_u1,next_w1,next_u2,next_w2,next_u3,next_w3
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: curr_u1,curr_w1,curr_u2,curr_w2,curr_u3,curr_w3
    REAL, DIMENSION(:), INTENT(INOUT) :: a,b,c,d,e,f,g,h
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: panel
    REAL, INTENT(IN) :: dx,dt,MaxCP,t,lx,ly,lz
    INTEGER, INTENT(IN) :: Nx,Ny,Nz,ir
    INTEGER, INTENT(OUT) :: cells
    INTEGER :: i,j,iy,ireach,idist
    REAL :: ux_xx,ux_yy,ux_zz,ux_xy,ux_xz,ux_yz
    REAL :: uy_xx,uy_yy,uy_zz,uy_xy,uy_xz,uy_yz
    REAL :: uz_xx,uz_yy,uz_zz,uz_xy,uz_xz,uz_yz
    REAL :: wx_xx,wx_yy,wx_zz,wx_xy,wx_xz,wx_yz
    REAL :: wy_xx,wy_yy,wy_zz,wy_xy,wy_xz,wy_yz
    REAL :: wz_xx,wz_yy,wz_zz,wz_xy,wz_xz,wz_yz
    REAL :: z1,z2,z3,z4,z5,z6,z7,z8,z9,ah,bh,ch,dh,eh,fh,gh,hh
    REAL :: dx2,dt2,ddt,aux1,aux3,aux4,aux5

    ! Obs: next_UU has retained the value of UU(t-1)

    dx2=dx*dx
    dt2=dt*dt
    ddt=2.*dt

    ireach = CEILING((t*MaxCP/dx)+ir)
    cells = 0

    DO iy=2,Ny-1
       DO i=2,Nx-1
          DO j=2,Nz-1

             idist = FLOOR(SQRT( (REAL(lx-i))**2 +  (REAL(ly-iy))**2 + (REAL(lz-j))**2 ))

             IF (idist <= ireach) THEN
	        cells=cells+1
                ux_xx = (curr_u1(i+1,iy,j) - 2.*curr_u1(i,iy,j) + curr_u1(i-1,iy,j))/dx2
                ux_yy = (curr_u1(i,iy+1,j) - 2.*curr_u1(i,iy,j) + curr_u1(i,iy-1,j))/dx2
                ux_zz = (curr_u1(i,iy,j+1) - 2.*curr_u1(i,iy,j) + curr_u1(i,iy,j-1))/dx2
                ux_xy = (curr_u1(i+1,iy+1,j) - curr_u1(i-1,iy+1,j) - curr_u1(i+1,iy-1,j) + &
                     curr_u1(i-1,iy-1,j))/(4.*dx2)
                ux_xz = (curr_u1(i+1,iy,j+1) - curr_u1(i-1,iy,j+1) - curr_u1(i+1,iy,j-1) + &
                     curr_u1(i-1,iy,j-1))/(4.*dx2)
                ux_yz = (curr_u1(i,iy+1,j+1) - curr_u1(i,iy-1,j+1) - curr_u1(i,iy+1,j-1) + &
                     curr_u1(i,iy-1,j-1))/(4.*dx2)
                uy_xx = (curr_u2(i+1,iy,j) - 2.*curr_u2(i,iy,j) + curr_u2(i-1,iy,j))/dx2
                uy_yy = (curr_u2(i,iy+1,j) - 2.*curr_u2(i,iy,j) + curr_u2(i,iy-1,j))/dx2
                uy_zz = (curr_u2(i,iy,j+1) - 2.*curr_u2(i,iy,j) + curr_u2(i,iy,j-1))/dx2
                uy_xy = (curr_u2(i+1,iy+1,j) - curr_u2(i-1,iy+1,j) - curr_u2(i+1,iy-1,j) + &
                     curr_u2(i-1,iy-1,j))/(4.*dx2)
                uy_xz = (curr_u2(i+1,iy,j+1) - curr_u2(i-1,iy,j+1) - curr_u2(i+1,iy,j-1) + &
                     curr_u2(i-1,iy,j-1))/(4.*dx2)
                uy_yz = (curr_u2(i,iy+1,j+1) - curr_u2(i,iy-1,j+1) - curr_u2(i,iy+1,j-1) + &
                     curr_u2(i,iy-1,j-1))/(4.*dx2)
                uz_xx = (curr_u3(i+1,iy,j) - 2.*curr_u3(i,iy,j) + curr_u3(i-1,iy,j))/dx2
                uz_yy = (curr_u3(i,iy+1,j) - 2.*curr_u3(i,iy,j) + curr_u3(i,iy-1,j))/dx2
                uz_zz = (curr_u3(i,iy,j+1) - 2.*curr_u3(i,iy,j) + curr_u3(i,iy,j-1))/dx2
                uz_xy = (curr_u3(i+1,iy+1,j) - curr_u3(i-1,iy+1,j) - curr_u3(i+1,iy-1,j) + &
                     curr_u3(i-1,iy-1,j))/(4.*dx2)
                uz_xz = (curr_u3(i+1,iy,j+1) - curr_u3(i-1,iy,j+1) - curr_u3(i+1,iy,j-1) + &
                     curr_u3(i-1,iy,j-1))/(4.*dx2)
                uz_yz = (curr_u3(i,iy+1,j+1) - curr_u3(i,iy-1,j+1) - curr_u3(i,iy+1,j-1) + &
                     curr_u3(i,iy-1,j-1))/(4.*dx2)

                wx_xx = (curr_w1(i+1,iy,j) - 2.*curr_w1(i,iy,j) + curr_w1(i-1,iy,j))/dx2
                wx_yy = (curr_w1(i,iy+1,j) - 2.*curr_w1(i,iy,j) + curr_w1(i,iy-1,j))/dx2
                wx_zz = (curr_w1(i,iy,j+1) - 2.*curr_w1(i,iy,j) + curr_w1(i,iy,j-1))/dx2
                wx_xy = (curr_w1(i+1,iy+1,j) - curr_w1(i-1,iy+1,j) - curr_w1(i+1,iy-1,j) + &
                     curr_w1(i-1,iy-1,j))/(4.*dx2)
                wx_xz = (curr_w1(i+1,iy,j+1) - curr_w1(i-1,iy,j+1) - curr_w1(i+1,iy,j-1) + &
                     curr_w1(i-1,iy,j-1))/(4.*dx2)
                wx_yz = (curr_w1(i,iy+1,j+1) - curr_w1(i,iy-1,j+1) - curr_w1(i,iy+1,j-1) + &
                     curr_w1(i,iy-1,j-1))/(4.*dx2)
                wy_xx = (curr_w2(i+1,iy,j) - 2.*curr_w2(i,iy,j) + curr_w2(i-1,iy,j))/dx2
                wy_yy = (curr_w2(i,iy+1,j) - 2.*curr_w2(i,iy,j) + curr_w2(i,iy-1,j))/dx2
                wy_zz = (curr_w2(i,iy,j+1) - 2.*curr_w2(i,iy,j) + curr_w2(i,iy,j-1))/dx2
                wy_xy = (curr_w2(i+1,iy+1,j) - curr_w2(i-1,iy+1,j) - curr_w2(i+1,iy-1,j) + &
                     curr_w2(i-1,iy-1,j))/(4.*dx2)
                wy_xz = (curr_w2(i+1,iy,j+1) - curr_w2(i-1,iy,j+1) - curr_w2(i+1,iy,j-1) + &
                     curr_w2(i-1,iy,j-1))/(4.*dx2)
                wy_yz = (curr_w2(i,iy+1,j+1) - curr_w2(i,iy-1,j+1) - curr_w2(i,iy+1,j-1) + &
                     curr_w2(i,iy-1,j-1))/(4.*dx2)
                wz_xx = (curr_w3(i+1,iy,j) - 2.*curr_w3(i,iy,j) + curr_w3(i-1,iy,j))/dx2
                wz_yy = (curr_w3(i,iy+1,j) - 2.*curr_w3(i,iy,j) + curr_w3(i,iy-1,j))/dx2
                wz_zz = (curr_w3(i,iy,j+1) - 2.*curr_w3(i,iy,j) + curr_w3(i,iy,j-1))/dx2
                wz_xy = (curr_w3(i+1,iy+1,j) - curr_w3(i-1,iy+1,j) - curr_w3(i+1,iy-1,j) + &
                     curr_w3(i-1,iy-1,j))/(4.*dx2)
                wz_xz = (curr_w3(i+1,iy,j+1) - curr_w3(i-1,iy,j+1) - curr_w3(i+1,iy,j-1) + &
                     curr_w3(i-1,iy,j-1))/(4.*dx2)
                wz_yz = (curr_w3(i,iy+1,j+1) - curr_w3(i,iy-1,j+1) - curr_w3(i,iy+1,j-1) + &
                     curr_w3(i,iy-1,j-1))/(4.*dx2)

                z1 = ux_xx + uy_xy + uz_xz 
                z2 = wx_xx + wy_xy + wz_xz
                z3 = -(uy_xy - ux_yy - (ux_zz - uz_xz))
                z4 = ux_xy + uy_yy + uz_yz
                z5 = wx_xy + wy_yy + wz_yz
                z6 = -(uz_yz - uy_zz - (uy_xx - ux_xy))
                z7 = ux_xz + uy_yz + uz_zz
                z8 = wx_xz + wy_yz + wz_zz
                z9 = -(ux_xz - uz_xx - (uz_yy - uy_yz))

                ah=a(panel(i,iy,j))
                bh=b(panel(i,iy,j))
                ch=c(panel(i,iy,j))
                dh=d(panel(i,iy,j))
                eh=e(panel(i,iy,j))
                fh=f(panel(i,iy,j))
                gh=g(panel(i,iy,j))
                hh=h(panel(i,iy,j))

                aux5=dh*dt/2.
                aux4=hh/ddt
                aux3=1./dt2 + aux4

                aux1=next_w1(i,iy,j)
                next_w1(i,iy,j)=(eh*z1+fh*z2-gh*z3+aux4*aux1 + &
                     (2.*curr_w1(i,iy,j)-aux1)/dt2)/aux3
                next_u1(i,iy,j)=(ah*z1+bh*z2+ch*z3)*dt2+aux5*(next_w1(i,iy,j)-aux1) + &
                     2.*curr_u1(i,iy,j)-next_u1(i,iy,j)

                aux1=next_w2(i,iy,j)
                next_w2(i,iy,j)=(eh*z4+fh*z5-gh*z6+aux4*aux1 + &
                     (2.*curr_w2(i,iy,j)-aux1)/dt2)/aux3
                next_u2(i,iy,j)=(ah*z4+bh*z5+ch*z6)*dt2+aux5*(next_w2(i,iy,j)-aux1) + &
                     2.*curr_u2(i,iy,j)-next_u2(i,iy,j)

                aux1=next_w3(i,iy,j)
                next_w3(i,iy,j)=(eh*z7+fh*z8-gh*z9+aux4*aux1 + &
                     (2.*curr_w3(i,iy,j)-aux1)/dt2)/aux3
                next_u3(i,iy,j)=(ah*z7+bh*z8+ch*z9)*dt2+aux5*(next_w3(i,iy,j)-aux1) + &
                     2.*curr_u3(i,iy,j)-next_u3(i,iy,j)

             END IF

          END DO
       END DO
    END DO

  END SUBROUTINE prop_porous3D

END MODULE mdle_prop3D
