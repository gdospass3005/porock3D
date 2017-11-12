MODULE mdle_source3D

  IMPLICIT NONE

CONTAINS

  SUBROUTINE source_dervgauss_init(trs,dt,ns,fpeak,tdelay)
    IMPLICIT NONE
    REAL, INTENT(INOUT), DIMENSION(ns) :: trs
    INTEGER, INTENT(IN) :: ns
    REAL, INTENT(IN) :: dt, fpeak
    REAL, INTENT(OUT) :: tdelay
    INTEGER :: i
    REAL :: t,pi=3.1415926536

    tdelay = 4.0 /(3.*fpeak*SQRT(2.*pi))

    do i=1,ns
       t=(i-1)*dt - tdelay
       trs(i) = dervgauss(t,fpeak,1.)
    end do

    OPEN(81,FILE='sourcememory.ad',STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
    WRITE(81, REC=1) trs(:)
    CLOSE(81)

  END SUBROUTINE source_dervgauss_init



  SUBROUTINE source_init(trs,dt,ns,fpeak,tdelay)
    USE mdle_fft
    IMPLICIT NONE
    REAL, INTENT(INOUT), DIMENSION(ns) :: trs
    INTEGER, INTENT(IN) :: ns
    REAL, INTENT(IN) :: dt, fpeak
    REAL, INTENT(OUT) :: tdelay
    INTEGER :: i
    REAL :: t
    REAL, DIMENSION(:), ALLOCATABLE  :: trig
    INTEGER, DIMENSION(:), ALLOCATABLE  :: ifax
    INTEGER  :: nfac
    REAL, DIMENSION(:), ALLOCATABLE  :: rkt
    REAL, DIMENSION(2*Ns) :: datadouble
    REAL, DIMENSION(2*Ns) :: cccc
    COMPLEX, DIMENSION(Ns) :: cdata
    COMPLEX :: j

    tdelay = 1.0 /fpeak

    do i=1,ns
       t=(i-1)*dt - tdelay
       trs(i) = ricker(t,fpeak,1.)
    end do

    ! apply half derivative filter

    ALLOCATE(trig(2*Ns),ifax(Ns),rkt(Ns))

    call fft_init(ns,trig,ifax,nfac,dt,rkt)

    datadouble(1:2*Ns:2) = trs
    datadouble(2:2*Ns:2) = 0.0
    call fft(datadouble,cccc,Ns,trig,ifax,nfac,1,1)
    cdata=cmplx(datadouble(1:2*Ns:2),datadouble(2:2*Ns:2))
    j = cmplx(0.,1.)
    cdata=sqrt(rkt*j)*cdata  ! half derivative
    datadouble(1:2*Ns:2)=real(cdata)
    datadouble(2:2*Ns:2)=aimag(cdata)
    call fft(datadouble,cccc,Ns,trig,ifax,nfac,1,-1)
    trs = datadouble(1:2*Ns:2)/Ns

    trs((CEILING(2.*tdelay/dt)):)=0.

    OPEN(81,FILE='sourcememory.ad',STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
    WRITE(81, REC=1) trs(:)
    CLOSE(81)

    DEALLOCATE(trig,ifax,rkt)

  END SUBROUTINE source_init


  SUBROUTINE source_UU1_3D(t,it,trs,s,r,tdelay)
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, tdelay
    INTEGER, INTENT(IN) :: it
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT), DIMENSION(:,:,:) :: s
    REAL, INTENT(IN), DIMENSION(:) :: trs
    REAL :: stemp,lx,ly,lz
    INTEGER :: ix,iy,iz,ir,ixlow,ixhigh,iylow,iyhigh,izlow,izhigh
    REAL :: ttot,tmin
    REAL :: PI=3.1415926536,r2,d2,wd,aux,dist

    s=0.

    IF (t <= (2.0*tdelay)) THEN
       ir = CEILING(r)
       lx = ir + 0.5
       ly = ir + 0.5
       lz = ir + 0.5
       stemp = trs(it)

       IF (r /= 0) THEN
          ixlow =  1
          ixhigh = 2*ir + 1
          iylow =  1
          iyhigh = 2*ir + 1
          izlow =  1
          izhigh = 2*ir + 1
          r2 = (REAL(r))**2
          aux = LOG(0.25) / r2 
          DO iy=iylow,iyhigh
             DO ix=ixlow,ixhigh
                DO iz=izlow,izhigh
                   d2 = (( REAL(lx-ix) )**2 + ( REAL(ly-iy) )**2 + ( REAL(lz-iz) )**2)
                   dist=sqrt(d2)
                   IF (d2 <= r2) THEN
                      ! wd = EXP( d2 * aux  )
                      wd = (r - dist)/r
                      IF (dist /= 0) THEN
                         s(ix,iy,iz) = ((ix-lx)/dist)*stemp*wd
                      ELSE
                         s(ix,iy,iz) = 0
                      END IF
                   END IF
                END DO
             END DO
          END DO
       ELSE
          STOP 'source radius must be greater or equal to 1.'
          !          s(ABS(lx),ABS(ly),ABS(lz)) = stemp
       END IF
    END IF

  END SUBROUTINE source_UU1_3D


  SUBROUTINE source_UU2_3D(t,it,trs,s,r,tdelay)
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, tdelay
    INTEGER, INTENT(IN) :: it
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT), DIMENSION(:,:,:) :: s
    REAL, INTENT(IN), DIMENSION(:) :: trs
    REAL :: stemp,lx,ly,lz
    INTEGER :: ix,iy,iz,ir,ixlow,ixhigh,iylow,iyhigh,izlow,izhigh
    REAL :: ttot,tmin
    REAL :: PI=3.1415926536,r2,d2,wd,aux,dist

    s=0.

    IF (t <= (2.0*tdelay)) THEN
       ir = CEILING(r)
       lx = ir + 0.5
       ly = ir + 0.5
       lz = ir + 0.5
       stemp = trs(it)
       IF (r /= 0) THEN
          ixlow =  1
          ixhigh = 2*ir + 1
          iylow =  1
          iyhigh = 2*ir + 1
          izlow =  1
          izhigh = 2*ir + 1
          r2 = (REAL(r))**2
          aux = LOG(0.25) / r2 
          DO iy=iylow,iyhigh
             DO ix=ixlow,ixhigh
                DO iz=izlow,izhigh
                   d2 = (( REAL(lx-ix) )**2 + ( REAL(ly-iy) )**2 + ( REAL(lz-iz) )**2)
                   dist=sqrt(d2)
                   IF (d2 <= r2) THEN
                      ! wd = EXP( d2 * aux  )
                      wd = (r - dist)/r
                      IF (dist /= 0) THEN
                         s(ix,iy,iz) = ((iy-ly)/dist)*stemp * wd
                      ELSE
                         s(ix,iy,iz) = 0
                      END IF
                   END IF
                END DO
             END DO
          END DO
       ELSE
          STOP 'source radius must be greater or equal to 1.'
          !          s(ABS(lx),ABS(ly),ABS(lz)) = stemp
       END IF
    END IF

  END SUBROUTINE source_UU2_3D


  SUBROUTINE source_UU3_3D(t,it,trs,s,r,tdelay)
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, tdelay
    INTEGER, INTENT(IN) :: it
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT), DIMENSION(:,:,:) :: s
    REAL, INTENT(IN), DIMENSION(:) :: trs
    REAL :: stemp,lx,ly,lz
    INTEGER :: ix,iy,iz,ir,ixlow,ixhigh,iylow,iyhigh,izlow,izhigh
    REAL :: ttot,tmin
    REAL :: PI=3.1415926536,r2,d2,wd,aux,dist

    s=0.

    IF (t <= (2.0*tdelay)) THEN
       ir = CEILING(r)
       lx = ir + 0.5
       ly = ir + 0.5
       lz = ir + 0.5
       stemp = trs(it)
       IF (r /= 0) THEN
          ixlow =  1
          ixhigh = 2*ir + 1
          iylow =  1
          iyhigh = 2*ir + 1
          izlow =  1
          izhigh = 2*ir + 1
          r2 = (REAL(r))**2
          aux = LOG(0.25) / r2 
          DO iy=iylow,iyhigh
             DO ix=ixlow,ixhigh
                DO iz=izlow,izhigh
                   d2 = (( REAL(lx-ix) )**2 + ( REAL(ly-iy) )**2 + ( REAL(lz-iz) )**2)
                   dist=sqrt(d2)
                   IF (d2 <= r2) THEN
                      ! wd = EXP( d2 * aux  )
                      wd = (r - dist)/r
                      IF (dist /= 0) THEN
                         s(ix,iy,iz) = ((iz-lz)/dist)*stemp * wd
                      ELSE
                         s(ix,iy,iz) = 0
                      END IF
                   END IF
                END DO
             END DO
          END DO
       ELSE
          STOP 'source radius must be greater or equal to 1.'
          !          s(ABS(lx),ABS(ly),ABS(lz)) = stemp
       END IF
    END IF

  END SUBROUTINE source_UU3_3D



  REAL FUNCTION ricker(t,fpeak,bombweight)
    REAL, INTENT(IN) :: t,fpeak
    REAL :: x,xx,pi=3.1415926536
    REAL, INTENT(IN) :: bombweight

    x=pi*fpeak*t
    xx=x*x

    ricker= bombweight*EXP(-xx)*(1.0 - 2.0*xx)

  END FUNCTION ricker



  REAL FUNCTION dervgauss(t,fpeak,bombweight)
    REAL, INTENT(IN) :: t,fpeak
    REAL :: x,xx,pi=3.1415926536
    REAL, INTENT(IN) :: bombweight

    x=3*fpeak*t
    xx=x*x

    dervgauss= bombweight*(-2.*pi*t)*EXP(-pi*xx)

  END FUNCTION dervgauss



END MODULE mdle_source3D
