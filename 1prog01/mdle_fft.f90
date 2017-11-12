MODULE mdle_fft

INTERFACE
SUBROUTINE ffttrig(N,trig)
real, DIMENSION(2*N), INTENT(OUT) :: trig
INTEGER, INTENT(IN) :: N
END SUBROUTINE ffttrig
SUBROUTINE factork(ifax,nfac,N)
INTEGER, INTENT(IN) :: N
INTEGER, DIMENSION(N), INTENT(OUT) :: ifax
INTEGER, INTENT(OUT) :: nfac
END SUBROUTINE factork
SUBROUTINE fft(data,cccc,N,trig,ifax,nfac,skip,isign)
REAL, DIMENSION(2*N), INTENT(INOUT) :: data
REAL, DIMENSION(2*N), INTENT(INOUT) :: cccc
REAL, DIMENSION(2*N), INTENT(IN) :: trig
INTEGER, INTENT(IN) :: skip
INTEGER, INTENT(IN) :: N,isign
INTEGER, DIMENSION(N), INTENT(IN) :: ifax
INTEGER, INTENT(IN) :: nfac
END SUBROUTINE fft
!SUBROUTINE fourt(data,N,ndim,isign,iform,work)
!INTEGER, INTENT(IN) :: N
!REAL, DIMENSION(2*N), INTENT(INOUT) :: data
!REAL, DIMENSION(2*N), INTENT(INOUT) :: work
!INTEGER, INTENT(IN) :: isign,ndim,iform
!END SUBROUTINE fourt
END INTERFACE

CONTAINS

SUBROUTINE fft_init(N,trig,ifax,nfac,dx,rkx)
INTEGER, INTENT(IN) :: N
REAL, DIMENSION(2*N), INTENT(OUT) :: trig
INTEGER, DIMENSION(N), INTENT(OUT) :: ifax
INTEGER, INTENT(OUT) :: nfac
REAL, INTENT(IN) :: dx
REAL, DIMENSION(N), INTENT(OUT) :: rkx
REAL :: PI=3.141592653589793238462643383279502884197
REAL :: dkx
INTEGER :: inyq_kx

call ffttrig(N,trig)
call factork(ifax,nfac,N)

! Spatial frequency sampling interval
dkx = 2*PI/(N*dx)

! Location of the Nyquist frequency on the frequency vector rkx
inyq_kx = N/2 + 1

! Nyquist spatial frequency
!kx_nyq = inyq_kx * dkx

! Frequency vector
do i=1,inyq_kx
  rkx(i) = (i-1)*dkx
end do
do i=inyq_kx+1,N
  rkx(i) = -rkx(N+2-i)
end do

END SUBROUTINE fft_init


SUBROUTINE ddx_shift(data,N,trig,ifax,nfac,rkx,shift)
IMPLICIT NONE
REAL, DIMENSION(N), INTENT(INOUT) :: data
INTEGER, INTENT(IN) :: N
REAL, DIMENSION(2*N), INTENT(IN) :: trig
INTEGER, DIMENSION(N), INTENT(IN) :: ifax
INTEGER, INTENT(IN) :: nfac
REAL, INTENT(IN) :: shift
REAL, DIMENSION(N), INTENT(IN) :: rkx
REAL, DIMENSION(2*N) :: datadouble
REAL, DIMENSION(2*N) :: cccc
COMPLEX, DIMENSION(N) :: cdata
COMPLEX :: j

datadouble(1:2*N:2) = data
datadouble(2:2*N:2) = 0.0
call fft(datadouble,cccc,N,trig,ifax,nfac,1,1)
cdata=cmplx(datadouble(1:2*N:2),datadouble(2:2*N:2))
j = cmplx(0.,1.)
cdata=rkx*j*cdata*exp(rkx*j*shift)   ! derivacao c/ malha escalonada
datadouble(1:2*N:2)=real(cdata)
datadouble(2:2*N:2)=aimag(cdata)
call fft(datadouble,cccc,N,trig,ifax,nfac,1,-1)

data = datadouble(1:2*N:2)/N

END SUBROUTINE ddx_shift

END MODULE mdle_fft
