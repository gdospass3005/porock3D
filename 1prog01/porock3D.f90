PROGRAM porock3D

  !Author:                 Gino Francisco dos Passos
  !Last edit:              SSA, February 1st., 2001.
  !Status:                 Fully running.

  !Description:            Forward modeling of displacement wavefields
  !                        in a 3-D porous medium due to a P source.
  !                        Displacements are computed with the use of
  !                        the Biot (1962) homogeneous equation for
  !                        isotropic poroelastic media. The methodology
  !                        is based on the work of Zhu & McMechan (1991).
  !                        Same spatial sampling in all axis.

  USE mdle_io_utils3D
  USE mdle_source3D
  USE mdle_utils
  USE mdle_taper3D
  USE mdle_prop3D

  IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! VARIABLES DEFINITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER(LEN=20) :: infile='data.dat'
  CHARACTER(LEN=20) :: parfile='par.dat'
  INTEGER :: Nx, Ny, Nz       ! grid shape (INPUT)
  INTEGER :: nlayers      ! number or layers
  REAL ::    dx           ! space sample interval (INPUT)

  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: panel  ! Model with layer numbers
  REAL, DIMENSION(:), ALLOCATABLE :: cp     ! P wavespeed of dry rock [m/s]
  REAL, DIMENSION(:), ALLOCATABLE :: cs     ! S wavespeed of dry rock [m/s]
  REAL, DIMENSION(:), ALLOCATABLE :: rhos   ! density of rock grain [g/cm^3]
  REAL, DIMENSION(:), ALLOCATABLE :: phi    ! porosity [adim.]
  REAL, DIMENSION(:), ALLOCATABLE :: kap    ! permeability [md]
  REAL, DIMENSION(:), ALLOCATABLE :: aa     ! tortuosity [adim.] !calculated
  REAL, DIMENSION(:), ALLOCATABLE :: cpf    ! P wavespeed of the fluid [m/s]
  REAL, DIMENSION(:), ALLOCATABLE :: rhof   ! density of the fluid [g/cm^3]
  REAL, DIMENSION(:), ALLOCATABLE :: eta    ! viscosity of the fluid [cP]

  CHARACTER(LEN=20) :: infile1='panel.ad'

  REAL :: tmax            ! maximum time to compute per trace (INPUT)
  REAL :: dt              ! time sample interval
  REAL :: dtallow         ! minimum time sample interval allowed
  INTEGER :: nt           ! n. of time samples per shot
  REAL :: t               ! time
  INTEGER :: dt_effect    ! n. of time steps per output time step (INPUT)
  INTEGER :: ns           ! n. of samples per output trace
  REAL :: hh              ! minimum spatial sample interval
  REAL :: VPmax, VPmin    ! maximum, minimum P wavespeeds (dry rock)
  REAL :: VSmax, VSmin    ! maximum, minimum S wavespeeds (dry rock)

  ! Variables related to the source
  REAL :: GG     ! points per wavelenght for the dispersion criterium (INPUT)
  REAL :: fmax            ! maximum temporal frequency allowable
  REAL :: fpeak           ! peak frequency of the Ricker wavelet
  REAL :: fpeakallow      ! maximum peak frequency allowed to the wavelet
  REAL :: tdelay          ! total time of the source to account
  REAL :: lx, ly, lz      ! location of the current source in grid units
  REAL :: r               ! radius of the source (in grid units) (INPUT)
  REAL, DIMENSION(:,:,:), ALLOCATABLE  :: s_panel
  REAL, DIMENSION(:), ALLOCATABLE :: src      ! source memory function 

  ! Variables related to the boundary tapers
  REAL, DIMENSION(:), ALLOCATABLE :: taper
  INTEGER :: nb           ! n. of boundary points (INPUT)
  REAL :: FF              ! absorption factor (INPUT)

  ! Variables related to the survey geometry (INPUT)
  CHARACTER(LEN=20) :: geometryfile='geometry.dat'
  INTEGER :: nshots       ! n. of shots
  INTEGER :: nsmin_x      ! x position of 1st shot in grid units
  INTEGER :: nsmin_y      ! y position of 1st shot in grid units
  INTEGER :: nsmin_z      ! z position of 1st shot in grid units
  INTEGER :: dns_x        ! x interval between two shots in grid units
  INTEGER :: dns_y        ! y interval between two shots in grid units
  INTEGER :: dns_z        ! z interval between two shots in grid units
  INTEGER :: nphones      ! n. of geophones (x and z displacements)
  INTEGER :: npmin_x      ! x position of 1st geophone in grid units
  INTEGER :: npmin_y      ! y position of 1st geophone in grid units
  INTEGER :: npmin_z      ! z position of 1st geophone in grid units
  INTEGER :: dnp_x        ! x interval between two geophones in grid units
  INTEGER :: dnp_y        ! y interval between two geophones in grid units
  INTEGER :: dnp_z        ! z interval between two geophones in grid units
  INTEGER :: shot_to_save       ! shot to take snapshots
  INTEGER :: nsnaps       ! number of snapshots
  INTEGER :: snapmin      ! 1st. snapshot (in time samples)
  INTEGER :: dsnap        ! time samples step between snapshots
  INTEGER :: yslice       ! slice in y to take snapshots

  ! Biot Equations' parameters 
  REAL, DIMENSION(:), ALLOCATABLE :: rho,rhor,rhom,vp2,vs2
  REAL, DIMENSION(:), ALLOCATABLE :: rks,rkf,rk0,rkc,beta,eme

  ! Poroelastic displacement fields 
  REAL, DIMENSION(:,:,:), POINTER :: next_u1,curr_u1,swap1
  REAL, DIMENSION(:,:,:), POINTER :: next_w1,curr_w1,swap2
  REAL, DIMENSION(:,:,:), POINTER :: next_u2,curr_u2,swap3
  REAL, DIMENSION(:,:,:), POINTER :: next_w2,curr_w2,swap4
  REAL, DIMENSION(:,:,:), POINTER :: next_u3,curr_u3,swap5
  REAL, DIMENSION(:,:,:), POINTER :: next_w3,curr_w3,swap6

  ! Maximum and minimum wavespeeds
  REAL :: MaxCP, MinCS

  ! Auxiliary variables
  REAL, DIMENSION(:), ALLOCATABLE :: gama,be,a,b,c,d,e,f,g,h
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: wr,ws,wf

  INTEGER :: i,iy,j,ishot,it,is,icsg,isnap,ip,iyp,jp     ! counters
  INTEGER :: aux,istat,ir,cells

  ! Auxiliary panels for output common shot gathers
  REAL, DIMENSION(:,:), ALLOCATABLE :: csg_u1,csg_w1,csg_u2,csg_w2,csg_u3,csg_w3 

  ! Output files
  CHARACTER(LEN=20) :: outfile1='csg_u1.ad'
  CHARACTER(LEN=20) :: outfile2='csg_w1.ad'
  CHARACTER(LEN=20) :: outfile3='csg_u2.ad'
  CHARACTER(LEN=20) :: outfile4='csg_w2.ad'
  CHARACTER(LEN=20) :: outfile5='csg_u3.ad'
  CHARACTER(LEN=20) :: outfile6='csg_w3.ad'
  CHARACTER(LEN=20) :: snapfile1='snapshot_u1.ad'
  CHARACTER(LEN=20) :: snapfile2='snapshot_w1.ad'
  CHARACTER(LEN=20) :: snapfile3='snapshot_u2.ad'
  CHARACTER(LEN=20) :: snapfile4='snapshot_w2.ad'
  CHARACTER(LEN=20) :: snapfile5='snapshot_u3.ad'
  CHARACTER(LEN=20) :: snapfile6='snapshot_w3.ad'

  ! Spent time handle
  INTEGER :: tmstart,tmcurr,tmlast,tmdt,tmest,secnds
  REAL :: t_per_cell
  CHARACTER(LEN=20) :: timefile='itertime.txt'
  CHARACTER(LEN=20) :: timefile2='t_per_cell.txt'

  ! Report File
  CHARACTER(LEN=20) :: reportfile='report.txt'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Uploading general parameters
  CALL input1_3D(infile,Nx,Ny,Nz,nlayers,dx,tmax,dt_effect,GG,r,nb,FF,&
       dt,fpeak) 

  ! Uploading acquisition geometry parameters
  CALL input2_3D(geometryfile,nshots,nsmin_x,nsmin_y,nsmin_z,dns_x,dns_y,dns_z,&
       nphones,npmin_x,npmin_y,npmin_z,dnp_x,dnp_y,dnp_z,shot_to_save,nsnaps,&
       snapmin,dsnap,yslice)

  ! Uploading rock and fluid parameters
  ALLOCATE(cp(nlayers),cs(nlayers),rhos(nlayers),phi(nlayers),&
       kap(nlayers),aa(nlayers),cpf(nlayers),rhof(nlayers),eta(nlayers))
  CALL input3_3D(parfile,nlayers,cp,cs,rhos,phi,kap,cpf,rhof,eta)

  ! Upload position of layers (3D model)
  ALLOCATE(panel(Nx,Ny,Nz))
  CALL startspanel3D(panel,  Nx,Ny,Nz,infile1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Get Biot parameters from rock and fluid parameters
  phi=phi/100.
  kap=kap*1.e-15
  rhos=rhos*1000.
  rhof=rhof*1000.
  eta=eta/1000.
  aa=1.-0.5*(1. - 1./phi)

  ALLOCATE(rho(nlayers),rhor(nlayers),rhom(nlayers),vp2(nlayers),vs2(nlayers),&
       rks(nlayers),rkf(nlayers),rkc(nlayers),rk0(nlayers),beta(nlayers),&
       eme(nlayers),gama(nlayers),be(nlayers))

  rho = phi*rhof + (1.-phi)*rhos  ! medium density
  rhor = rhof/rho
  rhom = aa*rhof/phi              ! fluid effective density
  vp2=cp*cp
  vs2=cs*cs
  rks = rhos*(vp2 - (4./3)*vs2)   ! solid bulk modulus
  rkf = rhof*cpf**2               ! fluid bulk modulus
  rk0 = rks*(1. - phi)            ! frame bulk modulus
  rkc = rks*(1. - phi) + rkf*phi  ! bulk modulus, saturated rock
  beta = 1. - rk0/rks
  eme = (phi/rkf + (beta-phi)/rks)**(-1)  ! Biot modulus (1941)
  gama = rhom - rhof*rhor
  be = eta/kap

  OPEN(50,FILE=reportfile,STATUS='UNKNOWN',ACTION='WRITE')
  WRITE(50,*) 'POROCK3D: Poroelastic 3-D Modeling the Biot Homogeneous'
  WRITE(50,*) '          Equations with 2nd order finite differences'
  WRITE(50,*) ' '
  WRITE(50,*) 'Variable -----  maximum / minimum'
  WRITE(50,*) ' phi          ',MAXVAL(phi),MINVAL(phi)
  WRITE(50,*) ' kap          ',MAXVAL(kap),MINVAL(kap)
  WRITE(50,*) ' rhos         ',MAXVAL(rhos),MINVAL(rhos)
  WRITE(50,*) ' rhof         ',MAXVAL(rhof),MINVAL(rhof)
  WRITE(50,*) ' eta          ',MAXVAL(eta),MINVAL(eta)
  WRITE(50,*) ' vp2          ',MAXVAL(vp2),MINVAL(vp2)
  WRITE(50,*) ' vs2          ',MAXVAL(vs2),MINVAL(vs2)
  WRITE(50,*) ' beta         ',MAXVAL(beta),MINVAL(beta)
  WRITE(50,*) ' eme          ',MAXVAL(eme),MINVAL(eme)
  WRITE(50,*) ' rhor         ',MAXVAL(rhor),MINVAL(rhor)
  WRITE(50,*) ' gama         ',MAXVAL(gama),MINVAL(gama)
  WRITE(50,*) ' aa           ',MAXVAL(aa),MINVAL(aa)
  WRITE(50,*) ' rho          ',MAXVAL(rho),MINVAL(rho)
  WRITE(50,*) ' rhom         ',MAXVAL(rhom),MINVAL(rhom)
  WRITE(50,*) ' rks          ',MAXVAL(rks),MINVAL(rks)
  WRITE(50,*) ' rkf          ',MAXVAL(rkf),MINVAL(rkf)
  WRITE(50,*) ' rk0          ',MAXVAL(rk0),MINVAL(rk0)
  WRITE(50,*) ' rkc          ',MAXVAL(rkc),MINVAL(rkc)
  WRITE(50,*) ' be           ',MAXVAL(be),MINVAL(be)

  ! Evaluate maximum and minimum velocities
  MaxCP = MAXVAL(SQRT(cp**2 + cs**2))
  MinCS = MINVAL((cpf/(aa**(1./2.))))  ! Slow-P minimum velocity
  PRINT*,'MaxCp,MinCs  ',MaxCp,MinCs

  DEALLOCATE(rks,rk0,rkf,rkc,kap,aa,cpf,eta,rhos)

  ALLOCATE(a(nlayers),b(nlayers),c(nlayers),d(nlayers),e(nlayers),&
       f(nlayers),g(nlayers),h(nlayers))
  a=(rhom*vp2 - beta*eme*rhor)/gama
  b=eme*(beta*rhom/rho-rhor)/gama
  c=rhom*vs2/gama
  d=rhor*be/gama
  e=(beta*eme-rhof*vp2)/gama
  f=eme*(1. - beta*rhor)/gama
  g=rhof*vs2/gama
  h=be/gama

  WRITE(50,*) ' '
  WRITE(50,*) ' a            ',MAXVAL(a),MINVAL(a)
  WRITE(50,*) ' b            ',MAXVAL(b),MINVAL(b)
  WRITE(50,*) ' c            ',MAXVAL(c),MINVAL(c)
  WRITE(50,*) ' d            ',MAXVAL(d),MINVAL(d)
  WRITE(50,*) ' e            ',MAXVAL(e),MINVAL(e)
  WRITE(50,*) ' f            ',MAXVAL(f),MINVAL(f)
  WRITE(50,*) ' g            ',MAXVAL(g),MINVAL(g)
  WRITE(50,*) ' h            ',MAXVAL(h),MINVAL(h)
  WRITE(50,*) ' ws           ',MAXVAL(ws),MINVAL(ws)
  WRITE(50,*) ' wr           ',MAXVAL(wr),MINVAL(wr)
  WRITE(50,*) ' '

  DEALLOCATE(rho,rhof,rhor,rhom,vp2,vs2,beta,eme,gama,be,&
       cp,cs)

  ir = CEILING(r)
  ALLOCATE(ws(2*ir+1,2*ir+1,2*ir+1),wf(2*ir+1,2*ir+1,2*ir+1),&
       wr(2*ir+1,2*ir+1,2*ir+1))
  DO i=1,2*ir+1
     DO iy=1,2*ir+1
        DO j=1,2*ir+1
           ip = lx+((i -1)-ir)
           iyp= ly+((iy-1)-ir)
           jp = lz+((j -1)-ir)
           ws(i,iy,j)=1. - phi(panel(ip,iyp,jp))
           wf(i,iy,j)=phi(panel(ip,iyp,jp))
           wr(i,iy,j)=phi(panel(ip,iyp,jp)) * ABS(wf(i,iy,j)-ws(i,iy,j))
        END DO
     END DO
  END DO
  DEALLOCATE(wf,phi)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Minimum spatial sampling interval
  hh = dx

  ! Time sampling interval to ensure stability
  dtallow = 0.61*hh/MaxCP
  IF (dt > dtallow) STOP 'dt not allowed'

  ! Total number of time steps
  nt = 1 + FLOOR(tmax/dt)

  ! Total number of time samples
  ns = nt/dt_effect

  ! Boundary taper
  ALLOCATE(taper(nb))
  CALL bt_exp_create(taper,nb,FF)

  ! Source frequency
  fmax = MinCS/(GG*hh)
  fpeakallow = 0.5*fmax  ! better with 0.2*fmax
  !  IF (fpeak > fpeakallow) &
  !       STOP 'fpeak not allowed'
  PRINT*,'fpeak,fpeakallow',fpeak,fpeakallow
  ! Source memory function vector initialization
  CALL quads(nt,aux)
  ALLOCATE(src(aux))
  !CALL source_init(src,dt,aux,fpeak,tdelay)
  CALL source_dervgauss_init(src,dt,aux,fpeak,tdelay)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Last settings before run

  ALLOCATE(next_u1(Nx,Ny,Nz),curr_u1(Nx,Ny,Nz)) 
  ALLOCATE(next_w1(Nx,Ny,Nz),curr_w1(Nx,Ny,Nz))
  ALLOCATE(next_u2(Nx,Ny,Nz),curr_u2(Nx,Ny,Nz))
  ALLOCATE(next_w2(Nx,Ny,Nz),curr_w2(Nx,Ny,Nz))
  ALLOCATE(next_u3(Nx,Ny,Nz),curr_u3(Nx,Ny,Nz))
  ALLOCATE(next_w3(Nx,Ny,Nz),curr_w3(Nx,Ny,Nz))

  ALLOCATE(s_panel(2*ir+1,2*ir+1,2*ir+1))
  ALLOCATE(csg_u1(ns,nphones),csg_w1(ns,nphones),csg_u2(ns,nphones),&
       csg_w2(ns,nphones),csg_u3(ns,nphones),csg_w3(ns,nphones))

  OPEN(51,FILE=outfile1,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
  OPEN(52,FILE=outfile2,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
  OPEN(53,FILE=outfile3,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
  OPEN(54,FILE=outfile4,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
  OPEN(55,FILE=outfile5,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
  OPEN(56,FILE=outfile6,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)

  OPEN(61,FILE=snapfile1,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(62,FILE=snapfile2,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(63,FILE=snapfile3,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(64,FILE=snapfile4,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(65,FILE=snapfile5,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(66,FILE=snapfile6,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

  OPEN(82,FILE=timefile,STATUS='UNKNOWN',ACTION='WRITE')
  OPEN(83,FILE=timefile2,STATUS='UNKNOWN',ACTION='WRITE')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Report
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(50,*) 'Grid shape (Nx,Ny,Nz):             ',Nx,Ny,Nz
  WRITE(50,*) 'Spatial sampling interval (dx):    ',dx
  WRITE(50,*) 'Total time per trace (tmax):       ',tmax
  WRITE(50,*) 'Number of time samples (ns):       ',ns
  WRITE(50,*) 'Output dt (=dt*dt_effect):         ',dt*dt_effect
  WRITE(50,*) 'Time Interval (dt):                ',dt
  WRITE(50,*) 'Max. Time Interval (dtallow):      ',dtallow
  WRITE(50,*) 'Maximum and minimum wavespeeds'
  WRITE(50,*) '(MaxCP, MinCS):                    ',MaxCP,MinCS
  WRITE(50,*) 'Maximum source frequency (fmax):   ',fmax
  WRITE(50,*) 'Peak frequency (fpeak)             ',fpeak
  WRITE(50,*) 'Max. Peak frequency (fpeakallow)   ',fpeakallow
  WRITE(50,*) 'size of source vector, in samples: ',CEILING(2.*tdelay/dt)

  CLOSE(50)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  icsg=0;
  ! Loop over shots
  DO ishot=1,nshots
     lx = dns_x * (ishot-1) + nsmin_x
     ly = dns_y * (ishot-1) + nsmin_y
     lz = dns_z * (ishot-1) + nsmin_z

     IF ((lx<=ir).OR.(lx>(Nx-ir))) STOP 'lx not allowed'
     IF ((ly<=ir).OR.(ly>(Ny-ir))) STOP 'ly not allowed'
     IF ((lz<=ir).OR.(lz>(Nz-ir))) STOP 'lz not allowed'

     next_u1=0; curr_u1=0
     next_w1=0; curr_w1=0
     next_u2=0; curr_u2=0
     next_w2=0; curr_w2=0
     next_u3=0; curr_u3=0
     next_w3=0; curr_w3=0

     ! Loop over time steps
     t=0.; is=0; isnap=0; cells=1
     tmstart=secnds(0.0); tmcurr=0; tmlast=0; tmdt=0; tmest=0
     DO it=1,nt
        t = (it-1)*dt

	tmlast = tmcurr
	tmcurr=secnds(REAL(tmstart)) 
	tmdt = tmcurr-tmlast
        WRITE(82,*)it,tmdt
	tmest = (REAL((nt-it)*tmdt)/cells)*(Nx*Ny*Nz)
        t_per_cell = REAL(tmdt)/cells
        WRITE(83,*)it,t_per_cell

        ! Insert source function
        CALL source_UU1_3D(t,it,src,s_panel,r,tdelay)
        curr_u1(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) = &
             curr_u1(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) + ws*s_panel
        curr_w1(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) = &
             curr_w1(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) + wr*s_panel
        CALL source_UU2_3D(t,it,src,s_panel,r,tdelay)
        curr_u2(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) = &
             curr_u2(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) + ws*s_panel
        curr_w2(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) = &
             curr_w2(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) + wr*s_panel
        CALL source_UU3_3D(t,it,src,s_panel,r,tdelay)
        curr_u3(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) = &
             curr_u3(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) + ws*s_panel    
        curr_w3(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) = &
             curr_w3(lx-ir:lx+ir,ly-ir:ly+ir,lz-ir:lz+ir) + wr*s_panel    

        PRINT*, 's',ishot,' it',it,'/',nt,&
             ' u1()=',curr_u1(Nx/2,Ny/2,Nz/2),' tmdt,tmest,',tmdt,tmest,&
             ' cells',cells,'/',(Nx*Ny*Nz)

        ! Do one time step
        CALL prop_porous3D(next_u1,next_w1,next_u2,next_w2,next_u3,next_w3,&
             curr_u1,curr_w1,curr_u2,curr_w2,curr_u3,curr_w3,&
             a,b,c,d,e,f,g,h,panel,&
             dx,dt,Nx,Ny,Nz,MaxCP,t,lx,ly,lz,ir,cells) 

        ! Apply boundary tapers
        CALL bt_apply_multiple3D(next_u1,next_w1,next_u2,next_w2,next_u3,next_w3,&
             curr_u1,curr_w1,curr_u2,curr_w2,curr_u3,curr_w3,Nx,Ny,Nz,nb,taper,&
             MaxCP,t,dx,lx,ly,lz)

        ! Output 
        CALL save_shotgathers_n_snapshots3D(next_u1,next_w1,next_u2,next_w2,next_u3,&
             next_w3,csg_u1,csg_w1,csg_u2,csg_w2,csg_u3,csg_w3,&
             it,is,icsg,isnap,dt_effect,ns,nphones,npmin_x,npmin_y,npmin_z,&
             dnp_x,dnp_y,dnp_z,51,52,53,54,55,56,ishot,shot_to_save,&
             nsnaps,snapmin,dsnap,61,62,63,64,65,66,Nx,Ny,Nz,yslice)

        ! Swap displacement panels
        swap1 => curr_u1
        curr_u1 => next_u1
        next_u1 => swap1

        swap2 => curr_w1
        curr_w1 => next_w1
        next_w1 => swap2

        swap3 => curr_u2
        curr_u2 => next_u2
        next_u2 => swap3

        swap4 => curr_w2
        curr_w2 => next_w2
        next_w2 => swap4

        swap5 => curr_u3
        curr_u3 => next_u3
        next_u3 => swap5

        swap6 => curr_w3
        curr_w3 => next_w3
        next_w3 => swap6

     END DO

  END DO


END PROGRAM porock3D

















