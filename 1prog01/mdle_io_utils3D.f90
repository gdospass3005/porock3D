MODULE mdle_io_utils3D

CONTAINS

  SUBROUTINE save_shotgathers_n_snapshots3D(next_u1,next_w1,next_u2,next_w2,next_u3,&
       next_w3,csg_u1,csg_w1,csg_u2,csg_w2,csg_u3,csg_w3,&
       it,is,icsg,isnap,dt_effect,ns,nphones,npmin_x,npmin_y,npmin_z,&
       dnp_x,dnp_y,dnp_z,n1,n2,n3,n4,n5,n6,ishot,shot_to_save,&
       nsnaps,snapmin,dsnap,nn1,nn2,nn3,nn4,nn5,nn6,Nx,Ny,Nz,yslice)

    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: csg_u1,csg_w1,csg_u2,csg_w2,csg_u3,csg_w3 
    REAL, DIMENSION(:,:,:), INTENT(IN) :: next_u1,next_w1,next_u2,next_w2,next_u3,next_w3
    INTEGER, INTENT(IN) :: it,dt_effect,ns,nphones
    INTEGER, INTENT(INOUT) :: is,icsg,isnap
    INTEGER, INTENT(IN) :: npmin_x,npmin_y,npmin_z,dnp_x,dnp_y,dnp_z,n1,n2,n3,n4,n5,n6
    INTEGER, INTENT(IN) :: ishot,shot_to_save,nsnaps
    INTEGER, INTENT(IN) :: snapmin,dsnap,nn1,nn2,nn3,nn4,nn5,nn6,Nx,Ny,Nz,yslice
    INTEGER :: i,j,ix,iy,iz

    if (MOD((it-1),dt_effect) == 0) then
       is=is+1

       do i=1,nphones
          ix = (i-1)*dnp_x + npmin_x
          iy = (i-1)*dnp_y + npmin_y
          iz = (i-1)*dnp_z + npmin_z
          csg_u1(is,i) = next_u1(ix,iy,iz)
          csg_w1(is,i) = next_w1(ix,iy,iz)
          csg_u2(is,i) = next_u2(ix,iy,iz)
          csg_w2(is,i) = next_w2(ix,iy,iz)
          csg_u3(is,i) = next_u3(ix,iy,iz)
          csg_w3(is,i) = next_w3(ix,iy,iz)
       end do

       if (is == ns) then  
          PRINT*,'Saving CSG'
          icsg = icsg+1
          do i=1,nphones
             WRITE(n1, REC=((icsg -1)*nphones + i)) csg_u1(:,i)
             WRITE(n2, REC=((icsg -1)*nphones + i)) csg_w1(:,i)
             WRITE(n3, REC=((icsg -1)*nphones + i)) csg_u2(:,i)
             WRITE(n4, REC=((icsg -1)*nphones + i)) csg_w2(:,i)
             WRITE(n5, REC=((icsg -1)*nphones + i)) csg_u3(:,i)
             WRITE(n6, REC=((icsg -1)*nphones + i)) csg_w3(:,i)
          end do
       end if

       if (ishot == shot_to_save) then
          if (is == (isnap * dsnap) + snapmin) then
             isnap=isnap+1
             if (isnap <= nsnaps) then
                PRINT*, 'Saving snapshot ',isnap,'/',nsnaps
                j=yslice
                do i=1,Nx
                   WRITE(nn1, REC=((isnap -1)*Nx + i )) next_u1(i,j,:)
                   WRITE(nn2, REC=((isnap -1)*Nx + i )) next_w1(i,j,:)
                   WRITE(nn3, REC=((isnap -1)*Nx + i )) next_u2(i,j,:)
                   WRITE(nn4, REC=((isnap -1)*Nx + i )) next_w2(i,j,:) 
                   WRITE(nn5, REC=((isnap -1)*Nx + i )) next_u3(i,j,:)
                   WRITE(nn6, REC=((isnap -1)*Nx + i )) next_w3(i,j,:)    
                end do

             end if
          end if
       end if
    end if

  END SUBROUTINE save_shotgathers_n_snapshots3D



  SUBROUTINE input1_3D(infile,Nx,Ny,Nz,nlayers,dx,tmax,dt_effect,G,r,nb,F,dt,fpeak)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: infile
    INTEGER, INTENT(OUT) :: Nx, Ny, Nz, nlayers
    REAL, INTENT(OUT) :: dx     
    REAL, INTENT(OUT) :: tmax         
    INTEGER, INTENT(OUT) :: dt_effect 
    REAL, INTENT(OUT) :: G     
    REAL, INTENT(OUT) :: r     
    INTEGER, INTENT(OUT) :: nb 
    REAL, INTENT(OUT) :: F,dt,fpeak    

    OPEN(20,FILE=infile,STATUS='UNKNOWN',ACTION='READ')
    READ(20,'(t57,4i8)') Nx,Ny,Nz,nlayers
    READ(20,'(t57,f8.1)') dx
    READ(20,'(t57,f8.2)') tmax
    READ(20,'(t57,i8)') dt_effect
    READ(20,'(t57,f8.1)') G
    READ(20,'(t57,f8.1)') r
    READ(20,'(t57,i8)') nb
    READ(20,'(t57,f8.4)') F
    READ(20,'(t57,f8.5)') dt
    READ(20,'(t57,f8.1)') fpeak
    CLOSE(20)

  END SUBROUTINE input1_3D



  SUBROUTINE input2_3D(geometryfile,nshots,nsmin_x,nsmin_y,nsmin_z,&
       dns_x,dns_y,dns_z,nphones,npmin_x,npmin_y,npmin_z,&
       dnp_x,dnp_y,dnp_z,shot_to_save,nsnaps,snapmin,dsnap,yslice)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: geometryfile
    INTEGER, INTENT(OUT) :: nshots,nsmin_x,nsmin_y,nsmin_z,&
         dns_x,dns_y,dns_z,nphones,npmin_x,npmin_y,npmin_z,&
         dnp_x,dnp_y,dnp_z,shot_to_save,nsnaps,snapmin,dsnap,yslice

    OPEN(21,FILE=geometryfile,STATUS='UNKNOWN',ACTION='READ')
    READ(21,'(t57,i5)') nshots
    READ(21,'(t57,i5)') nsmin_x
    READ(21,'(t57,i5)') nsmin_y
    READ(21,'(t57,i5)') nsmin_z
    READ(21,'(t57,i5)') dns_x
    READ(21,'(t57,i5)') dns_y
    READ(21,'(t57,i5)') dns_z
    READ(21,'(t57,i5)') nphones
    READ(21,'(t57,i5)') npmin_x
    READ(21,'(t57,i5)') npmin_y
    READ(21,'(t57,i5)') npmin_z
    READ(21,'(t57,i5)') dnp_x
    READ(21,'(t57,i5)') dnp_y
    READ(21,'(t57,i5)') dnp_z
    READ(21,'(t57,i5)') shot_to_save
    READ(21,'(t57,i5)') nsnaps
    READ(21,'(t57,i5)') snapmin
    READ(21,'(t57,i5)') dsnap
    READ(21,'(t57,i5)') yslice
    CLOSE(21)

  END SUBROUTINE input2_3D


  SUBROUTINE input3_3D(parfile,nlayers,cp,cs,rhos,phi,kap,cpf,rhof,eta)
    IMPLICIT NONE
    CHARACTER(LEN=20), INTENT(IN) :: parfile
    INTEGER, INTENT(IN) :: nlayers
    REAL, INTENT(OUT), DIMENSION(:) :: cp,cs,rhos,phi,kap,cpf,rhof,eta
    INTEGER :: i

    OPEN(91,FILE=parfile,STATUS='UNKNOWN',ACTION='READ')
    DO i=1,nlayers
       READ(91,*) cp(i)
       READ(91,*) cs(i)
       READ(91,*) rhos(i)
       READ(91,*) phi(i)
       READ(91,*) kap(i)
       READ(91,*) cpf(i)
       READ(91,*) rhof(i)
       READ(91,*) eta(i)
    END DO

  END SUBROUTINE input3_3D



  SUBROUTINE startspanel3D(A,Nx,Ny,Nz,infile)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: A
    CHARACTER(LEN=20), INTENT(IN) :: infile
    INTEGER, INTENT(IN) :: Nx,Ny,Nz
    INTEGER :: i,j

    OPEN(30,FILE=infile,STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO j=1,Ny 
       DO i=1,Nx 
          READ(30, REC=i+(j-1)*Nx) A(i,j,:)
       END DO
    END DO
    CLOSE(30)

  END SUBROUTINE startspanel3D


END MODULE mdle_io_utils3D



