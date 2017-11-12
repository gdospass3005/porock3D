	SUBROUTINE FACTORK(IFAX,NFAC,NFFT)
	INTEGER IFAX(1)
	nf=nfft
	nfac=0
C	FACTOR 6
	DO WHILE (NF/6*6.EQ.NF)
	 nfac=nfac+1
	 ifax(nfac)=6
	 nf=nf/6
	END DO
c	FACTOR 4
	DO WHILE (NF/4*4.EQ.NF)
	 nfac=nfac+1
	 ifax(nfac)=4
         nf=nf/4
	END DO
C	FACTOR 5
	DO WHILE (NF/5*5.EQ.NF)
	 nfac=nfac+1
	 ifax(nfac)=5
	 nf=nf/5
	END DO
c	FACTOR 3
	DO WHILE (NF/3*3.EQ.NF)
	 nfac=nfac+1
	 ifax(nfac)=3
         nf=nf/3
	END DO
c	FACTOR 2
	DO WHILE (NF/2*2.EQ.NF)
	 nfac=nfac+1
	 ifax(nfac)=2
         nf=nf/2
	END DO
	if(nf.eq.1) return
	 write(6,*)'***** factorization in factor failed'
	stop 999
	END
