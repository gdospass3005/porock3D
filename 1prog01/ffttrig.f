
       subroutine ffttrig(nfft,trig)
       complex trig(1),cc
         dinc=2.*3.14159265/nfft
         do n=1,nfft
          cc=dinc*(n-1)
          trig(n)=cexp((0.,1.)*cc)
         end do
         return
       end
