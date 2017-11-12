MODULE mdle_utils

CONTAINS

  SUBROUTINE quads(nt,aux1)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: aux1
    INTEGER, INTENT(IN)  :: nt
! Quadrando um n. inteiro (nt)
! funciona na faixa 2<nt<1.e09
    INTEGER :: i,pot

    pot=30; aux1=2**pot
    do i=1,30
       if (aux1>=nt) then
          pot = pot-1
          aux1 = 2**pot
       end if
    end do
    aux1=2*aux1

  END SUBROUTINE quads


END MODULE mdle_utils
