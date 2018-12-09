MODULE problema_vero

  ! Modulo contenente i flussi per Shallow-Water Equations in 1D
  ! W = vettore variabili conservative = [a; au]
  ! u = velocità orizzontale
  ! a = h - zb
  ! h = altezza colonna d'acqua dallo 0
  ! zb = altezza fondale

  use inversa

  implicit none

  real(kind = 8), parameter :: g = 9.81d0

  interface flusso_vero

    module procedure flusso_vero_s, flusso_vero_v

  endinterface flusso_vero

  interface sorgente

    module procedure sorgente_s, sorgente_v

  endinterface sorgente

  interface eigenstructure

    module procedure eigenstructure_s, eigenstructure_v

  endinterface eigenstructure

  interface jacobian

    module procedure jacobian_s, jacobian_v

  endinterface jacobian


contains

  subroutine flusso_vero_s(w, flusso)

    real(kind = 8), dimension(:), intent(in) :: w

    real(kind = 8), dimension(size(w)), intent(out) :: flusso

    flusso(1) =  w(2)  ! f(1) = a u

    flusso(2) = (w(2) ** 2) / w(1) + 0.5d0 * g * (w(1) ** 2) ! f(2) = a u^2 + 1/2 g a^2

  endsubroutine flusso_vero_s

  subroutine flusso_vero_v(ww, fflusso)

    real(kind = 8), dimension(:, :), intent(in) :: ww

    real(kind = 8), dimension(size(ww, 1), size(ww, 2)), intent(out) :: fflusso

    fflusso(1, :) = ww(2, :)  ! f(1) = a u

    fflusso(2, :) = (ww(2, :) ** 2) / ww(1, :) + 0.5d0 * g * (ww(1, :) ** 2) ! f(2) = a u^2 + 1/2 g a^2

  endsubroutine flusso_vero_v

  subroutine sorgente_s(w, cf, d_zb, S)

    real(kind = 8), dimension(:), intent(in) :: w

    real(kind = 8),               intent(in) :: cf, d_zb

    real(kind = 8), dimension(size(w)), intent(out) :: S

    S(1) = 0

    S(2) = cf * (w(2) / w(1)) * sqrt(w(1)**2 + w(2)**2) + g * w(1) * d_zb ! S(2) = cf u |w| + g a dx(zb)

  endsubroutine sorgente_s

  subroutine sorgente_v(ww, ccf, dd_zb, SS)

    real(kind = 8), dimension(:, :), intent(in) :: ww

    real(kind = 8), dimension(:),    intent(in) :: ccf, dd_zb

    real(kind = 8), dimension(size(ww, 1), size(ww, 2)), intent(out) :: SS

    SS(1, :) = 0

    SS(2, :) = ccf * (ww(2, :) / ww(1, :)) * sqrt(ww(1, :)**2 + ww(2, :)**2) +&
               g * ww(1, :) * dd_zb ! S(2) = cf u |w| + g a dx(zb)

  endsubroutine sorgente_v

  subroutine jacobian_s(w, J)

    ! Subroutine per il calcolo dello jacobiano in 1D

    real(kind = 8), dimension(:), intent(in) :: w

    real(kind = 8), dimension(size(w), size(w)), intent(out) :: J

    J(1, 1) = 0; J(1, 2) = 1.d0;  ! J(1, :) = [ 0, 1 ]

    J(2, 1) = -(w(2) / w(1))**2 + g * w(1) ! J(2, 1) = - u^2 + ga

    J(2, 2) = 2 * w(2) / w(1) ! J(2, 2) = 2u


  endsubroutine jacobian_s

  subroutine jacobian_v(ww, JJ)

    ! Subroutine per il calcolo dello jacobiano in 1D

    real(kind = 8), dimension(:, :), intent(in) :: ww

    real(kind = 8), dimension(size(ww, 1), size(ww, 1), size(ww, 2)), intent(out) :: JJ

    JJ(1, 1, :) = 0; JJ(1, 2, :) = 1.d0;  ! JJ(1, :, :) = [ 0, 1 ]

    JJ(2, 1, :) = -(ww(2, :) / ww(1, :))**2 + g * ww(1, :) ! JJ(2, 1, :) = - u^2 + ga

    JJ(2, 2, :) = 2 * ww(2, :) / ww(1, :) ! JJ(2, 2, :) = 2u




  endsubroutine jacobian_v

  subroutine eigenstructure_s(w, lambda, R, L)

    real(kind = 8), dimension(:), intent(in) :: w

    real(kind = 8), dimension(size(w)), intent(out) :: lambda

    real(kind = 8), dimension(size(w), size(w)), intent(out) :: R, L

    real(kind = 8), dimension(size(w), size(w)) :: B ! Matrice di appoggio per il calcolo dell'inversa

    lambda(1) = w(2) / w(1) - sqrt(g * w(1))  ! lambda(1) = u - sqrt(ga)
    lambda(2) = w(2) / w(1) + sqrt(g * w(1))  ! lambda(2) = u + sqrt(ga)

    R(1, 1) = 1.d0;   R(1, 2) = 1.d0;   ! R(1, :) = [1, 1]

    R(2, 1) = w(2) / w(1) - sqrt(g * w(1))  ! R(2, 1) = u - sqrt(ga)
    R(2, 2) = w(2) / w(1) + sqrt(g * w(1))  ! R(2, 2) = u + sqrt(ga)

    !R(2, 1) = - sqrt(g * w(1))  ! R(2, 1) = - sqrt(ga)
    !R(2, 2) =   sqrt(g * w(1))  ! R(2, 2) =  sqrt(ga)

    B = R ! Mi appoggio su B

    call inverse(B, L, size(w))

  endsubroutine eigenstructure_s

  subroutine eigenstructure_v(ww, llambda, RR, LL)

    real(kind = 8), dimension(:, :), intent(in) :: ww

    real(kind = 8), dimension(size(ww, 1), size(ww, 2)), intent(out) :: llambda

    real(kind = 8), dimension(size(ww, 1), size(ww, 1), size(ww, 2)), intent(out) :: RR, LL

    real(kind = 8), dimension(size(ww, 1), size(ww, 1)) :: B ! Matrice di appoggio per il calcolo dell'inversa

    integer :: j

    llambda(1, :) = ww(2, :) / ww(1, :) - sqrt(g * ww(1, :))  ! lambda(1) = u - sqrt(ga)
    llambda(2, :) = ww(2, :) / ww(1, :) + sqrt(g * ww(1, :))  ! lambda(2) = u + sqrt(ga)

    RR(1, 1, :) = 1.d0;   RR(1, 2, :) = 1.d0;   ! R(1, :) = [1, 1]

    RR(2, 1, :) = ww(2, :) / ww(1, :) - sqrt(g * ww(1, :))  ! R(2, 1) = u - sqrt(ga)
    RR(2, 2, :) = ww(2, :) / ww(1, :) + sqrt(g * ww(1, :))  ! R(2, 2) = u + sqrt(ga)

    !RR(2, 1, :) = - sqrt(g * ww(1, :))  ! R(2, 1) = - sqrt(ga)
    !RR(2, 2, :) =   sqrt(g * ww(1, :))  ! R(2, 2) =  sqrt(ga)

    do j = 1, size(ww, 2)

       B = RR(:, :, j)

       call inverse(B, LL(:, :, j), size(ww, 1))

    enddo

  endsubroutine eigenstructure_v

ENDMODULE problema_vero
