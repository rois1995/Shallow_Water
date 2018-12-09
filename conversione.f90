module conversione

  ! Modulo per la conversione da variabili conservative
  ! a variabili fisiche e viceversa
  ! w = variabili conservative = [a, au]
  ! p = variabili fisiche = [a, u]

  implicit none

  interface cons2fis

    module procedure cons2fis_s, cons2fis_v

  endinterface cons2fis

  interface fis2cons

    module procedure fis2cons_s, fis2cons_v

  endinterface fis2cons

contains

  subroutine cons2fis_s(w, p)

    real(kind = 8), dimension(:), intent(in) :: w

    real(kind = 8), dimension(size(w)), intent(out) :: p

    p(1) = w(1)

    p(2) = w(2) / w(1)

  endsubroutine cons2fis_s

  subroutine cons2fis_v(ww, pp)

    real(kind = 8), dimension(:, :), intent(in) :: ww

    real(kind = 8), dimension(size(ww, 1), size(ww, 2)), intent(out) :: pp

    pp(1, :) = ww(1, :)

    pp(2, :) = ww(2, :) / ww(1, :)

  endsubroutine cons2fis_v

  subroutine fis2cons_s(p, w)

    real(kind = 8), dimension(:), intent(in) :: p

    real(kind = 8), dimension(size(p)), intent(out) :: w

    w(1) = p(1)

    w(2) = p(2) * p(1)

  endsubroutine fis2cons_s

  subroutine fis2cons_v(pp, ww)

    real(kind = 8), dimension(:, :), intent(in) :: pp

    real(kind = 8), dimension(size(pp, 1), size(pp, 2)), intent(out) :: ww

    ww(1, :) = pp(1, :)

    ww(2, :) = pp(2, :) * pp(1, :)

  endsubroutine fis2cons_v

endmodule conversione
