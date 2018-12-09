module Boundary_conditions

  ! Verifica delle condizioni al contorno da applicare

  use problema_vero

  use conversione

  implicit none

contains

  subroutine BC_dx(stato_dx, BC)

    real(kind = 8), dimension(:), intent(in) :: stato_dx

    real(kind = 8), dimension(size(stato_dx), size(stato_dx)), intent(inout) :: BC

    real(kind = 8), dimension(size(stato_dx)) :: lambda

    real(kind = 8), dimension(size(stato_dx), size(stato_dx)) :: R, L

    real(kind = 8), dimension(size(stato_dx)) :: stato_ext,& ! stato_ext = [rho_ext, m_ext]
                                                 delta_w,& ! incremento delle variabili conservative
                                                 delta_v,& ! incremento variabili caratteristiche
                                                 BC_cons ! BC in forma conservativa


    call eigenstructure(stato_dx, lambda, R, L)

    !write(*,*) " BC_fis(:, 2) prima di trasformazioni stato dx = ", BC(:, 2)

    call fis2cons(BC(:, 2), BC_cons) ! Calcolo la m all' uscita dell'estremo destro
                                                  ! BC = [ altezza , velocità ]

    delta_w = BC_cons - stato_dx

    delta_v = matmul(L, delta_w) ! passo alle variabili caratteristiche

    delta_v = delta_v * 0.5d0 * (1 - sign(1.d0, lambda)) ! controllo sul segno delle velocità di propagazione
                                                      ! applico condizioni al contorno solo se
                                                      ! è minore o uguale a 0

    delta_w = matmul(R, delta_v) ! torno alle variabili conservative

    BC_cons = stato_dx + delta_w

    call cons2fis(BC_cons, BC(:, 2)) ! Torno alle variabili fisiche
                                                    ! per le condizioni al contorno
                                                    ! BC(1, 2) = h

    !write(*,*) " BC_fis(:, 2) dopo trasformazioni stato dx = ", BC(:, 2)


  endsubroutine BC_dx

  subroutine BC_sx(stato_sx, BC)

    real(kind = 8), dimension(:), intent(in) :: stato_sx

    real(kind = 8), dimension(size(stato_sx), size(stato_sx)), intent(inout) :: BC

    real(kind = 8), dimension(size(stato_sx)) :: lambda

    real(kind = 8), dimension(size(stato_sx), size(stato_sx)) :: R, L

    real(kind = 8), dimension(size(stato_sx)) :: stato_ext,& ! stato_ext = [rho_ext, m_ext]
                                                 delta_w,& ! incremento delle variabili conservative
                                                 delta_v,& ! incremento variabili caratteristiche
                                                 BC_cons ! BC in forma conservativa


    call eigenstructure(stato_sx, lambda, R, L)

    write(*,*) " BC_fis(:, 1) prima di trasformazioni stato sx = ", BC(:, 1)


    call fis2cons(BC(:, 1), BC_cons) ! Calcolo la m all' uscita dell'estremo sinistro
                                                  ! BC = [ velocità, pressione ]

    delta_w = BC_cons - stato_sx

    !write(*,*) "R = ", R

    !write(*,*) "prima del checking lambda delta_w = ", delta_w

    delta_v = matmul(L, delta_w) ! passo alle variabili caratteristiche

    !write(*,*) "L = ", L

    !write(*,*) "stato_sx = ", stato_sx

    !write(*,*) "lambda_sx = ", lambda

    where (lambda <= 0) ! controllo sul segno delle velocità di propagazione

        delta_v = 0  ! applico condizioni al contorno solo se è minore o uguale a 0

    endwhere

    delta_w = matmul(R, delta_v) ! torno alle variabili conservative

    BC_cons = stato_sx + delta_w

    !write(*,*) "delta_w = ", delta_w

    call cons2fis(BC_cons, BC(:, 1)) ! Torno alle variabili fisiche
                                                    ! per le condizioni al contorno
                                                    ! BC(1, 1) = u

    !write(*,*) " BC_fis(:, 1) dopo trasformazioni stato sx = ", BC(:, 1)


  endsubroutine BC_sx

endmodule Boundary_conditions
