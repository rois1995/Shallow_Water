module fluxes

  ! Modulo per il calcolo dei flussi

  use problema_vero

  use Boundary_conditions

  use conversione

  implicit none

contains

  subroutine flux_LW(w, Dt, Dx_Vec, BC, LW)

    real(kind = 8), dimension(:, :), intent(inout) :: w, BC

    real(kind = 8),                  intent(in) :: Dt

    real(kind = 8), dimension(:),    intent(in) :: Dx_vec

    real(kind = 8), dimension(size(w, 1), size(w, 2) + 1), intent(out) :: LW

    real(kind = 8), dimension(size(w, 1), size(w, 1), size(w, 2) - 1) :: A

    real(kind = 8), dimension(size(w, 1), size(w, 2)) :: flusso

    real(kind = 8), dimension(size(w, 1), size(w, 2) - 1) :: w_medio, flusso_medio,&
                                                             flusso_meno_medio!, A_appoggio

    !real(kind = 8), dimension(size(w, 1)) :: A_appoggio
    real(kind = 8), dimension(size(w, 1), size(w, 1)) :: A_appoggio

    real(kind = 8), dimension(size(w, 1), 2) :: BC_fis

    real(kind = 8), dimension(size(Dx_vec)) :: C
    integer :: j

    C = Dt / (2 * Dx_vec)

    call flusso_vero(w, flusso)

    ! BC brutte

    w_medio = (w(:, 1:size(w, 2) - 1) + w(:, 2:size(w, 2))) / 2

    call jacobian(w_medio, A)

    flusso_medio = (flusso(:, 1:size(w, 2) - 1) + flusso(:, 2:size(w, 2))) / 2
    flusso_meno_medio = (flusso(:, 2:size(w, 2)) - flusso(:, 1:size(w, 2) - 1)) / 2

    do j = 1, size(A, 3)

       A_appoggio = A(:, :, j)

       LW(:, j + 1) = flusso_medio(:, j) -&
                      C * matmul(A_appoggio, flusso_meno_medio(:, j))

    enddo

    !call cons2fis(BC, BC_fis)

    !write(*,*) "BC_fis = ", BC_fis

    !call BC_dx(w(:, size(w, 2)), BC_fis(:, 2))  !! ERRORE QUI, NON SO PERCHè
    !call BC_sx(w(:, 1), BC_fis(:, 1))

    !write(*,*) "BC_fis(:, 1) = ", BC_fis(:, 1)

    !call fis2cons(BC_fis, BC)

    !write(*,*) "BC = ", BC

    call cons2fis(BC, BC)

    !write(*,*) "BC_fis = ", BC

    call BC_dx(w(:, size(w, 2)), BC(:, 2))
    call BC_sx(w(:, 1), BC(:, 1))

    !write(*,*) "BC_fis(:, 1) = ", BC(:, 1)

    call fis2cons(BC, BC)

    !write(*,*) "BC = ", BC

    call flusso_vero(BC(:, 1), LW(:, 1))
    !print *, "BC(:, 1) = ", BC(:, 1)
    !print *, "BC(:, 2) = ", BC(:, 2)
    call flusso_vero(BC(:, 2), LW(:, size(LW, 2)))


  endsubroutine flux_LW



endmodule fluxes
