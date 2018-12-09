module slope_calc

  ! Calcolo delle pendenze di Upwind e di Lax-Wendroff

  use Boundary_conditions

  implicit none

contains

  subroutine Up_slope(Np, order, slope_vec)

    ! Calcolo del vettore delle pendenze col metodo Upwind

    integer, intent(in) :: Np, order

    real(kind = 8), dimension(order, order, Np), intent(inout) :: slope_vec

    slope_vec = 0


  endsubroutine Up_slope




  subroutine LW_slope(Dx_vec, stato, BC_fis, BC_cons, slope_vec_for, slope_vec_back, V, R)

    ! Calcolo del vettore delle pendenze col metodo Lax-Wendroff

    real(kind = 8), dimension(:), intent(in) :: Dx_vec

    real(kind = 8), dimension(:, :), intent(inout) :: BC_cons, BC_fis, stato

    real(kind = 8), dimension(size(BC_cons, 1), size(BC_cons, 2),&
    size(stato, 2)), intent(inout) :: slope_vec_for, slope_vec_back

    real(kind = 8), dimension(size(BC_cons, 1),&
    size(stato, 2) - 1) :: Beta_vec_for, Beta_vec_back

    real(kind = 8), dimension(size(BC_cons, 1), size(stato, 2)), intent(inout) :: V

    real(kind = 8), dimension(size(BC_cons, 1)) :: stato_dx, stato_sx

    real(kind = 8), dimension(size(BC_cons, 1), size(BC_cons, 1), size(stato, 2)), intent(inout) :: R

    real(kind = 8), dimension(size(BC_cons, 1), size(BC_cons, 1), size(stato, 2)) :: L

    real(kind = 8), dimension(size(BC_cons, 1), size(BC_cons, 1)) :: L_app, B

    real(kind = 8), dimension(size(BC_cons, 1)) :: R_app

    real(kind = 8), dimension(size(BC_cons, 1), size(stato, 2) - 1) :: delta_V_for, delta_V_back

    real(kind = 8), dimension(size(BC_cons, 1), size(stato, 2) - 1) :: delta_U


    integer :: dim, j, k


    delta_U = stato(:, 2:size(stato, 2)) - stato(:, 1:size(stato, 2) - 1)

    ! Mi porto nelle variabili caratteristiche

    do j=  1, size(stato, 2)

      call R_calc(stato(:, j), R(:, :, j))

      B = R(:, :, j)

      call L_calc(B, L(:, :, j))

      L_app = L(:, :, j)

      V(:, j) = matmul(L_app, stato(:, j))

    enddo

    do j=  1, size(stato, 2) - 1

      L_app = L(:, :, j + 1)

      delta_V_for(:, j) = matmul(L_app, delta_U(:, j))

      L_app = L(:, :, j)

      delta_V_back(:, j) = matmul(L_app, delta_U(:, j))

    enddo

    dim = size(stato, 2)

    do k = 1, size(stato, 1)

       Beta_vec_for(k, :) = delta_V_for(k, :)/Dx_vec(1:dim - 1)

       Beta_vec_back(k, :) = delta_V_back(k, :)/Dx_vec(2:dim)

       ! Torno alle pendenze nelle variabili conservative

       do j = 1, dim - 1

          R_app = R(:, k, j + 1)
          slope_vec_for(:, k, j) = R_app * Beta_vec_for(k, j)

          R_app = R(:, k, j)
          slope_vec_back(:, k, j + 1) = R_app * Beta_vec_back(k, j)

       enddo

    enddo

    slope_vec_for(:, :, size(slope_vec_for, 3)) = slope_vec_for(:, :, size(slope_vec_for, 3) - 1)
    slope_vec_back(:, :, 1) = slope_vec_back(:, :, 2)

    stato_dx = stato(:, size(stato, 2))

    call BC_dx(stato_dx, BC_fis)

    stato_sx = stato(:, 1)

    call BC_sx(stato_sx, BC_fis)

    call fis2cons(BC_fis, BC_cons)


    !write(*,*) "BC_dx = ", BC_cons(:, 2)

  endsubroutine LW_slope

endmodule slope_calc
