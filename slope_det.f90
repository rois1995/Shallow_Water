module slope_det

  ! Modulo per la determinazione della pendenza

  use slope_calc

  use checking1

  implicit none

contains

  subroutine slope_fin(Dx_vec, stato, BC_fis, sslope_def, BC_cons, V, R)

    real(kind = 8), dimension(:), intent(in) :: Dx_vec

    real(kind = 8), dimension(:, :), intent(inout) :: stato, BC_fis

    real(kind = 8), dimension(size(stato, 1), size(stato, 1), size(stato, 2)), intent(inout) :: sslope_def, R

    real(kind = 8), dimension(size(stato, 1), size(stato, 1), size(stato, 2)) :: sslope_LW_for,&
                                                                                 sslope_LW_back,&
                                                                                 sslope_Up

    real(kind = 8), dimension(size(BC_fis, 1), size(BC_fis, 2)), intent(inout) :: BC_cons

    real(kind = 8), dimension(size(stato, 1), size(stato, 2)), intent(inout) :: V

    integer :: j

     !-----------------------------------------------------------------------------
     ! Calcolo vettore delle pendenze per rho e per m
     !-----------------------------------------------------------------------------

     call Up_slope(size(stato, 2), size(stato, 1), sslope_Up) ! Vettore delle pendenze di Upwind (tutte 0)

     ! Vettore delle pendenze di Lax_Wendroff

     call LW_slope(Dx_vec, stato, BC_fis, BC_cons, sslope_LW_for, sslope_LW_back, V, R)

     do j = 1, size(Dx_vec)

       call slope(sslope_LW_for(:, :, j), sslope_LW_back(:, :, j), sslope_def(:, :, j))

     enddo

endsubroutine slope_fin

endmodule slope_det
