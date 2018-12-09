module RK

  ! Integratore nel tempo tramite Runge-Kutta a 2 ed a 3 passi

  use fluxes

  use slope_det

  implicit none


contains

  function cose(stato, k, flag, lambda, sslope, BC_fis, V, R, Dx_vec) result(U_next)

    real(kind = 8), dimension(:, :, :), intent(inout) :: R, sslope

    real(kind = 8), dimension(:, :), intent(inout) :: stato, V, lambda

    real(kind = 8),                  intent(in) :: k

    real(kind = 8), dimension(:),    intent(in)   :: Dx_vec

    integer,                         intent(in) :: flag

    real(kind = 8), dimension(:, :), intent(inout) :: BC_fis

    real(kind = 8), dimension(size(BC_fis,1), size(BC_fis,2)) :: BC_cons

    real(kind = 8), dimension(size(stato, 1), size(stato, 2)) :: U_next

    real(kind = 8), dimension(size(stato, 1), size(stato, 2)) :: U_app_1, U_app_2,&
                                                                 U_app_3, U_app_4,&
                                                                 U_app_1_2, U_app_2_2,&
                                                                 U_app_3_2, U_app_4_2,&
                                                                 d_flux, Dx_mat

    real(kind = 8), dimension(size(stato, 1), size(stato, 2) + 1) :: flusso

    real(kind = 8) :: a = 1.d0

    integer :: j

    ! Calcolo delle pendenze

    call slope_fin(Dx_vec, stato, BC_fis, sslope, BC_cons, V, R)

    ! Calcolo del vettore differenza dei flussi

    call fluxo(stato, k, Dx_vec, sslope, BC_cons, V, R, flusso)

    d_flux = flusso(:, 2:size(flusso, 2)) - flusso(:, 1:size(flusso, 2) - 1)

    do j = 1, size(Dx_mat, 1)

       Dx_mat(j, :) = Dx_vec

    enddo


    select case (flag)

    case (1)

      write(*,*) "-----------------------------------------------"
      write(*,*) "No Runge-Kutta"
      write(*,*) "-----------------------------------------------"

      U_next = stato - k * d_flux / Dx_mat


    case (2)

      write(*,*) "-----------------------------------------------"
      write(*,*) "Esecuzione del metodo di Runge-Kutta a 2 passi"
      write(*,*) "-----------------------------------------------"

        ! Primo passo

        U_app_1 = stato - 0.5d0 * k * d_flux / Dx_mat

        ! Calcolo delle pendenze

        call slope_fin(Dx_vec, U_app_1, BC_fis, sslope, BC_cons, V, R)

        ! Ricalcolo i flussi

        call fluxo(U_app_1, k, Dx_vec, sslope, BC_cons, V, R, flusso)

        d_flux = flusso(:, 2:size(flusso, 2)) - flusso(:, 1:size(flusso, 2) - 1)

        ! Secondo passo

        U_next = stato - 0.5d0 * k * d_flux / Dx_mat


    case (3)

      write(*,*) "-----------------------------------------------"
      write(*,*) "Esecuzione del metodo di Runge-Kutta a 3 passi"
      write(*,*) "-----------------------------------------------"

        ! Primo passo

        U_app_1 = stato - 0.5d0 * k * d_flux / Dx_mat

        ! Calcolo delle pendenze

        call slope_fin(Dx_vec, U_app_1, BC_fis, sslope, BC_cons, V, R)

        ! Ricalcolo i flussi

        call fluxo(U_app_1, k, Dx_vec, sslope, BC_cons, V, R, flusso)

        d_flux = flusso(:, 2:size(flusso, 2)) - flusso(:, 1:size(flusso, 2) - 1)

        ! Secondo passo

        U_app_2 = (3/4.d0) * stato + 0.25 * U_app_1 - (k / 8.d0) * d_flux / Dx_mat


        ! Calcolo delle pendenze

        call slope_fin(Dx_vec, U_app_2, BC_fis, sslope, BC_cons, V, R)

        ! Ricalcolo i flussi

        call fluxo(U_app_2, k, Dx_vec, sslope, BC_cons, V, R, flusso)

        d_flux = flusso(:, 2:size(flusso, 2)) - flusso(:, 1:size(flusso, 2) - 1)

        ! Secondo passo

        U_next = (1/3.d0) * stato + (2/3.d0) * U_app_2 - (1/3.d0) * k * d_flux / Dx_mat

        !!!!!!!!!!!!!!!!!!!!!!!!!! Controlloare se è giusto mettere i /Dx_mat anche in RK3 e RK4


     case (4)

      write(*,*) "-----------------------------------------------"
      write(*,*) "Esecuzione del metodo di Runge-Kutta a 4 passi"
      write(*,*) "-----------------------------------------------"

        ! Primo passo

        U_app_1 = - 0.5d0 * k * d_flux / Dx_mat

        U_app_1_2 = U_app_1 / 2.d0 + stato

        ! Calcolo delle pendenze

        call slope_fin(Dx_vec, U_app_1_2, BC_fis, sslope, BC_cons, V, R)

        ! Ricalcolo i flussi

        call fluxo(U_app_1_2, k, Dx_vec, sslope, BC_cons, V, R, flusso)

        d_flux = flusso(:, 2:size(flusso, 2)) - flusso(:, 1:size(flusso, 2) - 1)

        ! Secondo passo

        U_app_2 = - 0.5d0 * k * d_flux / Dx_mat

        U_app_2_2 = U_app_2 / 2.d0 + stato

        ! Calcolo delle pendenze

        call slope_fin(Dx_vec, U_app_2_2, BC_fis, sslope, BC_cons, V, R)

        ! Ricalcolo i flussi

        call fluxo(U_app_2_2, k, Dx_vec, sslope, BC_cons, V, R, flusso)

        d_flux = flusso(:, 2:size(flusso, 2)) - flusso(:, 1:size(flusso, 2) - 1)

        ! Terzo passo

        U_app_3 = - 0.5d0 * k * d_flux / Dx_mat

        U_app_3_2 = U_app_3 + stato

        ! Calcolo delle pendenze

        call slope_fin(Dx_vec, U_app_3_2, BC_fis, sslope, BC_cons, V, R)

        ! Ricalcolo i flussi

        call fluxo(U_app_3_2, k, Dx_vec, sslope, BC_cons, V, R, flusso)

        d_flux = flusso(:, 2:size(flusso, 2)) - flusso(:, 1:size(flusso, 2) - 1)

        ! Quarto passo

        U_app_4 = - 0.5d0 * k * d_flux / Dx_mat

        ! Update della soluzione

        U_next = stato + (1/6.d0) * (U_app_1 + 2 * U_app_2 + 2 * U_app_3 + U_app_4)

    endselect



  end function

endmodule RK
