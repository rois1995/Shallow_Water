PROGRAM SW_completo

  ! Solutore per le Shallow-Water Equations in 1D

  USE gnufor_mio

  USE LOGO

  use problema_vero

  use fluxes

  ! Discretizzazione nello spazio

  INTEGER, parameter :: Np = 100

  real(kind = 8), dimension(Np) :: x_vec

  real(kind = 8), dimension(Np - 2) :: Dx_vec

  real(kind = 8) :: Dx, t, Dt

  real(kind = 8), parameter :: t0 = 0, Tf = 10.d0, xs = 0.d0, xd = 10.d0

  real(kind = 8) :: xm = (xd + xs)/2.d0

  real(kind = 8), parameter :: CFL = 0.9

  ! Problema fisico

  real(kind = 8), dimension(Np) :: h_i, zb = 0.0d0, u_i = 0.0d0, d_zb = 0.0d0, cf = 0.1

  real(kind = 8) :: cf_wind = 0.001, W_wind = 10.d0

  real(kind = 8), dimension(2, Np) :: w, w_old, lambda, s

  real(kind = 8), dimension(2, 2, Np) :: R, L

  real(kind = 8), dimension(2, 2) :: BC

  real(kind = 8), dimension(2, Np + 1) :: LW, UP

  real(kind = 8), dimension(Np, 1000) :: h_mat, u_mat, x_mat

  ! Variabili utili per i plot

  CHARACTER(len = 50) :: title

  CHARACTER(len = 50), dimension(:), allocatable :: legend

  ! Interi utili

  INTEGER :: j

  !----------------------------------------------------------------------------
  ! LOGO FIGO
  !----------------------------------------------------------------------------

  call disegna_logo

  !----------------------------------------------------------------------------
  ! DEFINIZIONE DEL DOMINIO DI CALCOLO
  !----------------------------------------------------------------------------

  ! Costruzione del dominio di calcolo

  Dx = abs(xd - xs) / (Np - 1)

  x_vec(1) = xs;   x_vec(Np) = xd

  do j = 2, Np - 1

    x_vec(j) = x_vec(j - 1) + Dx

    if ( j == 2 .or. j == Np - 1 )  then

       Dx_vec(j - 1) = Dx/2

    else

       Dx_vec(j - 1) = Dx

    endif

  enddo

  !h_i = 1 + 0.02 * sin(2 * 3.14d0 * (xd - x_vec)/ (xd / 4)) ! Condizione iniziale sull'altezza sinusoidale
  h_i = 1 - 0.02 * exp(-(x_vec - xm) ** 2) ! Condizione iniziale sull'altezza gaussiana

  if(allocated(legend)) deallocate(legend);

  allocate(legend(1))

  legend(1) = 'Altezza pelo libero'

  title = 'boh'

  !call plot_function(x_vec, h_i, legend, title)

  w_old(1, :) = (h_i - zb);   w_old(2, :) = (h_i - zb) * u_i

  h_mat(:, 1) = h_i
  u_mat(:, 1) = u_i
  x_mat(:, 1) = x_vec

  !write(*,*) "w_old(1, :) = ", w_old(1, :)

  j = 2

  do while (t < 7.6)

    call eigenstructure(w_old, lambda, R, L)

    Dt = Dx * CFL / (max(maxval(lambda(1, :)), maxval(lambda(2, :))))

    t = t + Dt

    write(*,*) "---------------------------------------------------------------"
    write(*,*) "Istante t = ", t
    write(*,*) "---------------------------------------------------------------"

    ! Bordi entrambe aperti

    BC(1, 1) = w_old(1, 1);  BC(2, 1) = w_old(2, 1) / w_old(1, 1);  ! Sono scritte in variabili fisiche
    BC(1, 2) = w_old(1, Np);  BC(2, 2) = w_old(2, Np) / w_old(1, Np);

    ! Bordi entrambi chiusi, non funziona per ora

    !BC(1, 1) = w_old(1, 1);  BC(2, 1) = 0;
    !BC(1, 2) = w_old(1, Np);  BC(2, 2) = 0;

    call flux_LW(w_old, Dt, Dx_vec, BC, LW)
    call sorgente(w_old, cf, d_zb, s)

    w = w_old - (Dt / (2 * Dx)) * (LW(:, 2:size(LW, 2)) - LW(:, 1:size(LW, 2) - 1)) -&
        s * Dt


        write(*,*) "LW(:, 1) = ", LW(:, 1)
        write(*,*) "LW(:, 2) = ", LW(:, 2)
        write(*,*) "j = ", j

    h_mat(:, j) = w(1, :) - zb
    u_mat(:, j) = w(2, :) / w(1, :)
    x_mat(:, j) = x_vec

    j = j + 1
    w_old = w


  enddo

  if(allocated(legend)) deallocate(legend);

  allocate(legend(2))

  legend(1) = 'altezza'
  legend(2) = 'fondale'

  !call plot_2functions(x_vec, h_mat(:, 2), zb, legend, title)

  if(allocated(legend)) deallocate(legend);

  allocate(legend(1))

  legend(1) = 'altezza'

  !call plot_video(x_mat, h_mat, title, legend)


endprogram SW_completo
