module fix_BC

  ! modulo per fissare le BC (periodiche, stazionarie, muro, parete etc..)
  ! Le condizioni al contorno vengono assegnate in variabili fisiche

  use conversione

  implicit none

contains

  subroutine BCs(flag, stato_sx, stato_dx, BC, Bound_speed, t)

    real(kind = 8), dimension(:), intent(in) :: stato_sx, stato_dx

    integer, dimension(:), intent(in) :: flag

    logical, dimension(:), intent(in) :: Bound_speed

    real(kind = 8), dimension(:), optional :: t ! t = [Tf t_att]

    real(kind = 8), dimension(size(stato_dx), size(stato_dx)), intent(inout) :: BC



    select case (flag(1))

    case (1)

      write(*,*) "A sinistra ho condizioni al contorno Wall"

        ! Si suppone che gli stato abbiano m = 0 e che non ci sia la formazione
        ! di vuoto

        call cons2fis(stato_sx, BC(:, 1))

        if(Bound_speed(1)) then

          if(present(t)) then

            BC(2, 1) = 1 - sin(0.1 * 2 * 3.14d0 * (t(1) - t(2)) / t(1))

          else

            write(*,*) "Manca il tempo come input"

            return

          endif

        else

          BC(2, 1) = 0

        endif

        write(*,*) "BC_s = ", BC(:, 1)

    case(2)

       write(*,*) "A sinistra ho condizioni al contorno Open"

         ! Si suppone che gli stato abbiano m = 0 e che non ci sia la formazione
         ! di vuoto

        call cons2fis(BC(:, 1), BC(:, 1))

    endselect

    !!!!! FIX IL CASO DI WALL PERCHé NON DEVE RICALCOLARLE NELLE ALTRE FUNZIONI.
    !!!! LA U è FISSATA!!

    select case (flag(2))

    case (1)

      write(*,*) "A destra ho condizioni al contorno Wall"

        ! Si suppone che gli stato abbiano m = 0 e che non ci sia la formazione
        ! di vuoto

      call cons2fis(stato_dx, BC(:, 2))

      if(Bound_speed(2)) then

        if(present(t)) then

          BC(2, 2) = 1 - sin(0.1 * 2 * 3.14d0 * (t(1) - t(2)) / t(1))

        else

          write(*,*) "Manca il tempo come input"

          return

        endif

      else

        BC(2, 2) = 0

      endif

      write(*,*) "BC_d = ", BC(:, 2)

    case(2)

      write(*,*) "A destra ho condizioni al contorno Open"

        ! Si suppone che gli stato abbiano m = 0 e che non ci sia la formazione
        ! di vuoto

       call cons2fis(BC(:, 2), BC(:, 2))

    endselect

  endsubroutine



endmodule fix_BC
