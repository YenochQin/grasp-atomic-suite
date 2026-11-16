!***********************************************************************
!                                                                      *
      subroutine ypot(j)
!                                                                      *
!   This subroutine  tabulates the  potential function Y(r) (Eq (14)   *
!   in  I P Grant, B J McKenzie, P H Norrington, D F Mayers, and N C   *
!   Pyper, Computer  Phys  Commun  21 (1980) 211) for orbital J. The   *
!   function is tabulated in the COMMON array YP.                      *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW, YZK                                     *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: orbital_index_radix
      use debug_c
      use grid_c
      use npot_c, only: double, first_grid_index_beyond_fermi_c, max_radial_points, negative_radius_times_potential
      use orb_c
      use pote_c, only: double, max_radial_points, xp, xq, yp
      use scf_c
      use tatb_c, only: double, mtp, radial_workspace_capacity, ta, tb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yzk_i, only: yzk
      use draw_i, only: draw
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, index, label, k, ioy1, ioy2, nb3, nrows, ii, ii1, ii2, ii3
      real(real64) :: coeff
!-----------------------------------------------
!   Debug printout: composition

      if (ldbpr(29) .or. ldbpr(30)) write (99, 300) orbital_principal_n(j), orbital_symmetry_labels(j)
!
!   Initialize array YP with the nuclear potential piece
!
!   Since YA() below  contains contributions from THIS node only,
!   the initialization should be in consistence with that.

      yp(:radial_point_count) = negative_radius_times_potential(:radial_point_count)

      do index = 1, nycof

         ! Decode information in label
         label = nya(index)
         k = mod(label,key)
         label = label/key
         ioy1 = mod(label,key)
         ioy2 = label/key
         coeff = ya(index)

         if (ldbpr(29)) write (99, 301) k, coeff, orbital_principal_n(ioy1), orbital_symmetry_labels(ioy1), &
            & orbital_principal_n(ioy2)&
            , orbital_symmetry_labels(ioy2)

         call yzk (k, ioy1, ioy2)                ! Accumulate contributions
         yp(:radial_point_count) = yp(:radial_point_count) - coeff*tb(:radial_point_count)
      end do
!
!   Debug printout
!
      if (ldbpr(30)) then
         write (99, 302)
         nb3 = radial_point_count/3
         if (3*nb3 == radial_point_count) then
            nrows = nb3
         else
            nrows = nb3 + 1
         end if
         do ii = 1, nrows
            ii1 = ii
            ii2 = ii1 + nrows
            ii3 = ii2 + nrows
            if (ii3 <= radial_point_count) then
               write (99, 303) radial_grid(ii1), yp(ii1), radial_grid(ii2), yp(ii2), radial_grid(ii3), yp(ii3&
                  )
            else if (ii2 <= radial_point_count) then
               write (99, 303) radial_grid(ii1), yp(ii1), radial_grid(ii2), yp(ii2)
            else
               write (99, 303) radial_grid(ii1), yp(ii1)
            end if
         end do
         call draw (yp, 1.0_real64, yp, 0.0_real64, radial_point_count)
      end if
!
      return
!
  300 format(/,/," Direct potential for ",1i2,1a2," orbital :"/,/)
  301 format(/,25x,"(",1i2,")"/,1x,1p,d21.14,"* Y    (",1i2,1a2,",",1i2,1a2,")"&
         )
  302 format(/,/,3(" --------- r --------- ------- Y (r) -------"))
  303 format(1p,6(1x,1d21.14))
      return
!
      end subroutine ypot
