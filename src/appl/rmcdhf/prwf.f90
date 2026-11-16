!***********************************************************************
!                                                                      *
      subroutine prwf(j)
!                                                                      *
!   Makes a (debug) printout of wave functions. There are two modes:   *
!                                                                      *
!      (1) J>0  - Used as a debug option in SOLVE, wavefunctions for   *
!                 orbital J are printed                                *
!      (2) J=0  - A  printout  of the grid and all wave functions is   *
!                 made                                                 *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_radial_points
      use def_c
      use grid_c
      use int_c
      use orb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE draw_I
      use draw_i, only: draw
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nb2, nrows, ii, ii1, ii2, k, mfk, jgg
      real(real64) :: cbz
!-----------------------------------------------
!
!
      cbz = speed_of_light_au/nuclear_charge
!
      if (j > 0) then
!
!   Mode (1)
!
         write (99, 300) orbital_principal_n(j), orbital_symmetry_labels(j)
         nb2 = mtp0/2
         if (2*nb2 == mtp0) then
            nrows = nb2
         else
            nrows = nb2 + 1
         end if
         do ii = 1, nrows
            ii1 = ii
            ii2 = ii1 + nrows
            if (ii2 <= mtp0) then
               write (99, 301) radial_grid(ii1), p(ii1), q(ii1), radial_grid(ii2), p(ii2), q(ii2)
            else if (ii1 <= mtp0) then
               write (99, 301) radial_grid(ii1), p(ii1), q(ii1)
            end if
         end do
         call draw (p, 1.0_real64, q, cbz, mtp0)
!
      else
!
!   Mode (2)
!
         do k = 1, orbital_count
            write (99, 300) orbital_principal_n(k), orbital_symmetry_labels(k)
            mfk = orbital_last_grid_index(k)
            nb2 = mfk/2
            if (2*nb2 == mfk) then
               nrows = nb2
            else
               nrows = nb2 + 1
            end if
            do ii = 1, nrows
               ii1 = ii
               ii2 = ii1 + nrows
               if (ii2 <= mfk) then
                  write (99, 301) radial_grid(ii1), radial_large_components(ii1,k), radial_small_components(ii1,k), &
                     & radial_grid(ii2), radial_large_components(ii2,&
                     k), radial_small_components(ii2,k)
               else if (ii1 <= mfk) then
                  write (99, 301) radial_grid(ii1), radial_large_components(ii1,k), radial_small_components(ii1,k)
               end if
            end do
            call draw (radial_large_components(:max_radial_points,k), 1.0_real64, &
               & radial_small_components(:max_radial_points,k), cbz, orbital_last_grid_index(k))
!
         end do
!
      end if
!
      return
!
  300 format("1",1i2,1a2," orbital:"/,2(&
         " --------- r --------- ------- P (r) -------",&
         " ------- Q (r) -------"))
  301 format(1p,6(1x,1d21.14))
      return
!
      end subroutine prwf
