module nucpot_i
    implicit none (type, external)
    private
    public :: nucpot
contains
!***********************************************************************
!                                                                      *
      subroutine nucpot
!                                                                      *
!   Evaluate the nuclear potential for point and Fermi models.         *
!                                                                      *
!   Call(s) to: [LIB92] ES.                                            *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 05 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use def_c,           only: nuclear_charge, pi
      use grid_c, only: double, radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use npot_c,          only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use es_i, only: es
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nb3, nrows, ii, ii1, ii2, ii3
      real(real64) :: c, a, abc, tabc, abc2, thabc2, abc3, cba, pi2, hpiac2, &
         h3php, s2mcba, s3mcba, sabc3, dmsas, en, zbn, ri, rmc, rmcba, rbc, &
         s2rcba, s3rcba
      logical :: set
!-----------------------------------------------
!
!
!   Point nucleus
!
      if (nuclear_model_parameter_count == 0) then
!
         negative_radius_times_potential(:radial_point_count) = nuclear_charge
!
!   Fermi distribution
!
      else if (nuclear_model_parameter_count == 2) then
!
         c = nuclear_model_parameters(1)
         a = nuclear_model_parameters(2)
         abc = a/c
         tabc = 2.0_real64*abc
         abc2 = abc*abc
         thabc2 = 3.0_real64*abc2
         abc3 = abc2*abc
         cba = c/a
         pi2 = pi*pi
         hpiac2 = 0.5_real64*pi2*abc2
         h3php = 1.5_real64 + hpiac2
         call es ((-cba), s2mcba, s3mcba)
         sabc3 = 6.0_real64*abc3
         dmsas = -sabc3*s3mcba
         en = 1.0_real64 + abc2*pi2 + dmsas
         zbn = nuclear_charge/en
!
         set = .false.
         do i = 1, radial_point_count
            ri = radial_grid(i)
            rmc = ri - c
            rmcba = rmc/a
            rbc = ri/c
            if (rbc <= 1.0_real64) then
               call es (rmcba, s2rcba, s3rcba)
               negative_radius_times_potential(i) = zbn*(dmsas + sabc3*s3rcba + rbc*(h3php - thabc2*s2rcba&
                   - 0.5_real64*rbc*rbc))
            else
               if (.not.set) then
                  first_grid_index_beyond_fermi_c = i
                  set = .true.
               end if
               call es ((-rmcba), s2rcba, s3rcba)
               negative_radius_times_potential(i) = nuclear_charge*(1.0_real64 + thabc2*(rbc*s2rcba + tabc*s3rcba)/en)
            end if
         end do
      end if
!
      if (ldbpr(2)) then
         write (99, 300)
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
               write (99, 301) radial_grid(ii1), negative_radius_times_potential(ii1), radial_grid(ii2), &
                  & negative_radius_times_potential(ii2), radial_grid(ii3), negative_radius_times_potential(ii3&
                  )
            else if (ii2 <= radial_point_count) then
               write (99, 301) radial_grid(ii1), negative_radius_times_potential(ii1), radial_grid(ii2), &
                  & negative_radius_times_potential(ii2)
            else
               write (99, 301) radial_grid(ii1), negative_radius_times_potential(ii1)
            end if
         end do
      end if
!
      return
!
  300 format(/,"From SUBROUTINE NUCPOT:"/,3(&
         " -------- r -------- ----- -r*V(r) -----"))
  301 format(1p,6(1x,1d19.12))
      return
!
      end subroutine nucpot
end module nucpot_i
