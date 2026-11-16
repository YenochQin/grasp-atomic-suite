!
!***********************************************************************
!                                                                      *
      module def_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only:  max_radial_points
      implicit none (type, external)
      private :: real64
      real(real64) :: tenmax, expmax, expmin, precis
      real(real64) :: hartree_to_wavenumber_cm, hartree_to_ev, speed_of_light_cm_per_s, einstein_a_au_to_si, einstein_b_au_to_si
      real(real64) :: femtometer_to_bohr, electron_mass_u, b1
      integer :: ioncty, electron_count
      real(real64) :: emn, nuclear_charge
      integer :: ionctyff, nelecff
      real(real64) :: emnff, zff
      integer :: ionctyii, nelecii
      real(real64) :: emnii, zii
      integer :: nelecr
      real(real64) :: speed_of_light_au
      real(real64) :: proton_mass_u, bohr_radius_cm
      integer :: nscf, nsic, nsolv
      real(real64) :: accy
      real(real64), pointer :: wt(:), weight(:)
      integer :: ncmin, ncmax
      integer, pointer :: iccmin(:)
      real(real64) :: vacuum_speed_of_light_au, pi
      real(real64) :: dp(max_radial_points), dq(max_radial_points)
      end module def_c
