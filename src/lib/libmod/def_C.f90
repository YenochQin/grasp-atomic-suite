!
!***********************************************************************
!                                                                      *
      module def_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double
      use parameter_def,   only:  max_radial_points
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double) :: tenmax, expmax, expmin, precis
      real(double) :: hartree_to_wavenumber_cm, hartree_to_ev, speed_of_light_cm_per_s, einstein_a_au_to_si, einstein_b_au_to_si
      real(double) :: femtometer_to_bohr, electron_mass_u, b1
      integer :: ioncty, electron_count
      real(double) :: emn, nuclear_charge
      integer :: ionctyff, nelecff
      real(double) :: emnff, zff
      integer :: ionctyii, nelecii
      real(double) :: emnii, zii
      integer :: nelecr
      real(double) :: speed_of_light_au
      real(double) :: proton_mass_u, bohr_radius_cm
      integer :: nscf, nsic, nsolv
      real(double) :: accy
      real(double), dimension(:), pointer :: wt, weight
      integer :: ncmin, ncmax
      integer, dimension(:), pointer :: iccmin
      real(double) :: vacuum_speed_of_light_au, pi
      real(double), dimension(max_radial_points) :: dp, dq
      end module def_c
