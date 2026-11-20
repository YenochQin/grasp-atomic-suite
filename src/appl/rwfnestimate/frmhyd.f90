

!***********************************************************************
!                                                                      *
      subroutine frmhyd(index, nsubs, modify, nradial)
!                                                                      *
!   This  subroutine  is  used  to  produce  estimates  of  the wave   *
!   functions from the hydrogenic approximation. The screening cons-   *
!   tant is taken to be 0.4 (Changed, see comments below,XHH)          *
!                                                                      *
!   Call(s) to: [ERWF]: DCBSRW.                                        *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters, only: radial_workspace_capacity
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: nuclear_charge
      use initial_orbital_selection, only: max_orbitals, set
      use orbital_csf_basis, only: orbital_energy_parameters, orbital_principal_n, orbital_kappa, orbital_symmetry_labels
      use radial_wavefunctions, only: &
          radial_large_components, radial_small_components, orbital_origin_parameters, orbital_last_grid_index
      use initial_orbital_source, only: source
      use hydrogenic_screening, only: sigma
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use calculation_input_defaults, only: inpci, ndef, ndump
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dirac_coulomb_radial_orbital, only: dcbsrw
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nsubs
      logical , intent(in) :: modify
      integer , intent(in) :: index(max_orbitals)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, loc
      integer, intent(in) :: nradial
      real(real64) :: zeff, delz
      logical :: yes, getyn
!-----------------------------------------------
!
!
      delz = 0.0_real64
      if (nradial == 4)  then
         write (stderr_unit, *) "Enter increase in Z for correlation orbitals"
         read  (5, *) delz
      end if
      write (stderr_unit, *) "Orbital Z_eff for hydrogenic orbitals"
      do j = 1, nsubs
         loc = index(j)
         if (set(loc)) cycle
         if (modify) then
            write (stderr_unit, *) "Input new value >"
            read (5, *) sigma(loc)
         end if
         zeff = nuclear_charge - sigma(loc) + delz
         write (stderr_unit, "(I2,A2,F10.2)") orbital_principal_n(loc), orbital_symmetry_labels(loc), zeff
!           ...Calculate radial wavefunctions
         call dcbsrw (orbital_principal_n(loc), orbital_kappa(loc), zeff, orbital_energy_parameters(loc), &
            & orbital_origin_parameters(loc), radial_large_components(:,loc), &
            radial_small_components(:,loc), orbital_last_grid_index(loc))
         set(loc) = .true.
         source(loc) = "Hyd"
      end do
      write (stderr_unit, *)
      return
      end subroutine frmhyd
