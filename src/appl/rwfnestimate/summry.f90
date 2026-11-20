

!***********************************************************************
!                                                                      *
      subroutine summry(nunit)
!                                                                      *
!   Prints  a summary of the complete  list of subshell radial wave-   *
!   functions on NUNIT.                                                *
!                                                                      *
!   Call(s) to: [LIB92]:                                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Modified by Yanting Li, output <r> for each orbital, 19 Aug 2022
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_energy_parameters, orbital_origin_exponents, &
          orbital_principal_n, orbital_symmetry_labels
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use initial_orbital_source, only: source
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use weighted_radial_overlap, only: rint
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nunit
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, lenth
!-----------------------------------------------
!
!
      write (nunit, '("Shell",6x,"e",11x,"p0",8x,"gamma",8x,"<r>",6x,"MTP","  SRC"/)')
!
      do i = 1, orbital_count
         lenth = len_trim(source(i))
!--------Output P(2) and Q(2) for each orbital
!         WRITE (NUNIT, 301) NP(I), NH(I), E(I), PZ(I), GAMA(I), PF(2,I), QF(2,I&
!            ), MF(I), SOURCE(I)(1:LENTH)
!--------Output <r> (RINT(I,I,1)) for each orbital
         write (nunit, '(1x,i2,a2,4d12.4,i5,2x,a3)') &
            orbital_principal_n(i), orbital_symmetry_labels(i), orbital_energy_parameters(i), &
            & orbital_origin_parameters(i), orbital_origin_exponents(i), rint(i,i,1)&
             , orbital_last_grid_index(i), source(i)(1:lenth)
!         WRITE (NUNIT,302) SOURCE(I)(1:LENTH)
      end do
!
      return
!
!-Format for P(2) and Q(2)
!  300 FORMAT('Shell',6X,'e',11X,'p0',8X,'gamma',8X,'P(2)',7X,'Q(2)',6X,'MTP',&
!         '  SRC'/)
!  301 FORMAT(1X,I2,A2,5D12.4,I5,2X,A3)
!-Format for <r>
      return
!  302 FORMAT ('      Source: ',A)
!
      end subroutine summry
