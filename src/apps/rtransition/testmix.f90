!***********************************************************************
!                                                                      *
      subroutine testmix
!                                                                      *
!   This routine checks the mixing coefficients                        *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use atomic_calculation_state, only: nelecff, nelecii
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
!-----------------------------------------------
!
!

      write (*, *) " ****************"
      write (*, *) " Entering testmix"
      write (*, *) " ****************"
      write (*, *)
      write (*, *) "Initial state"
      write (*, *) "EVALII", (eavii + evalii(i),i=1,nvecii)
      write (*, *) "NELECII,NCFII,NWII,NVECMXII", nelecii, ncfii, nwii, nvecmxii
      write (*, *) nvecii
      write (*, *) (ivecii(i),i=1,nvecii)
      write (*, *) (iatjpoii(i),iasparii(i),i=1,nvecii)
      write (*, *) ((evecii(i + (j - 1)*ncfii),i=1,ncfii),j=1,nvecii)

      write (*, *) "Final state"
      write (*, *) "EVALFF", (eavff + evalff(i),i=1,nvecff)
      write (*, *) "NELECFF,NCFFF,NWFF,NVECMXFF", nelecff, ncfff, nwff, nvecmxff
      write (*, *) nvecff
      write (*, *) (ivecff(i),i=1,nvecff)
      write (*, *) (iatjpoff(i),iasparff(i),i=1,nvecff)
      write (*, *) ((evecff(i + (j - 1)*ncfff),i=1,ncfff),j=1,nvecff)
      write (*, *)
      write (*, *) " ***************"
      write (*, *) " Leaving testmix"
      write (*, *) " ***************"

      return
      end subroutine testmix
