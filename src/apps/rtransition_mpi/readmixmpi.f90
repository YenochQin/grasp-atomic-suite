!***********************************************************************
!                                                                      *
      subroutine readmix(name, inpci, init)
      use transition_matrix_workspace, only: final_mixing_unit, initial_mixing_unit
!                                                                      *
!   Open and read the mixing coefficent files                          *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use array_memory_management, only: alloc
      use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use orbital_csf_basis, only: max_orbitals
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: inpci
      integer , intent(in) :: init
      character  :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ib, iatjp, iaspa, i, j
      character :: mixing_file_signature*6
!-----------------------------------------------
!
      if (init == 1) then
!
!   Read the initial state mixing file
!
         read (initial_mixing_unit) ib, ncfii, nvecii, iatjp, iaspa
         call alloc (evalii, nvecii, "EVALII", "READMIX")
         call alloc (evecii, ncfii*nvecii, "EVECII", "READMIX")
         call alloc (ivecii, nvecii,"IVECII", "READMIX" )
         call alloc (iatjpoii, nvecii, "IATJPOII", "READMIX")
         call alloc (iasparii, nvecii, "ISPARII", "READMIX")
         read (initial_mixing_unit) (ivecii(i),i=1,nvecii)
         iatjpoii(:nvecii) = iatjp
         iasparii(:nvecii) = iaspa
         read (initial_mixing_unit) eavii, (evalii(i),i=1,nvecii)

         read (initial_mixing_unit) ((evecii(i + (j - 1)*ncfii),i=1,ncfii),j=1,nvecii)
!
!     CLOSE(68)
!
      else
!
!   Read the final state mixing file
!

         read (final_mixing_unit) ib, ncfff, nvecff, iatjp, iaspa
         call alloc (evalff, nvecff, "EVALFF", "READMIX")
         call alloc (evecff, ncfff*nvecff, "EVECFF", "READMIX")
         call alloc (ivecff, nvecff,"IVECFF", "READMIX" )
         call alloc (iatjpoff, nvecff, "IATJPOFF", "READMIX")
         call alloc (iasparff, nvecff, "ISPARFF", "READMIX")
         read (final_mixing_unit) (ivecff(i),i=1,nvecff)
         iatjpoff(:nvecff) = iatjp
         iasparff(:nvecff) = iaspa
         read (final_mixing_unit) eavff, (evalff(i),i=1,nvecff)

         read (final_mixing_unit) ((evecff(i + (j - 1)*ncfff),i=1,ncfff),j=1,nvecff)
!
!   Close the initial state mixing  file
!
!      CLOSE(78)
      end if

      return
      end subroutine readmix
