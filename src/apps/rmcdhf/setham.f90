!***********************************************************************
!                                                                      *
      subroutine setham(symmetry_block_index, mpi_rank, mpi_process_count)
      use mcp_file_metadata, only: mcp_input_units
      use, intrinsic :: iso_fortran_env, only: output_unit
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use grasp_parameters,    only: orbital_index_radix
      use array_memory_management, only: alloc, dalloc
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use mcp_file_metadata, only: max_mcp_file_multipole_rank
      use orbital_csf_basis, only: csf_count, orbital_count, orbital_kappa, orbital_two_j
      use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
          csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: iq
      use one_electron_radial_integral, only: rinti
      use fco_i, only: fco
      use slater_radial_integral, only: slater
      use gco_i, only: gco
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: mpi_rank, mpi_process_count
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
      character(len=6), parameter :: routine_name = "SETHAM"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncfpat, orbital_index, csf_index, itmp, idiag, km, k, k0, ib, nkjia, nkjib, kmin&
         , kmax, ndim, jblockt, ncft, ncoeff, io_status, lab, ncontr, iclmndum, i, &
         loc, nfile, id, ic
      integer(int64) :: lod      !  INTEGER*8 version of LOC
      real(real64) :: diaa, coef, f0aa, fkaa, f0ab, gkab, tegral
      logical :: set
      character :: mcp_channel_label*3
      real(real64), pointer :: coeff(:)
      integer, pointer :: indx(:)
!-----------------------------------------------
!
!     !*** needs EMT, IENDC.
!     POINTER (PNTEMT,EMT(1))
!     POINTER (PIENDC,IENDC(0:*))
!     POINTER (PNIROW,IROW(1))
!
!     !*** needs ncf, nw
!     POINTER (PNTRIQ,RIQDUMMY)
!
!
!     POINTER (pncfpast, ncfpast(1))
!     POINTER (pncminpast, ncminpast(1))
!     POINTER (pnevecpast, nevecpast(1))
!
!     POINTER (PCOEFF,COEFF(1))
!     !POINTER (PICLMN,ICLMNdum)
!     POINTER (PINDX,INDX(1))
!
!-----------------------------------------------------------------------
      ncfpat = csf_offsets_by_block(symmetry_block_index)
!=======================================================================
!   Accumulate diagonal terms that do not require MCP coefficients
!=======================================================================

!=======================================================================
!   Piece involving I(a,a) integrals
!=======================================================================

      do orbital_index = 1, orbital_count
         set = .false.
         do csf_index = mpi_rank + 1, csf_count, mpi_process_count
            itmp = iq(orbital_index,csf_index + ncfpat)
            if (itmp <= 0) cycle
            !*** Occupation number not zero ...
            if (.not.set) then
               diaa = rinti(orbital_index,orbital_index,0)
               set = .true.
            end if
            ! IDIAG = IENDC(IR-1)+1
            idiag = stored_segment_end_offsets(csf_index)                    ! lower-triangle-by-rows mode
            hamiltonian_values(idiag) = hamiltonian_values(idiag) + itmp*diaa
         end do
      end do

!=======================================================================
!                    0
!   Piece involving F (a,a) integrals
!=======================================================================

      do orbital_index = 1, orbital_count
         set = .false.
         do csf_index = mpi_rank + 1, csf_count, mpi_process_count
            coef = fco(0,csf_index + ncfpat,orbital_index,orbital_index)
            if (coef == 0.0_real64) cycle
            !*** Angular coefficient not zero ...
            if (.not.set) then
               f0aa = slater(orbital_index,orbital_index,orbital_index,orbital_index,0)
               set = .true.
            end if
            ! IDIAG = IENDC(IR-1)+1
            idiag = stored_segment_end_offsets(csf_index)
            hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*f0aa
         end do
      end do

!=======================================================================
!                    k
!   Piece involving F (a,a) integrals
!=======================================================================

      km = 0
      k = 0
    6 continue
      k = k + 2
      do orbital_index = 1, orbital_count
         k0 = orbital_two_j(orbital_index) - 1
         km = max0(k0,km)
         if (k > k0) cycle
         set = .false.

         do csf_index = mpi_rank + 1, csf_count, mpi_process_count
            coef = fco(k,csf_index + ncfpat,orbital_index,orbital_index)
            if (coef == 0.0_real64) cycle
            if (.not.set) then
               fkaa = slater(orbital_index,orbital_index,orbital_index,orbital_index,k)
               set = .true.
            end if
            ! IDIAG = IENDC(IR-1)+1
            idiag = stored_segment_end_offsets(csf_index)
            hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*fkaa
         end do
      end do
      if (k < km) goto 6

!=======================================================================
!                    0
!   Piece involving F (a,b) integrals
!=======================================================================

      do orbital_index = 1, orbital_count - 1
         do ib = orbital_index + 1, orbital_count
            set = .false.
            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
               coef = fco(0,csf_index + ncfpat,orbital_index,ib)
               if (coef == 0.0_real64) cycle
               if (.not.set) then
                  f0ab = slater(orbital_index,ib,orbital_index,ib,0)
                  set = .true.
               end if
            ! IDIAG = IENDC(IR-1)+1
               idiag = stored_segment_end_offsets(csf_index)
               hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*f0ab
            end do
         end do
      end do

!=======================================================================
!                    k
!   Piece involving G (a,b) integrals
!=======================================================================

      km = 0
      k = -1
   12 continue
      k = k + 1
      do orbital_index = 1, orbital_count - 1
         nkjia = orbital_two_j(orbital_index)
         do ib = orbital_index + 1, orbital_count
            nkjib = orbital_two_j(ib)
            set = .false.
            if (orbital_kappa(orbital_index)*orbital_kappa(ib) > 0) then
               kmin = abs((nkjia - nkjib)/2)
            else
               kmin = abs((nkjia - nkjib)/2) + 1
            end if
            if (mod(k - kmin,2) /= 0) cycle

            kmax = (nkjia + nkjib)/2
            km = max0(kmax,km)
            if (k<kmin .or. k>kmax) cycle

            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
               coef = gco(k,csf_index + ncfpat,orbital_index,ib)
               if (coef == 0.0_real64) cycle
               if (.not.set) then
                  gkab = slater(orbital_index,ib,ib,orbital_index,k)
                  set = .true.
               end if
               ! IDIAG = IENDC(IR-1)+1
               idiag = stored_segment_end_offsets(csf_index)
               hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*gkab
            end do
         end do
      end do
      if (k < km) goto 12

!=======================================================================
!   Local storage for reading mcpXXX files
!=======================================================================

      ndim = 1
      call alloc (coeff, ndim, "COEFF", "SETHAM" )
      call alloc (indx, ndim, "INDX", "SETHAM")

!=======================================================================
!   Accumulate one-electron terms that require MCP coefficients
!=======================================================================

      read (mcp_input_units(31)) mcp_channel_label, jblockt, ncft, ncoeff

      if (jblockt /= symmetry_block_index) then
         write (stderr_unit, *) routine_name, ": blk1=", jblockt, " blk2=", symmetry_block_index
         error stop
      end if
      if (ncft /= csf_count) then
         write (stderr_unit, *) routine_name, ": ncf1 = ", ncft, " ncf2 = ", csf_count
         error stop
      end if

!=======================================================================
!   Loop over non-zero labels which have non-zero elements
!=======================================================================

      read (mcp_input_units(31), iostat=io_status) lab, ncontr
      if (io_status /= 0) error stop "IOS .NE. 0 when reading LAB, NCONTR"
      do while(lab/=0 .or. ncontr/=0)

         !*** decode the label of I(ab)
         orbital_index = mod(lab,orbital_label_radix)
         ib = lab/orbital_label_radix

         !*** Compute radial integral I(ab)
         tegral = rinti(orbital_index,ib,0)

         ! Read column index, sparse matrix index, and coefficient
         ! for all contributions from this integral.
         if (ncontr > ndim) then
            call dalloc (coeff, "COEFF", "SETHAM")
            call dalloc (indx, "INDX", "SETHAM")
            ndim = ncontr
            call alloc (coeff, ndim, "COEFF", "SETHAM")
            call alloc (indx, ndim, "INDX", "SETHAM")
         end if
         read (mcp_input_units(31)) (iclmndum,indx(i),coeff(i),i=1,ncontr)

         !*** Store all the contributions from this integral
         do i = 1, ncontr
            loc = indx(i)
            lod = loc   !  (Convert type)
            if (lod > stored_matrix_element_count) then
               write (output_unit, *) "  Error in computing 1-e contribution"
               write (output_unit, *) "  LOC = ", loc , "  NELMNT = ", stored_matrix_element_count
               error stop
            end if
            hamiltonian_values(loc) = hamiltonian_values(loc) + tegral*coeff(i)
         end do

         read (mcp_input_units(31), iostat=io_status) lab, ncontr
         if (io_status == 0) cycle
         error stop "IOS .NE. 0 when reading LAB, NCONTR"
      end do

!=======================================================================
!   Accumulate two-electron terms that require MCP coefficients
!=======================================================================

      do nfile = 32, 32 + max_mcp_file_multipole_rank
         k = nfile - 32

         read (mcp_input_units(nfile)) mcp_channel_label, jblockt, ncft, ncoeff

         if (jblockt /= symmetry_block_index) then
            write (stderr_unit, *) routine_name, ": blk3=", jblockt, " blk4=", symmetry_block_index
            error stop
         end if
         if (ncft /= csf_count) then
            write (stderr_unit, *) routine_name, ": ncf3 = ", ncft, " ncf4 = ", csf_count
            error stop
         end if

!=======================================================================
!   Loop over non-zero labels which have non-zero elements
!=======================================================================

         read (mcp_input_units(nfile), iostat=io_status) lab, ncontr
         if (io_status /= 0) error stop "IOS .NE. 0 when reading LAB, NCONTR 2"
         do while(lab/=0 .or. ncontr/=0)

            !                         k
            !*** decode the label of R (abcd)
            id = mod(lab,orbital_label_radix)
            lab = lab/orbital_label_radix
            ib = mod(lab,orbital_label_radix)
            lab = lab/orbital_label_radix
            ic = mod(lab,orbital_label_radix)
            orbital_index = lab/orbital_label_radix

            !*** Compute radial integral
            tegral = slater(orbital_index,ib,ic,id,k)

            ! Read column index, sparse matrix index, and coefficient
            ! for all contributions from this integral.
            if (ncontr > ndim) then
               call dalloc (coeff, "COEFF", "SETHAM")
               call dalloc (indx, "INDX", "SETHAM")
               ndim = ncontr
               call alloc (coeff, ndim, "COEFF", "SETHAM")
               call alloc (indx, ndim, "INDX", "SETHAM")
            end if
            read (mcp_input_units(nfile)) (iclmndum,indx(i),coeff(i),i=1,ncontr)

            !*** Store all the contributions from this integral
            do i = 1, ncontr
               loc = indx(i)
               lod = loc
               if (lod > stored_matrix_element_count) then   !! NELMENT (from Hmat_C) is INT*8
                  write (output_unit, *) "  Error in computing 2-e contribution"
                  write (output_unit, *) "  LOC= ", loc , "  NELMNT = ", stored_matrix_element_count
                  error stop
               end if
               hamiltonian_values(loc) = hamiltonian_values(loc) + tegral*coeff(i)
            end do

            read (mcp_input_units(nfile), iostat=io_status) lab, ncontr
            if (io_status == 0) cycle
            error stop "IOS .NE. 0 when reading LAB, NCONTR 2"
         end do
      end do

!=======================================================================
!   Deallocate local storage
!=======================================================================

      call dalloc (coeff, "COEFF", "SETHAM")
      call dalloc (indx, "INDX", "SETHAM")

      return
      end subroutine setham
