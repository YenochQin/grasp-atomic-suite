!***********************************************************************
!                                                                      *
      subroutine setham(jblock, myid, nprocs)
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use parameter_def,    only: orbital_index_radix
      use memory_man, only: alloc, dalloc, double, long
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      use hmat_c, only: double, long, stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use mcpa_c, only: double, max_mcp_file_multipole_rank
      use orb_c, only: csf_count, double, orbital_count, orbital_kappa, orbital_two_j
      use pos_c, only: double, total_csf_count, mixing_coefficient_storage_size, &
          csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i, only: iq
      use rinti_i, only: rinti
      use fco_i, only: fco
      use slater_i, only: slater
      use gco_i, only: gco
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jblock
      integer , intent(in) :: myid
      integer , intent(in) :: nprocs
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
      character(len=6), parameter :: myname = "SETHAM"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncfpat, ia, ir, itmp, idiag, km, k, k0, ib, nkjia, nkjib, kmin&
         , kmax, ndim, jblockt, ncft, ncoeff, ios, lab, ncontr, iclmndum, i, &
         loc, nfile, id, ic
      integer(int64) :: lod
      real(real64) :: diaa, coef, f0aa, fkaa, f0ab, gkab, tegral
      logical :: set
      character :: mcplab*3
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
      ncfpat = csf_offsets_by_block(jblock)
!=======================================================================
!   Accumulate diagonal terms that do not require MCP coefficients
!=======================================================================

!=======================================================================
!   Piece involving I(a,a) integrals
!=======================================================================

      do ia = 1, orbital_count
         set = .false.
         do ir = myid + 1, csf_count, nprocs
            itmp = iq(ia,ir + ncfpat)
            if (itmp <= 0) cycle
            !*** Occupation number not zero ...
            if (.not.set) then
               diaa = rinti(ia,ia,0)
               set = .true.
            end if
            ! IDIAG = IENDC(IR-1)+1
            idiag = stored_segment_end_offsets(ir)                    ! lower-triangle-by-rows mode
            hamiltonian_values(idiag) = hamiltonian_values(idiag) + itmp*diaa
         end do
      end do

!=======================================================================
!                    0
!   Piece involving F (a,a) integrals
!=======================================================================

      do ia = 1, orbital_count
         set = .false.
         do ir = myid + 1, csf_count, nprocs
            coef = fco(0,ir + ncfpat,ia,ia)
            if (coef == 0.0_real64) cycle
            !*** Angular coefficient not zero ...
            if (.not.set) then
               f0aa = slater(ia,ia,ia,ia,0)
               set = .true.
            end if
            ! IDIAG = IENDC(IR-1)+1
            idiag = stored_segment_end_offsets(ir)
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
      do ia = 1, orbital_count
         k0 = orbital_two_j(ia) - 1
         km = max0(k0,km)
         if (k > k0) cycle
         set = .false.

         do ir = myid + 1, csf_count, nprocs
            coef = fco(k,ir + ncfpat,ia,ia)
            if (coef == 0.0_real64) cycle
            if (.not.set) then
               fkaa = slater(ia,ia,ia,ia,k)
               set = .true.
            end if
            ! IDIAG = IENDC(IR-1)+1
            idiag = stored_segment_end_offsets(ir)
            hamiltonian_values(idiag) = hamiltonian_values(idiag) + coef*fkaa
         end do
      end do
      if (k < km) goto 6

!=======================================================================
!                    0
!   Piece involving F (a,b) integrals
!=======================================================================

      do ia = 1, orbital_count - 1
         do ib = ia + 1, orbital_count
            set = .false.
            do ir = myid + 1, csf_count, nprocs
               coef = fco(0,ir + ncfpat,ia,ib)
               if (coef == 0.0_real64) cycle
               if (.not.set) then
                  f0ab = slater(ia,ib,ia,ib,0)
                  set = .true.
               end if
            ! IDIAG = IENDC(IR-1)+1
               idiag = stored_segment_end_offsets(ir)
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
      do ia = 1, orbital_count - 1
         nkjia = orbital_two_j(ia)
         do ib = ia + 1, orbital_count
            nkjib = orbital_two_j(ib)
            set = .false.
            if (orbital_kappa(ia)*orbital_kappa(ib) > 0) then
               kmin = abs((nkjia - nkjib)/2)
            else
               kmin = abs((nkjia - nkjib)/2) + 1
            end if
            if (mod(k - kmin,2) /= 0) cycle

            kmax = (nkjia + nkjib)/2
            km = max0(kmax,km)
            if (k<kmin .or. k>kmax) cycle

            do ir = myid + 1, csf_count, nprocs
               coef = gco(k,ir + ncfpat,ia,ib)
               if (coef == 0.0_real64) cycle
               if (.not.set) then
                  gkab = slater(ia,ib,ib,ia,k)
                  set = .true.
               end if
               ! IDIAG = IENDC(IR-1)+1
               idiag = stored_segment_end_offsets(ir)
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

      read (31) mcplab, jblockt, ncft, ncoeff

      if (jblockt /= jblock) then
         write (stderr_unit, *) myname, ": blk1=", jblockt, " blk2=", jblock
         error stop
      end if
      if (ncft /= csf_count) then
         write (stderr_unit, *) myname, ": ncf1 = ", ncft, " ncf2 = ", csf_count
         error stop
      end if

!=======================================================================
!   Loop over non-zero labels which have non-zero elements
!=======================================================================

      read (31, iostat=ios) lab, ncontr
      if (ios /= 0) error stop "IOS .NE. 0 when reading LAB, NCONTR"
      do while(lab/=0 .or. ncontr/=0)

         !*** decode the label of I(ab)
         ia = mod(lab,key)
         ib = lab/key

         !*** Compute radial integral I(ab)
         tegral = rinti(ia,ib,0)

         ! Read column index, sparse matrix index, and coefficient
         ! for all contributions from this integral.
         if (ncontr > ndim) then
            call dalloc (coeff, "COEFF", "SETHAM")
            call dalloc (indx, "INDX", "SETHAM")
            ndim = ncontr
            call alloc (coeff, ndim, "COEFF", "SETHAM")
            call alloc (indx, ndim, "INDX", "SETHAM")
         end if
         read (31) (iclmndum,indx(i),coeff(i),i=1,ncontr)

         !*** Store all the contributions from this integral
         do i = 1, ncontr
            loc = indx(i)
            lod = loc
            if (lod > (stored_matrix_element_count)) then
               write (6, *) "  Error in computing 1-e contribution"
               write (6, *) "  LOC = ", loc, "  NELMNT = ", (stored_matrix_element_count)
               error stop
            end if
            hamiltonian_values(loc) = hamiltonian_values(loc) + tegral*coeff(i)
         end do

         read (31, iostat=ios) lab, ncontr
         if (ios == 0) cycle
         error stop "IOS .NE. 0 when reading LAB, NCONTR"
      end do

!=======================================================================
!   Accumulate two-electron terms that require MCP coefficients
!=======================================================================

      do nfile = 32, 32 + max_mcp_file_multipole_rank
         k = nfile - 32

         read (nfile) mcplab, jblockt, ncft, ncoeff

         if (jblockt /= jblock) then
            write (stderr_unit, *) myname, ": blk3=", jblockt, " blk4=", jblock
            error stop
         end if
         if (ncft /= csf_count) then
            write (stderr_unit, *) myname, ": ncf3 = ", ncft, " ncf4 = ", csf_count
            error stop
         end if

!=======================================================================
!   Loop over non-zero labels which have non-zero elements
!=======================================================================

         read (nfile, iostat=ios) lab, ncontr
         if (ios /= 0) error stop "IOS .NE. 0 when reading LAB, NCONTR 2"
         do while(lab/=0 .or. ncontr/=0)

            !                         k
            !*** decode the label of R (abcd)
            id = mod(lab,key)
            lab = lab/key
            ib = mod(lab,key)
            lab = lab/key
            ic = mod(lab,key)
            ia = lab/key

            !*** Compute radial integral
            tegral = slater(ia,ib,ic,id,k)

            ! Read column index, sparse matrix index, and coefficient
            ! for all contributions from this integral.
            if (ncontr > ndim) then
               call dalloc (coeff, "COEFF", "SETHAM")
               call dalloc (indx, "INDX", "SETHAM")
               ndim = ncontr
               call alloc (coeff, ndim, "COEFF", "SETHAM")
               call alloc (indx, ndim, "INDX", "SETHAM")
            end if
            read (nfile) (iclmndum,indx(i),coeff(i),i=1,ncontr)

            !*** Store all the contributions from this integral
            do i = 1, ncontr
               loc = indx(i)
               lod = loc
               if (lod > stored_matrix_element_count) then
                  write (6, *) "  Error in computing 2-e contribution"
                  write (6, *) "  LOC = ", loc, "  NELMNT = ", stored_matrix_element_count
                  error stop
               end if
               hamiltonian_values(loc) = hamiltonian_values(loc) + tegral*coeff(i)
            end do

            read (nfile, iostat=ios) lab, ncontr
            if (ios == 0) cycle
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
