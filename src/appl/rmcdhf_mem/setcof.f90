!***********************************************************************
!                                                                      *
      subroutine setcof(eol, j)
!                                                                      *
!   This  subroutine  sets  up the coefficients and orbital pointers   *
!   for the direct and exchange potentials for orbital  J .  It also   *
!   sets  up  the  coefficients  and  pointers for the inhomogeneous   *
!   terms arising from off-diagonal I (a,b) integrals.                 *
!                                                                      *
!   Call(s) to: [LIB92]: alloc, dalloc.                                *
!               [RSCF92]: alcsca, dsubrs.                              *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 21 Dec 1992   *
!   Modified by Xinghong He                 Last update: 21 Dec 1997   *
!   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
!      Spin-angular coefficiants have been put on memory               *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use grasp_parameters,    only: orbital_index_radix
      use array_memory_management, only: alloc, dalloc, ralloc
      use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa, orbital_two_j
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use mcp_file_metadata, only: max_mcp_file_multipole_rank
      use mpi_s, only: mpi_process_count, mpi_rank
      use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
          csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use sparse_hamiltonian_storage,           only: stored_matrix_element_count
      use rmcdhf_mem_c,     only: indx, iclmn, coeff, imin2_30, imax2_30
!      USE cons_C,           ONLY: EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dsubrs_i, only: dsubrs
      use fco_i, only: fco
      use gco_i, only: gco
      use read1_mem_i, only: read1_mem
      use read3_mem_i, only: read3_mem
      use read02_mem_i, only: read02_mem
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      logical  :: eol
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
      character(len=6), parameter :: myname = "SETCOF"
      real(real64), parameter :: eps = 1.0e-10_real64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer  :: indexs(4)
      integer :: ndim, nakj, nkjj, ilabel, nwterm, ib, k, nb, ir, nkjib, kmax, &
         kmin, nfile, i, jblock, jblockt, iendcdum, ncoeff, ios, lab, ncontr, &
         ia, ic, ithis, irank, iind, il, iorb, iyo1, iyo2, ifound, loc1, loc2, &
         ithis2, indind, nelmntgg
      real(real64) :: ucfj, sumr, ykab, xkab, sum,  contr
      character :: mcplab*3, idstring*3, msg*128
      integer(int64) :: ncontr_tot, ncontr_tot2
      integer :: iclmndum
!-----------------------------------------------
!** Locals
!     POINTER (PCOEFF,COEFF(1))
!     POINTER (PICLMN,ICLMN(1))
!     POINTER (PINDEX,indx(1))
!
!     POINTER (PNIROW,IROW(1))
!     POINTER (PNTRDA,DA(1))
!     POINTER (PNTRXA,XA(1))
!     POINTER (PNTRYA,YA(1))
!     POINTER (PNTNDA,NDA(1))
!     POINTER (PNTNXA,NXA(1))
!     POINTER (PNTNYA,NYA(1))
!
!
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (pnevblk, nevblk(1))
!     POINTER (pncmaxblk, ncmaxblk(1))
!
!
!     POINTER (pncfpast, ncfpast(1))
!     POINTER (pncminpast, ncminpast(1))
!     POINTER (pnevecpast, nevecpast(1))
!
!     INTEGER, DIMENSION(:), pointer :: nda, nxa, nya
!     INTEGER, DIMENSION(:), pointer :: ncfblk, nevblk, ncmaxblk, &
!                                       ncfpast, ncminpast, nevecpast

!=======================================================================
!   Initializations
!=======================================================================

      ndim = 1
      call alloc (coeff, ndim, "COEFF", "SETCOF" )
      call alloc (iclmn, ndim, "ICLMN", "SETCOF")
      call alloc (indx, ndim, "INDX", "SETCOF")

      ndcof = 0
      nxcof = 0
      nycof = 0

      nakj = orbital_kappa(j)
      nkjj = orbital_two_j(j)
      ucfj = ucf(j)

!=======================================================================
!   Generate YA coefficients that do not require MCP output list.
!   Computation distributed and then collected.
!=======================================================================

      ilabel = 0
      nwterm = key*(key + 1)
      do ib = 1, orbital_count
         ilabel = ilabel + nwterm
         if (ib == j) then
            kmax = nkjj - 1
         else
            kmax = 0
         end if
         do k = 0, kmax, 2

            !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
            sumr = 0.0_real64
            do nb = 1, symmetry_block_count
               do ir = mpi_rank + 1, csf_counts_by_block(nb), mpi_process_count
                  sumr = sumr + dsubrs(eol,ir,ir,nb)*fco(k,ir + csf_offsets_by_block(nb),j,&
                     ib)
               end do
            end do
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            if (ib == j) then
               ykab = 2.0_real64*sumr/ucfj
            else
               ykab = sumr/ucfj
            end if

           !*** The following IF has to be removed since YKAB
           !*** is incomplete in multiprocessor case
            !IF (ABS (YKAB) .GT. EPS) THEN
            nycof = nycof + 1
            if (nycof > nydim ) then
               if (nydim > 0) then
                 nydim = 2*nydim
!cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                 call ralloc(nya, nydim, "NYA", "SETCOF")
                 call ralloc(ya, nydim, "YA", "SETCOF")
!                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
               else
!cjb similar here
                 nydim = 64
                 call alloc(nya, nydim, "NYA", "SETCOF")
                 call alloc(ya, nydim, "YA", "SETCOF")
!                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
!                 NYDIM = NYCOF
!cjb
               end if
            end if
            ya(nycof) = ykab
            nya(nycof) = k + ilabel
            !ENDIF

         end do
      end do

!=======================================================================
!   Generate XA coefficients that do not require MCP output list
!   Computation distributed and then collected.
!=======================================================================

      ilabel = key*j
      nwterm = key*key*(key + 1)
      do ib = 1, orbital_count
         ilabel = ilabel + nwterm
         if (ib == j) cycle
         nkjib = orbital_two_j(ib)
         if (nakj*orbital_kappa(ib) > 0) then
            kmin = abs((nkjj - nkjib)/2)
         else
            kmin = abs((nkjj - nkjib)/2) + 1
         end if
         kmax = (nkjj + nkjib)/2

         do k = kmin, kmax, 2

            !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
            sumr = 0.0_real64
            do nb = 1, symmetry_block_count
               do ir = mpi_rank + 1, csf_counts_by_block(nb), mpi_process_count
                  sumr = sumr + dsubrs(eol,ir,ir,nb)*gco(k,ir + csf_offsets_by_block(nb),j,&
                     ib)
               end do
            end do
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            xkab = sumr/ucfj

           !*** The following IF has to be removed since YKAB
           !*** is incomplete in multiprocessor case
            !IF (ABS (XKAB) .GT. EPS) THEN
            nxcof = nxcof + 1
            if (nxcof > nxdim) then
               if (nxdim > 0) then
                 nxdim = 2*nxdim
!cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                 call ralloc (xa, nxdim, "XA", "SETCOF")
                 call ralloc (nxa, nxdim, "NXA", "SETCOF")
!                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
               else
                  nxdim = 64
                 call alloc (xa, nxdim, "XA", "SETCOF")
                 call alloc (nxa, nxdim, "NXA", "SETCOF")
!                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
!                  NXDIM = NXCOF
!cjb
               end if
            end if
            xa(nxcof) = xkab
            nxa(nxcof) = k + ilabel
            !ENDIF
         end do
      end do

!=======================================================================
!   Generate DA coefficients; these arise from the one-electron
!   integrals
!=======================================================================

      ncontr_tot = 0
      ncontr_tot2 = 0
      do jblock = 1, symmetry_block_count

         !*** Read in IROW from file mcp.30 ***
         nelmntgg = imax2_30(jblock) - imin2_30(jblock)
         stored_matrix_element_count = int(nelmntgg,kind=int64)
         call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "SETCOF")
         call read02_mem(jblock)

         !*** Loop over labels having non-zero coefficients
         !*** it exits when no more labels for the block

         l123: do while(.true.)
            ncontr_tot2 = ncontr_tot2 + 1
            call read1_mem(31,ncontr_tot2,lab, ncontr)

            ! 0, 0 marks the end of a block. This is the normal exit

            if (lab==0 .and. ncontr==0) then
               call dalloc (matrix_partner_indices, "IROW", "SETCOF")
               exit                              ! Actually to next block
            end if

            !*** Decode the labels of I(ab) ***
            ia = mod(lab,key)
            ib = lab/key

            ! At least one orbital should be J in order to have
            ! non-zero value; otherwise, goto next label.

            if (ia/=j .and. ib/=j) then
               ncontr_tot = ncontr_tot + ncontr  ! No contributions from this integral; skip
               cycle                             ! to next label
            end if

            if (ncontr > ndim) then
               call dalloc (coeff, "COEFF", "SETCOF")
               call dalloc (iclmn, "ICLMN", "SETCOF")
               call dalloc (indx, "INDX", "SETCOF")
               ndim = ncontr
               call alloc (coeff, ndim, "COEFF", "SETCOF" )
               call alloc (iclmn, ndim, "ICLMN", "SETCOF")
               call alloc (indx, ndim, "INDX", "SETCOF")
            end if

            ! Read the column index, the sparse matrix index, and the
            ! coefficient for all contributions from this integral

            call read3_mem(31,ncontr_tot,ncontr)
            ncontr_tot = ncontr_tot + ncontr

            ! Add up all the contributions from this integral;
            ! off-diagonal contributions have double the weight

            sum = 0.0_real64
            do i = 1, ncontr
               ir = matrix_partner_indices(indx(i))
               ic = iclmn(i)

               contr = dsubrs(eol,ir,ic,jblock)*coeff(i)

               if (ir /= ic) contr = contr + contr
               sum = sum + contr
            end do
            sum = 0.5_real64*sum/ucfj

            ! Put coefficients in the list. Since there is always
            ! (almost) some repetition from different blocks, a check
            ! and merge is performed. This will significantly reduce
            ! the NDCOF and thus the number of calls to YZK later.

            !*** Find the right counting parameter ***
            if (ia == j) then
               ithis = ib
            else
               ithis = ia
            end if

            !*** Check it against the previously recorded ***
            if (jblock > 1) then
               do i = 1, ndcof
                  if (nda(i) /= ithis) cycle     ! found, add the value
                  da(i) = da(i) + sum
                  cycle  l123
               end do
            end if

            !*** Not found in the record, add an item ***
            ndcof = ndcof + 1
            if (ndcof > nddim) then
               if (nddim > 0) then
                  nddim = 2*nddim
                  call ralloc (da, nddim, "DA", "SETCOF")
                  call ralloc (nda, nddim, "NDA", "SETCOF")
!                  CALL RALLOC (DA, NDCOF, 'DA', 'SETCOF')
!                  CALL RALLOC (NDA, NDCOF, 'NDA', 'SETCOF')
               else
                  nddim = 64
                 call alloc (da, nddim, "DA", "SETCOF")
                  call alloc (nda, nddim, "NDA", "SETCOF")
!                  CALL ALLOC (DA, NDCOF, 'DA', 'SETCOF')
!                  CALL ALLOC (NDA, NDCOF, 'NDA', 'SETCOF')
!                  NDDIM = NDCOF
               end if
            end if
            da(ndcof) = sum                      ! print*, DA(NDCOF), ndcof, myid, 'myid'
            nda(ndcof) = ithis

         end do l123                             ! For labels
      end do                                     ! For blocks

!=======================================================================
!   Generate YA and XA coefficients; these arise from the two-electron
!   integrals
!=======================================================================

      do nfile = 32, 32 + max_mcp_file_multipole_rank

!=======================================================================
!   Loop over blocks again, this time, for V-coefficients
!=======================================================================

         ncontr_tot = 0
         ncontr_tot2 = 0
         do jblock = 1, symmetry_block_count
!           ...Read in IROW from file mcp.30
            nelmntgg = imax2_30(jblock) - imin2_30(jblock)
            stored_matrix_element_count = int(nelmntgg,kind=int64)
            call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "SETCOF")
            call read02_mem(jblock)

            k = nfile - 32                       ! multipolarity of the integral

!=======================================================================
!   Attempt to read another block of data
!=======================================================================

  999       continue
            ncontr_tot2 = ncontr_tot2 + 1
            call read1_mem(nfile,ncontr_tot2,lab, ncontr)
!
            if (lab==0 .and. ncontr==0) then
               call dalloc (matrix_partner_indices, "IROW", "SETCOF")
               cycle
            end if
            !***                       k
            !*** Decode the labels of R (abcd)
            indexs(4) = mod(lab,key)
            lab = lab/key
            indexs(2) = mod(lab,key)
            lab = lab/key
            indexs(3) = mod(lab,key)
            indexs(1) = lab/key

            !*** Determine the number of indices that match
            irank = 0
            irank = irank + count(indexs==j)

            if (irank == 0) then
               ncontr_tot = ncontr_tot + ncontr
               goto 999
            end if

            !*** At least one subshell index matches; allocate storage
            !*** for reading in the rest of this block
            if (ncontr > ndim) then
               call dalloc (coeff, "COEFF", "SETCOF")
               call dalloc (iclmn, "ICLMN", "SETCOF")
               call dalloc (indx, "INDX", "SETCOF")
               ndim = ncontr
               call alloc (coeff, ndim, "COEFF", "SETCOF")
               call alloc (iclmn, ndim, "ICLMN", "SETCOF")
               call alloc (indx, ndim, "INDX", "SETCOF")
            end if

            !*** Read column index, sparse matrix index, and
            !*** coefficient for all contributions from this integral
            call read3_mem(nfile,ncontr_tot,ncontr)
            ncontr_tot = ncontr_tot + ncontr

            !*** Add up all the contributions from this integral;
            !*** off-diagonal contributions have double the weight
            sum = 0.0_real64
            do i = 1, ncontr
               ir = matrix_partner_indices(indx(i))
               ic = iclmn(i)
               contr = dsubrs(eol,ir,ic,jblock)*coeff(i)
               if (ir /= ic) contr = contr + contr
               sum = sum + contr
            end do
            sum = 0.5_real64*sum/ucfj

            select case (irank)
            case (1)

!=======================================================================
!   One matching index: exchange potential contribution
!=======================================================================

               !*** Similar to DA, find ithis ***
               ithis = -911                      ! initialize to an impossible value
                              ! though not necessary
               do iind = 1, 4
                  if (indexs(iind) /= j) cycle   ! at least one
                  il = iind + 2
                  if (il > 4) il = il - 4
                  iorb = indexs(il)
                  il = iind + 1
                  if (il > 4) il = il - 4
                  iyo1 = indexs(il)
                  il = iind + 3
                  if (il > 4) il = il - 4
                  iyo2 = indexs(il)
                  ithis = ((iorb*key + iyo2)*key + iyo1)*key + k
                  exit
               end do

               if (ithis == (-911)) stop "ithis .EQ. -911"

               !*** Check ithis against the previously recorded ***
               if (jblock > 1) then
                  do i = 1, nxcof
                     if (nxa(i) /= ithis) cycle
                     xa(i) = xa(i) + sum
                     goto 999
                  end do
               end if

               !*** Not found in records, add an item ***
               nxcof = nxcof + 1
               if (nxcof > nxdim) then
                  if (nxdim > 0) then
                    nxdim = 2*nxdim
!cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                    call ralloc (xa, nxdim, "XA", "SETCOF")
                    call ralloc (nxa, nxdim, "NXA", "SETCOF")
!                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                 else
                    nxdim = 64
                    call alloc (xa, nxdim, "XA", "SETCOF")
                    call alloc (nxa, nxdim, "NXA", "SETCOF")
!                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
!                  NXDIM = NXCOF
!cjb
                  end if
               end if
               xa(nxcof) = sum
               nxa(nxcof) = ithis

            case (2)

!=======================================================================
!   Two matching indices: either direct or exchange potential
!   contribution
!=======================================================================

               ifound = 0
               do iind = 1, 4
                  if (indexs(iind) /= j) cycle
                  if (ifound == 0) then
                     loc1 = iind
                     ifound = ifound + 1
                  else if (ifound == 1) then
                     loc2 = iind
                     exit
                  end if
               end do

               if (loc2 - loc1 == 2) then
!
!   Direct contribution
!
                  !*** Find ithis ***
                  il = loc1 + 3
                  if (il > 4) il = il - 4
                  iyo2 = indexs(il)
                  il = loc1 + 1
                  if (il > 4) il = il - 4
                  iyo1 = indexs(il)
                  ithis = (iyo2*key + iyo1)*key + k

                  !*** Check it against the previously recorded ***
                  if (jblock > 1) then
                     do i = 1, nycof
                        if (nya(i) /= ithis) cycle
                        ya(i) = ya(i) + sum + sum
                        goto 999
                     end do
                  end if

                  !*** Not found, add an item ***
                  nycof = nycof + 1
                  if (nycof > nydim) then
                     if (nydim > 0) then
                       nydim = 2*nydim
!cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                       call ralloc(nya, nydim, "NYA", "SETCOF")
                       call ralloc(ya, nydim, "YA", "SETCOF")
!                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
                     else
                       nydim = 64
                       call alloc(nya, nydim, "NYA", "SETCOF")
                       call alloc(ya, nydim, "YA", "SETCOF")
!                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
!                 NYDIM = NYCOF
!cjb
                     end if
                  end if
                  ya(nycof) = sum + sum
                  nya(nycof) = ithis

               else
!
!   Exchange contribution
!
                  !*** Find ithis ***
                  il = loc1 + 2
                  if (il > 4) il = il - 4
                  iorb = indexs(il)
                  il = loc1 + 1
                  if (il > 4) il = il - 4
                  iyo1 = indexs(il)
                  il = loc1 + 3
                  if (il > 4) il = il - 4
                  iyo2 = indexs(il)
                  ithis = ((iorb*key + iyo2)*key + iyo1)*key + k

                  !*** Check it against the previously recorded ***
                  if (jblock > 1) then
                     do i = 1, nxcof
                        if (nxa(i) /= ithis) cycle
                        xa(i) = xa(i) + sum + sum
                        goto 999
                     end do
                  end if

                  !*** Not found, add an item ***
                  nxcof = nxcof + 1
                  if (nxcof > nxdim) then
                     if (nxdim > 0) then
                       nxdim = 2*nxdim
!cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                       call ralloc (xa, nxdim, "XA", "SETCOF")
                       call ralloc (nxa, nxdim, "NXA", "SETCOF")
!                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                     else
                       nxdim = 64
                       call alloc (xa, nxdim, "XA", "SETCOF")
                       call alloc (nxa, nxdim, "NXA", "SETCOF")
!                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
!                  NXDIM = NXCOF
!cjb
                     end if
                  end if
                  xa(nxcof) = sum + sum
                  nxa(nxcof) = ithis

               end if

            case (3)

!=======================================================================
!   Three matching indices: direct and exchange potential contributions
!=======================================================================

               !*** Find ithis AND ithis2
               ithis = -911
               ithis2 = -911
               do iind = 1, 4
                  if (indexs(iind) == j) cycle
                  indind = indexs(iind)
                  iyo2 = indind
                  iyo1 = j
                  ithis = (iyo2*key + iyo1)*key + k
                  iorb = indind
                  iyo1 = j
                  iyo2 = j
                  ithis2 = ((iorb*key + iyo2)*key + iyo1)*key + k
               end do

               if (ithis==(-911) .or. ithis2==(-911)) stop "ithis2"

               !*** Check the previously recorded for YA
               if (jblock > 1) then
                  do i = 1, nycof
                     if (nya(i) /= ithis) cycle
                     ya(i) = ya(i) + sum + sum
                     goto 456
                  end do
               end if

               ! Not found, add an item ***
               nycof = nycof + 1
               if (nycof > nydim) then
                  if (nydim > 0) then
                    nydim = 2*nydim
                    call ralloc (nya, nydim, "NYA", "SETCOF")
                    call ralloc (ya, nydim, "YA", "SETCOF")
!                     CALL RALLOC (NYA, NYCOF, 'NYA', 'SETCOF')
!                     CALL RALLOC (YA, NYCOF, 'YA', 'SETCOF')
                  else
                    nydim = 64
                     call alloc (nya, nydim, "NYA", "SETCOF")
                     call alloc (ya, nydim, "YA", "SETCOF")
!                     CALL ALLOC (NYA, NYCOF, 'NYA', 'SETCOF')
!                     CALL ALLOC (YA, NYCOF, 'YA', 'SETCOF')
!                     NYDIM = NYCOF
                  end if
               end if
               ya(nycof) = sum + sum
               nya(nycof) = ithis

  456          continue

               !*** Check the previously recorded for XA
               if (jblock > 1) then
                  do i = 1, nxcof
                     if (nxa(i) /= ithis2) cycle
                     xa(i) = xa(i) + sum
                     goto 999
                  end do
               end if

               ! Not found, add an item ***
               nxcof = nxcof + 1
               if (nxcof > nxdim) then
                  if (nxdim > 0) then
                    nxdim = 2*nxdim
!cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                   call ralloc (xa, nxdim, "XA", "SETCOF")
                   call ralloc (nxa, nxdim, "NXA", "SETCOF")
!                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                 else
                   call alloc (xa, nxdim, "XA", "SETCOF")
                   call alloc (nxa, nxdim, "NXA", "SETCOF")
!                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
!                  NXDIM = NXCOF
!cjb
                  end if
               end if
               xa(nxcof) = sum
               nxa(nxcof) = ithis2

            case (4)

!=======================================================================
!   Four matching indices: direct potential contribution
!=======================================================================

               !*** Find ithis AND ithis2
               iyo2 = j
               iyo1 = j
               ithis = (iyo2*key + iyo1)*key + k

               !*** Check the previously recorded for YA
               if (jblock > 1) then
                  do i = 1, nycof
                     if (nya(i) /= ithis) cycle
                     ya(i) = ya(i) + 4.0_real64*sum
                     goto 999
                  end do
               end if

               ! Not found, add an item ***
               nycof = nycof + 1
               if (nycof > nydim) then
                  if (nydim > 0) then
                    nydim = 2*nydim
!cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                    call ralloc(nya, nydim, "NYA", "SETCOF")
                    call ralloc(ya, nydim, "YA", "SETCOF")
!                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
                 else
!cjb similar here
                   nydim = 64
                   call alloc(nya, nydim, "NYA", "SETCOF")
                   call alloc(ya, nydim, "YA", "SETCOF")
!                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
!                 NYDIM = NYCOF
!cjb
                  end if
               end if
               ya(nycof) = 4.0_real64*sum
               nya(nycof) = ithis

            end select

            goto 999
         end do                                  ! loop for V-Coefficients
      end do

!=======================================================================
!   Deallocate storage for arrays local to this routine
!=======================================================================
      call dalloc (coeff, "COEFF", "SETCOF")
      call dalloc (iclmn, "ICLMN", "SETCOF")
      call dalloc (indx, "INDX", "SETCOF")

      return
      end subroutine setcof
