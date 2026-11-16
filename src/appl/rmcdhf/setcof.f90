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
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:21:02   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double, long
      use parameter_def,    only: keyorb
      use memory_man
      use orb_c
      use hblock_c
      use hmat_c
      use iounit_c
      use mcpa_c
      use mpi_s
      use pos_c
      use scf_c
!      USE cons_C,           ONLY: EPS
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dsubrs_i
      use fco_i
      use gco_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
      logical  :: eol
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
      character*6, parameter :: myname = 'SETCOF'
      real, parameter :: eps = 1.0d-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(4) :: indexs
      integer :: ndim, nakj, nkjj, ilabel, nwterm, ib, k, nb, ir, nkjib, kmax, &
         kmin, nfile, i, jblock, jblockt, iendcdum, ncoeff, ios, lab, ncontr, &
         ia, ic, ithis, irank, iind, il, iorb, iyo1, iyo2, ifound, loc1, loc2, &
         ithis2, indind, nelmntgg
      real(double) :: ucfj, sumr, ykab, xkab, sum,  contr
      character :: mcplab*3, idstring*3, msg*128
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
      real(double), dimension(:), pointer :: coeff
      integer, dimension(:), pointer :: iclmn, indx
!     INTEGER, DIMENSION(:), pointer :: nda, nxa, nya
!     INTEGER, DIMENSION(:), pointer :: ncfblk, nevblk, ncmaxblk, &
!                                       ncfpast, ncminpast, nevecpast

!=======================================================================
!   Initializations
!=======================================================================

      ndim = 1
      call alloc (coeff, ndim, 'COEFF', 'SETCOF' )
      call alloc (iclmn, ndim, 'ICLMN', 'SETCOF')
      call alloc (indx, ndim, 'INDX', 'SETCOF')

      ndcof = 0
      nxcof = 0
      nycof = 0

      nakj = nak(j)
      nkjj = nkj(j)
      ucfj = ucf(j)

!=======================================================================
!   Generate YA coefficients that do not require MCP output list.
!   Computation distributed and then collected.
!=======================================================================

      ilabel = 0
      nwterm = key*(key + 1)
      do ib = 1, nw
         ilabel = ilabel + nwterm
         if (ib == j) then
            kmax = nkjj - 1
         else
            kmax = 0
         endif
         do k = 0, kmax, 2

            !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
            sumr = 0.d0
            do nb = 1, nblock
               do ir = myid + 1, ncfblk(nb), nprocs
                  sumr = sumr + dsubrs(eol,ir,ir,nb)*fco(k,ir + ncfpast(nb),j,&
                     ib)
               end do
            end do
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            if (ib == j) then
               ykab = 2.d0*sumr/ucfj
            else
               ykab = sumr/ucfj
            endif

           !*** The following IF has to be removed since YKAB
           !*** is incomplete in multiprocessor case
            !IF (ABS (YKAB) .GT. EPS) THEN
            nycof = nycof + 1
            if (nycof .gt. nydim ) then
               if (nydim .gt. 0) then
                 nydim = 2*nydim
!cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                 call ralloc(nya, nydim, 'NYA', 'SETCOF')
                 call ralloc(ya, nydim, 'YA', 'SETCOF')
!                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
               else
!cjb similar here
                 nydim = 64
                 call alloc(nya, nydim, 'NYA', 'SETCOF')
                 call alloc(ya, nydim, 'YA', 'SETCOF')
!                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
!                 NYDIM = NYCOF
!cjb
               endif
            endif
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
      do ib = 1, nw
         ilabel = ilabel + nwterm
         if (ib == j) cycle
         nkjib = nkj(ib)
         if (nakj*nak(ib) > 0) then
            kmin = abs((nkjj - nkjib)/2)
         else
            kmin = abs((nkjj - nkjib)/2) + 1
         endif
         kmax = (nkjj + nkjib)/2

         do k = kmin, kmax, 2

            !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
            sumr = 0.d0
            do nb = 1, nblock
               do ir = myid + 1, ncfblk(nb), nprocs
                  sumr = sumr + dsubrs(eol,ir,ir,nb)*gco(k,ir + ncfpast(nb),j,&
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
                 call ralloc (xa, nxdim, 'XA', 'SETCOF')
                 call ralloc (nxa, nxdim, 'NXA', 'SETCOF')
!                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
               else
                  nxdim = 64
                 call alloc (xa, nxdim, 'XA', 'SETCOF')
                 call alloc (nxa, nxdim, 'NXA', 'SETCOF')
!                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
!                  NXDIM = NXCOF
!cjb
               endif
            endif
            xa(nxcof) = xkab
            nxa(nxcof) = k + ilabel
            !ENDIF
         end do
      end do

!=======================================================================
!   Subroutine setham called from matrix had gone through the mcp
!   files once. setcof will do it again. (Go to setmcp to see
!   contents of these records.
!=======================================================================

      do nfile = 30, 32 + kmaxf
         rewind (nfile)
         if (nfile /= 30) then
            do i = 1, 3
               read (nfile)
               cycle
            end do
         else
            do i = 1, 3
               read (nfile)
               read (nfile)
            end do
         endif
      end do

!=======================================================================
!   Generate DA coefficients; these arise from the one-electron
!   integrals
!=======================================================================

      do jblock = 1, nblock

         !*** Read in IROW from file mcp.30 ***
         read (30) mcplab, jblockt, ncf
         if (jblockt /= jblock) then
            write (istde, *) myname, '1: jblockt .NE. jblock'
            error stop
         endif
         read (30) nelmntgg
         nelmnt = int(nelmntgg,kind=long)
         call alloc (irow, nelmnt, 'IROW', 'SETCOF')
         read (30) (iendcdum,i=myid + 1,ncf,nprocs), (irow(i),i=1,nelmnt)

         !*** Block info file mcp.31 ***
         read (31) mcplab, jblockt, ncf, ncoeff
         if (jblockt /= jblock) then
            write (istde, *) myname, '2: jblockt .NE. jblock'
            error stop
         endif

         !*** Loop over labels having non-zero coefficients
         !*** it exits when no more labels for the block

         l123: do while(.true.)
            read (31, iostat=ios) lab, ncontr

            ! 0, 0 marks the end of a block. This is the normal exit

            if (lab==0 .and. ncontr==0) then
               call dalloc (irow, 'IROW', 'SETCOF')
               exit                              ! Actually to next block
            endif

            !*** Decode the labels of I(ab) ***
            ia = mod(lab,key)
            ib = lab/key

            ! At least one orbital should be J in order to have
            ! non-zero value; otherwise, goto next label.

            if (ia/=j .and. ib/=j) then
               read (31)                         ! No contributions from this integral; skip
               cycle                             ! to next label
            endif

            if (ncontr > ndim) then
               call dalloc (coeff, 'COEFF', 'SETCOF')
               call dalloc (iclmn, 'ICLMN', 'SETCOF')
               call dalloc (indx, 'INDX', 'SETCOF')
               ndim = ncontr
               call alloc (coeff, ndim, 'COEFF', 'SETCOF' )
               call alloc (iclmn, ndim, 'ICLMN', 'SETCOF')
               call alloc (indx, ndim, 'INDX', 'SETCOF')
            endif

            ! Read the column index, the sparse matrix index, and the
            ! coefficient for all contributions from this integral

            read (31) (iclmn(i),indx(i),coeff(i),i=1,ncontr)

            ! Add up all the contributions from this integral;
            ! off-diagonal contributions have double the weight

            sum = 0.d0
            do i = 1, ncontr
               ir = irow(indx(i))
               ic = iclmn(i)

               contr = dsubrs(eol,ir,ic,jblock)*coeff(i)

               if (ir /= ic) contr = contr + contr
               sum = sum + contr
            end do
            sum = 0.5d0*sum/ucfj

            ! Put coefficients in the list. Since there is always
            ! (almost) some repetition from different blocks, a check
            ! and merge is performed. This will significantly reduce
            ! the NDCOF and thus the number of calls to YZK later.

            !*** Find the right counting parameter ***
            if (ia == j) then
               ithis = ib
            else
               ithis = ia
            endif

            !*** Check it against the previously recorded ***
            if (jblock > 1) then
               do i = 1, ndcof
                  if (nda(i) /= ithis) cycle     ! found, add the value
                  da(i) = da(i) + sum
                  cycle  l123
               end do
            endif

            !*** Not found in the record, add an item ***
            ndcof = ndcof + 1
            if (ndcof > nddim) then
               if (nddim > 0) then
                  nddim = 2*nddim
                  call ralloc (da, nddim, 'DA', 'SETCOF')
                  call ralloc (nda, nddim, 'NDA', 'SETCOF')
!                  CALL RALLOC (DA, NDCOF, 'DA', 'SETCOF')
!                  CALL RALLOC (NDA, NDCOF, 'NDA', 'SETCOF')
               else
                  nddim = 64
                 call alloc (da, nddim, 'DA', 'SETCOF')
                  call alloc (nda, nddim, 'NDA', 'SETCOF')
!                  CALL ALLOC (DA, NDCOF, 'DA', 'SETCOF')
!                  CALL ALLOC (NDA, NDCOF, 'NDA', 'SETCOF')
!                  NDDIM = NDCOF
               endif
            endif
            da(ndcof) = sum                      ! print*, DA(NDCOF), ndcof, myid, 'myid'
            nda(ndcof) = ithis

         end do l123                             ! For labels
      end do                                     ! For blocks

!=======================================================================
!   Generate YA and XA coefficients; these arise from the two-electron
!   integrals
!=======================================================================

      do nfile = 32, 32 + kmaxf

!         ...Re-position file mcp.30

         rewind (30)
         do i = 1, 6
            read (30)
         end do

!=======================================================================
!   Loop over blocks again, this time, for V-coefficients
!=======================================================================

         do jblock = 1, nblock
!           ...Read in IROW from file mcp.30
            read (30) mcplab, jblockt, ncf
            if (jblockt /= jblock) then
               write (istde, *) myname, ':3 jblockt .NE. jblock'
               error stop
            endif
            read (30) nelmntgg
            nelmnt = int(nelmntgg,kind=long)
            call alloc (irow, nelmnt, 'IROW', 'SETCOF')
            read (30) (iendcdum,i=myid + 1,ncf,nprocs), (irow(i),i=1,nelmnt)

            read (nfile) mcplab, jblockt, ncf, ncoeff
            if (jblockt /= jblock) then
               write (istde, *) myname, ':4 jblockt .NE. jblock'
               error stop
            endif

            k = nfile - 32                       ! multipolarity of the integral

!=======================================================================
!   Attempt to read another block of data
!=======================================================================

  999       continue
            read (nfile, iostat=ios) lab, ncontr
!
            if (lab==0 .and. ncontr==0) then
               call dalloc (irow, 'IROW', 'SETCOF')
               cycle
            endif
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
               read (nfile)
               go to 999
            endif

            !*** At least one subshell index matches; allocate storage
            !*** for reading in the rest of this block
            if (ncontr > ndim) then
               call dalloc (coeff, 'COEFF', 'SETCOF')
               call dalloc (iclmn, 'ICLMN', 'SETCOF')
               call dalloc (indx, 'INDX', 'SETCOF')
               ndim = ncontr
               call alloc (coeff, ndim, 'COEFF', 'SETCOF')
               call alloc (iclmn, ndim, 'ICLMN', 'SETCOF')
               call alloc (indx, ndim, 'INDX', 'SETCOF')
            endif

            !*** Read column index, sparse matrix index, and
            !*** coefficient for all contributions from this integral
            read (nfile) (iclmn(i),indx(i),coeff(i),i=1,ncontr)

            !*** Add up all the contributions from this integral;
            !*** off-diagonal contributions have double the weight
            sum = 0.d0
            do i = 1, ncontr
               ir = irow(indx(i))
               ic = iclmn(i)
               contr = dsubrs(eol,ir,ic,jblock)*coeff(i)
               if (ir /= ic) contr = contr + contr
               sum = sum + contr
            end do
            sum = 0.5d0*sum/ucfj

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

               if (ithis == (-911)) error stop 'ithis .EQ. -911'

               !*** Check ithis against the previously recorded ***
               if (jblock > 1) then
                  do i = 1, nxcof
                     if (nxa(i) /= ithis) cycle
                     xa(i) = xa(i) + sum
                     go to 999
                  end do
               endif

               !*** Not found in records, add an item ***
               nxcof = nxcof + 1
               if (nxcof > nxdim) then
                  if (nxdim > 0) then
                    nxdim = 2*nxdim
!cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                    call ralloc (xa, nxdim, 'XA', 'SETCOF')
                    call ralloc (nxa, nxdim, 'NXA', 'SETCOF')
!                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                 else
                    nxdim = 64
                    call alloc (xa, nxdim, 'XA', 'SETCOF')
                    call alloc (nxa, nxdim, 'NXA', 'SETCOF')
!                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
!                  NXDIM = NXCOF
!cjb
                  endif
               endif
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
                  endif
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
                        go to 999
                     end do
                  endif

                  !*** Not found, add an item ***
                  nycof = nycof + 1
                  if (nycof > nydim) then
                     if (nydim > 0) then
                       nydim = 2*nydim
!cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                       call ralloc(nya, nydim, 'NYA', 'SETCOF')
                       call ralloc(ya, nydim, 'YA', 'SETCOF')
!                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
                     else
                       nydim = 64
                       call alloc(nya, nydim, 'NYA', 'SETCOF')
                       call alloc(ya, nydim, 'YA', 'SETCOF')
!                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
!                 NYDIM = NYCOF
!cjb
                     endif
                  endif
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
                        go to 999
                     end do
                  endif

                  !*** Not found, add an item ***
                  nxcof = nxcof + 1
                  if (nxcof > nxdim) then
                     if (nxdim > 0) then
                       nxdim = 2*nxdim
!cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                       call ralloc (xa, nxdim, 'XA', 'SETCOF')
                       call ralloc (nxa, nxdim, 'NXA', 'SETCOF')
!                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                     else
                       nxdim = 64
                       call alloc (xa, nxdim, 'XA', 'SETCOF')
                       call alloc (nxa, nxdim, 'NXA', 'SETCOF')
!                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
!                  NXDIM = NXCOF
!cjb
                     endif
                  endif
                  xa(nxcof) = sum + sum
                  nxa(nxcof) = ithis

               endif

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

               if (ithis==(-911) .or. ithis2==(-911)) error stop 'ithis2'

               !*** Check the previously recorded for YA
               if (jblock > 1) then
                  do i = 1, nycof
                     if (nya(i) /= ithis) cycle
                     ya(i) = ya(i) + sum + sum
                     go to 456
                  end do
               endif

               ! Not found, add an item ***
               nycof = nycof + 1
               if (nycof > nydim) then
                  if (nydim > 0) then
                    nydim = 2*nydim
                    call ralloc (nya, nydim, 'NYA', 'SETCOF')
                    call ralloc (ya, nydim, 'YA', 'SETCOF')
!                     CALL RALLOC (NYA, NYCOF, 'NYA', 'SETCOF')
!                     CALL RALLOC (YA, NYCOF, 'YA', 'SETCOF')
                  else
                    nydim = 64
                     call alloc (nya, nydim, 'NYA', 'SETCOF')
                     call alloc (ya, nydim, 'YA', 'SETCOF')
!                     CALL ALLOC (NYA, NYCOF, 'NYA', 'SETCOF')
!                     CALL ALLOC (YA, NYCOF, 'YA', 'SETCOF')
!                     NYDIM = NYCOF
                  endif
               endif
               ya(nycof) = sum + sum
               nya(nycof) = ithis

  456          continue

               !*** Check the previously recorded for XA
               if (jblock > 1) then
                  do i = 1, nxcof
                     if (nxa(i) /= ithis2) cycle
                     xa(i) = xa(i) + sum
                     go to 999
                  end do
               endif

               ! Not found, add an item ***
               nxcof = nxcof + 1
               if (nxcof > nxdim) then
                  if (nxdim > 0) then
                    nxdim = 2*nxdim
!cjb RALLOC ALLOC NXA, XA, (NXDIM NXCOF)
                   call ralloc (xa, nxdim, 'XA', 'SETCOF')
                   call ralloc (nxa, nxdim, 'NXA', 'SETCOF')
!                  CALL RALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL RALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
                 else
                   call alloc (xa, nxdim, 'XA', 'SETCOF')
                   call alloc (nxa, nxdim, 'NXA', 'SETCOF')
!                  CALL ALLOC (XA, NXCOF, 'XA', 'SETCOF')
!                  CALL ALLOC (NXA, NXCOF, 'NXA', 'SETCOF')
!                  NXDIM = NXCOF
!cjb
                  endif
               endif
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
                     ya(i) = ya(i) + 4.d0*sum
                     go to 999
                  end do
               endif

               ! Not found, add an item ***
               nycof = nycof + 1
               if (nycof > nydim) then
                  if (nydim > 0) then
                    nydim = 2*nydim
!cjb RALLOC ALLOC NYA, YA, (NYDIM NYCOF)
                    call ralloc(nya, nydim, 'NYA', 'SETCOF')
                    call ralloc(ya, nydim, 'YA', 'SETCOF')
!                 CALL RALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL RALLOC(YA, NYCOF, 'YA', 'SETCOF')
                 else
!cjb similar here
                   nydim = 64
                   call alloc(nya, nydim, 'NYA', 'SETCOF')
                   call alloc(ya, nydim, 'YA', 'SETCOF')
!                 CALL ALLOC(NYA, NYCOF, 'NYA', 'SETCOF')
!                 CALL ALLOC(YA, NYCOF, 'YA', 'SETCOF')
!                 NYDIM = NYCOF
!cjb
                  endif
               endif
               ya(nycof) = 4.d0*sum
               nya(nycof) = ithis

            end select

            go to 999
         end do                                  ! loop for V-Coefficients
      end do

!=======================================================================
!   Deallocate storage for arrays local to this routine
!=======================================================================
      call dalloc (coeff, 'COEFF', 'SETCOF')
      call dalloc (iclmn, 'ICLMN', 'SETCOF')
      call dalloc (indx, 'INDX', 'SETCOF')

      return
      end subroutine setcof
