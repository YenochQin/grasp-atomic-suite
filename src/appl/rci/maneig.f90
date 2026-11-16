!***********************************************************************
!                                                                      *
      subroutine maneig(iatjpo, iaspar)
!                                                                      *
!   This module  manages the  operation of the  eigensolvers and the   *
!   storage of the eigenpairs.  There are two principal branches:      *
!                                                                      *
!      (1) Matrix of order 1: the trivial case                         *
!      (2) Matrix of order greater than 1: there are two branches      *
!             (i) Matrices of order less than or equal to IOLPCK:      *
!                 eigenpairs are found using LAPACK SUBROUTINEs        *
!            (ii) Matrices of order greater than IOLPCK: eigenpairs    *
!                 are found using DVDSON; this involves up to three    *
!                 steps:                                               *
!                    (a) The matrix is analysed to determine its       *
!                        block structure (only irreducibe matrices     *
!                        are correctly treated by DVDSON)              *
!                    (b) Eigenpairs are extracted for each block       *
!                    (c) The appropriate eigenpairs are selected and   *
!                        stored                                        *
!                 Different methods of storage and different           *
!                 versions of the matrix-vector multiply are used      *
!                 depending upon the order and density of the matrix   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, posfile,         *
!                        RALLOC.                                       *
!               [RCI92]: DNICMV, SPICMV2, SPODMV.                      *
!               [DVDSON]: DVDSON.                                      *
!               [AUXBLAS]: DINIT/SINIT.                                *
!               [BLAS]: DCOPY/SCOPY, DSWAP/SSWAP.                      *
!               [LAPACK]: DSPEVX/SSPEVX.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
!   Modified Misha Saparov                                  Feb 1997   *
!            Charlotte F. Fischer                           May 1997   *
!                 Except for the disk version, all matrices have       *
!                 diagonals  shifted by EAV                            *
!  All arrays allocated here are de-allocated except pneval and pnevec *
!  which will be de-allocated in matrix.                               *
!  Block Version  By Xinghong He          Last revision: 18 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:38   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double, long
      use memory_man
      use eigv_c
      use fposition_c
      use hmat_c
      use orb_c,           only: ncf, nw, iqa
      use prnt_c
      use where_c
      use wchblk_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dnicmv_i
      use spicmv2_i
      use spodmv_i
      use posfile_i
      use dinit_i
      use dspevx_i
      use iniestsd_i
      use gdvd_i
      use iniest2_i
      use iniestdm_i
      use itjpo_i
      use ispar_i
      use suite_linear_algebra, only: dcopy, dspevx, dswap
      use suite_linear_algebra, only: dlamch
      implicit none
!-----------------------------------------------
!  E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: iatjpo
      integer, intent(out) :: iaspar
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!cjb  INTEGER, PARAMETER      :: IOLPCK = 1000
      integer, parameter      :: iolpck = 2000
! GG      REAL(DOUBLE), PARAMETER :: ABSTOL = 1.0D-10
!cjb  NINCOR
!cjb  INTEGER, PARAMETER      :: NINCOR = 1         ! To enforce DISK
      integer, parameter      :: nincor = 268435456 ! = 2 GB or  memory
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer(long) :: nstore
      integer :: myid, nprocs, nrows, i, ndense, ncfdum, iccutdum, myiddum, &
         nprocsdum, iofset, nelc, ir, nvecmn, nvex, m, info, loc, nbrkev,   &
         imv, ndense_l, lim, lwork, liwork, maxitr, mblock, nend, &
         ilow, ihigh, niv, nloops, nmv, ierr, j, ia
      real(double) :: elsto,  dummy, &
          diatmp, crite, critc, critr, ortho, dmungo, amax, wa, abstol
! GG          DIATMP, CRITE, CRITC, CRITR, ORTHO, DMUNGO, AMAX, WA
      logical :: hiend, ldisc, sparse
      character(len=8) :: cnum
      real(double), dimension(:), pointer :: w, z, work, diag
      integer, dimension(:), pointer :: iwork, ifail, jwork
!-----------------------------------------------------------------------
      abstol = 2*dlamch('S')
      myid = 0
      nprocs = 1
      !IF (MYID == 0) WRITE (6, *) 'Calling maneig...'

!  (nrows+1) is the number of records of the present block's .res file

      nrows = (ncf - myid - 1 + nprocs)/nprocs
      if (ncf < nprocs) nrows = ncf/(myid + 1)
      !CALL posfile (1, imcdf, nrows+1)
      call posfile (0, imcdf, nposition)

      if (ncf == 1) then
!-----------------------------------------------------------------------
!
! (1) - Trivial    ncf = 1
!
!-------------------------------------------------------
         write (24, *) 'Trivial eigenvalue problem.'

!   Matrix of order 1: the trivial case; we assume that the value
!   of EAV is available

         call alloc (eval, 1,'EVAL', 'MANEIG' )
         call alloc (evec, 1, 'EVECO', 'MANEIG')
         eval(1) = 0.d0
         evec(1) = 1.d0

         ! Still read through the .res file
!GG
!GG Gediminas NIST 2005.11.03
!GG         READ (imcdf)
         do i = 1, nrows + 1
            read (imcdf)
         end do

      else                                       !if-2
!-----------------------------------------------------------------------
!
! (2) - Non trivial
!
!-------------------------------------------------------
!
!   Matrix of order greater than 1; how many elements in a triangle?
!
         ndense = (ncf*(ncf + 1))/2

         if (ncf <= iolpck) then
!-----------------------------------------------------------------------
!
! (2.1) - LAPACK    Dense, Memory,
!
!-------------------------------------------------------
            if (myid == 0) then
               write (6, *) &
                  'LAPACK routine DSPEVX selected for eigenvalue problem.'
               write (24, *) &
                  'LAPACK routine DSPEVX selected for eigenvalue problem.'
            endif

!   Allocate storage for the dense representation of the matrix
!   and initialize emt

            call alloc (emt, ndense, 'EMT', 'MANEIG')
            call dinit (ndense, 0.0d00, emt, 1)

!   Read the matrix into position from the disc file; it's already
!   been properly positioned.

            call alloc (work, ncf,'WORK', 'MANEIG' )
            call alloc (irow, ncf, 'IROW', 'MANEIG')
            read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
            if (ncf/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) stop &
               'maneig:1'

            do i = myid + 1, ncf, nprocs
               iofset = (i*(i - 1))/2
               read (imcdf) nelc, elsto, (work(ir),ir=1,nelc), (irow(ir),ir=1,&
                  nelc)
               ! In the row-mode of the lower triangle,
               ! diagonal is the last one
               do ir = 1, nelc - 1
                  emt(iofset+irow(ir)) = work(ir)
               end do
               emt(iofset+irow(nelc)) = work(nelc) - eav

            end do

            call dalloc (work, 'WORK', 'MANEIG')
            call dalloc (irow, 'IROW', 'MANEIG')

!   Find the eigenpairs
!
!    ivec() - serial numbers of eigenstates of the current block
!    iccmin() - serial numbers of eigenstates of all blocks.
!    nvecmn - minimum serial number of the eigenstates of the block
!    nvecmx - maximum .............
!    nvex - clear from def: NVECMX-NVECMN+1

            nvecmn = ncf
            do i = 1, nvec
               nvecmn = min(nvecmn,ivec(i))
            end do
            nvex = nvecmx - nvecmn + 1
            call alloc (w, nvex, 'W', 'MANEIG')
            call alloc (z, ncf*nvex,'Z', 'MANEIG' )
            call alloc (work, ncf*8,'WORK', 'MANEIG' )
            call alloc (iwork, ncf*5,'IWORK', 'MANEIG' )
! GG            CALL ALLOC (IFAIL, NVEX, 'IFAIL', 'MANEIG')
            call alloc (ifail, ncf, 'IFAIL', 'MANEIG')
            call dspevx ('V', 'I', 'U', ncf, emt, dummy, dummy, nvecmn, nvecmx&
               , abstol, m, w, z, ncf, work, iwork, ifail, info)
            if (info /= 0) stop 'maneig: Failure in DSPEVX [LAPACK]'
            call dalloc (work, 'WORK', 'MANEIG')
            call dalloc (iwork, 'IWORK', 'MANEIG')
            call dalloc (ifail, 'IFAIL', 'MANEIG')
            call dalloc (emt, 'EMT', 'MANEIG')

!   Store the eigenpairs in their proper positions EVAL() and EVEC()

            call alloc (eval, nvec,'EVAL', 'MANEIG' )
            call alloc (evec, ncf*nvec, 'EVEC', 'MANEIG')

            do i = 1, nvec
               loc = ivec(i)
               eval(i) = w(loc - nvecmn + 1)
               iofset = ncf*(i - 1)
               loc = ncf*(loc - nvecmn)
               call dcopy (ncf, z(loc + 1:loc+(ncf)), 1, evec(iofset+1:iofset+(ncf)), 1)
            end do
            call dalloc (w, 'W', 'MANEIG')
            call dalloc (z, 'Z', 'MANEIG')

         else
!-----------------------------------------------------------------------

! (2.2) - DVDSON --- preparation work
!
!-------------------------------------------------------
           write (24,*)'DVDSON routine selected for eigenvalue problem;'

!   Sparse or dense matrix multiply? On disc or in core?

!GG            NBRKEV = (NCF + 1)*(NCF + 1)/3       ! Normal
            !NBRKEV =  1                      ! To enforce DENSE
            !NBRKEV = (NCF*(NCF+1)) / 2 + 1   ! To enforde SPARSE
!--------------------------------------------------------------
! Uncomment out the following to force sparse storage
!
!
!--------------------------------------------------------------
               sparse = .true.
               nstore = nelmnt + nelmnt/2 + (ncf + 1)/2

            call alloc (diag, ncf, 'DIAG', 'MANEIG')

            if (nstore > nincor) then
!-----------------------------------------------------------------------
!
! (2.2.1) - DVDSON --- Disk, load diagonal
!
!-------------------------------------------------------
               write (24, *) ' matrix stored on disc;'

!   Disk storage; necessarily sparse; one column of the matrix in
!   memory

               ldisc = .true.
               sparse = .true.
               imv = 1

!   Load diagonal - Each node will have the same, complete copy
!   after this if block

               read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
               if (ncf/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) stop &
                  'maneig:2'

               do i = myid + 1, ncf, nprocs
                  read (imcdf) nelc, elsto, (dummy,ir=2,nelc), diatmp
                  diag(i) = diatmp - eav
               end do

            else
!-----------------------------------------------------------------------
!
! (2.2.2) - DVDSON --- Memory, load all
!
!-------------------------------------------------------
!
!   Core storage; load matrix into memory
!
               ldisc = .false.
               if (sparse) then
!-----------------------------------------------------------------------

! (2.2.2.1) - DVDSON --- Memory, load all, sparse
!
!-------------------------------------------------------
                  if (myid == 0) write (24, *) &
                     ' matrix stored in sparse representation in core;'

                  imv = 2
                  write (6, *) 'nelmnt = ', nelmnt
                  call alloc (emt, nelmnt, 'EMT', 'MANEIG')
                  call alloc (irow, nelmnt, 'IROW', 'MANEIG')
                  call alloc (iendc, 0, ncf , 'IENDC', 'MANEIG')
                  iofset = 0
                  iendc(0) = 0

                  read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
                  if (ncf/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) &
                     stop 'maneig:3'
                  do i = myid + 1, ncf, nprocs
                     read (imcdf) nelc, elsto, (emt(ir+iofset),ir=1,nelc), (&
                        irow(ir + iofset),ir=1,nelc)
                     emt(nelc+iofset) = emt(nelc+iofset) - eav
                     diag(i) = emt(nelc+iofset)
                     iofset = iofset + nelc
                     iendc(i) = iofset
!                     WRITE (31 + MYID, *) I, IENDC(I), DIAG(I)

                  end do
               else
!-----------------------------------------------------------------------
!
! (2.2.2.2) - DVDSON --- Memory, load all, dense
!
!-------------------------------------------------------
                  write (24, *) &
                     ' matrix stored in full representation in core;'

                  imv = 3

! Find NDENSE_L, the number of elements on the node (dense form)

                  ndense_l = 0
                  do i = myid + 1, ncf, nprocs
                     ndense_l = ndense_l + i
                  end do

                  call alloc (emt, ndense_l, 'EMT', 'MANEIG')
                  call dinit (ndense_l, 0.0d00, emt, 1)
                  call alloc (work, ncf, 'WORK', 'MANEIG')
                  call alloc (irow, ncf, 'IROW', 'MANEIG')

                  read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
                  if (ncf/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) &
                     stop 'maneig:4'

                  iofset = 0
                  do i = myid + 1, ncf, nprocs
                     read (imcdf) nelc, elsto, (work(ir),ir=1,nelc), (irow(ir),&
                        ir=1,nelc)
                     work(nelc) = work(nelc) - eav
                     diag(i) = work(nelc)
                     do ir = 1, nelc
                        emt(iofset+irow(ir)) = work(ir)
                     end do
                     iofset = iofset + i
                  end do
                  call dalloc (work, 'WORK', 'MANEIG')
                  call dalloc (irow, 'IROW', 'MANEIG')

               endif
!               ...Memory mode - sparse or dense
!-----------------------------------------------------------------------
!  (2.2.2.3e)          *** E n d   o f   D V D S O N   m e m o r y
!-----------------------------------------------------------------------
            endif
!            ...Disk or Memory
!-----------------------------------------------------------------------
!  (2.2.3e)      *** E n d   o f   D V D S O N
!-----------------------------------------------------------------------
!
!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECMX
!   plus a smaller number if storage is severely constrained
!
            lim = min(ncf,2*nvecmx + 60)
!           lwork = 2*ncf*lim + lim*lim + (nvecmx+10)*lim + nvecmx
            lwork = 2*ncf*lim + lim*lim*2 + 11*lim + nvecmx
            call alloc (work, lwork, 'WORK', 'MANEIG')
            liwork = 6*lim + nvecmx
            call alloc (iwork, liwork, 'IWORK', 'MANEIG')
!*changed by Misha 02/12/97
            crite = 1.0d-17
            critc = 1.0d-09
            critr = 1.0d-09
            ortho = max(1d-9,critr)
! end of changes

!            maxitr = MAX (nvecmx*100, ncf/10)
            maxitr = max(nvecmx*200,ncf/10)
            !maxitr = MIN (nvect*100, ncf)  ! FROM RSCFVU !!!
            call alloc (jwork, lim,'JWORK', 'MANEIG' )

            call alloc (eval, nvecmx, 'EVAL', 'MANEIG')
            call alloc (evec, ncf*nvecmx, 'EVEC', 'MANEIG')

            dmungo = 10.d99
            call dinit (nvecmx, dmungo, eval, 1)

!   Compute the eigenpairs in each block

            nvex = nvecmx
            if (ldisc) then
               mblock = nvex
            else
               mblock = 1
            endif
            nend = ncf*nvex

            ilow = 1
            ihigh = nvex
            niv = nvex
!************************************************************************
!
!   Call Davidson eigensolver
!
            select case (imv)
            case (1)
!******************** sparse and matrix on disk **********************
               write (6, *) ' Sparse - Disk,  iniestsd'
               call posfile (0, imcdf, nposition)! was within iniestsd before
               call iniestsd (2000, ncf, myid, nprocs, niv, work, imcdf, eav)

               !NIV = 0   ! Why equal 0 ???
               !WRITE (6, *) ' Calling gdvd(spodmv,...'
               call gdvd (spodmv, ncf, lim, diag, ilow, ihigh, jwork, niv, &
                  mblock, crite, critc, critr, ortho, maxitr, work, lwork, &
                  iwork, liwork, hiend, nloops, nmv, ierr)

            case (2)
!******************** sparse and matrix in memory ********************
               write (6, *) ' Sparse - Memory,  iniest2'
!            CALL INIEST2 (1000, NCF,NIV,WORK,EMT,IENDC,IROW)
               call iniest2 (2000, ncf, niv, work, emt, iendc, irow)
               !WRITE (*, *) NCF, NIV, (WORK(I),I=NCF*NIV + 1,NCF*NIV + NIV)
               !WRITE (*, *) LIM, ILOW, IHIGH, MBLOCK, MAXITR, LWORK, LIWORK
               !WRITE (*, *) IERR
               call gdvd (spicmv2, ncf, lim, diag, ilow, ihigh, jwork, niv, &
                  mblock, crite, critc, critr, ortho, maxitr, work, lwork, &
                  iwork, liwork, hiend, nloops, nmv, ierr)
               !WRITE (*, *) 'after gdvd...'
               !WRITE (*, *) NCF, NIV, (WORK(I),I=NCF*NIV + 1,NCF*NIV + NIV)
               !WRITE (*, *) LIM, ILOW, IHIGH, MBLOCK, MAXITR, LWORK, LIWORK
               !WRITE (*, *) HIEND, NLOOPS, NMV, IERR

               call dalloc (emt, 'EMT', 'MANEIG')
               call dalloc (irow, 'IROW', 'MANEIG')
               call dalloc (iendc, 'IENDC', 'MANEIG')

            case (3)
!*************************** dense and in memory **********************
               write (6, *) ' Dense - Memory,  iniestdm'
!            CALL INIESTDM (1000,NCF,NIV,WORK,EMT)
               call iniestdm (2000, ncf, niv, work, emt)
               call gdvd (dnicmv, ncf, lim, diag, ilow, ihigh, jwork, niv, &
                  mblock, crite, critc, critr, ortho, maxitr, work, lwork, &
                  iwork, liwork, hiend, nloops, nmv, ierr)
               call dalloc (emt, 'EMT', 'MANEIG')
            end select
!************************************************************************
            call dalloc (diag, 'DIAG', 'MANEIG')
            call dalloc (iwork, 'IWORK', 'MANEIG')
            call dalloc (jwork, 'JWORK', 'MANEIG')

            write (24, *) ' ', nloops, ' iterations;'
            write (24, *) ' ', nmv, ' matrix-vector multiplies.'

            if (ierr /= 0) then
               write (istde, *) 'MANEIG: Returned from DVDSON with'
               write (istde, *) ' IERR = ', ierr, '.'
               stop 'maneig: DVDSON wrong'
            endif

!   Put the eigenpairs in order, overwriting as necessary

            call dcopy (nvex, work(nend+1:nend+(nvex)), 1, eval, 1)
            call dcopy (ncf*nvex, work(1:ncf*nvex), 1, evec, 1)
            call dalloc (work, 'WORK', 'MANEIG')

!   Rearrange and reallocate storage for the eigenpairs
!   as necessary

            if (nvec < nvecmx) then
               call alloc (iwork, nvecmx, 'IWORK', 'MANEIG')
               do i = 1, nvecmx
                  iwork(i) = i
               end do
               do i = 1, nvec
                  iofset = ivec(i)
                  loc = iwork(i)
                  if (iofset == loc) cycle
                  call dswap (1, eval(iofset:(iofset)+((1)-1)), 1, eval(i:(i)+((1)-1)), 1)
                  iwork(i) = iwork(iofset)
                  iwork(iofset) = loc
                  iofset = ncf*(iofset - 1)
                  loc = ncf*(i - 1)
                  call dswap (ncf, evec(iofset+1:iofset+(ncf)), 1, evec(loc+1:loc+(ncf)), 1)
               end do
               call dalloc (iwork, 'IWORK', 'MANEIG')
               call ralloc (eval,  nvec, 'EVAL', 'MANEIG')
               call ralloc (evec, ncf*nvec, 'EVEC', 'MANEIG' )

            endif

         endif
! (2.3e)              *** E N D   O F    N O N - T R I V I A L   C A S E

      endif
! (3e)               *** E N D   O F    A L L

!--------------------------------------------------------------------
! Only the following quantities are needed after this routine is
! finished:
!   eval(), evec(), iatjpo, iaspar
!--------------------------------------------------------------------
!
!   Clean up eigenvectors; determine their J/P values
!
      do j = 1, nvec

!         Find the dominant component of each eigenvector

         iofset = (j - 1)*ncf

         amax = 0.d0
         do i = 1, ncf
            wa = abs(evec(i+iofset))
            if (wa <= amax) cycle
            amax = wa
            ia = i
         end do

!          Find the angular momentum and parity of the dominant component

         iatjpo = itjpo(ia)
         iaspar = ispar(ia)

!          Change sign of eigenvactor if dominant component is negative

         if (evec(ia+iofset) >= 0.d0) cycle
         evec(1+iofset:ncf+iofset) = -evec(1+iofset:ncf+iofset)
      end do

      return
      end subroutine maneig
