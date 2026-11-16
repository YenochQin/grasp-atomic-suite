!***********************************************************************
!                                                                      *
      subroutine maneig(iatjpo, iaspar, nelmnt_a)
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
!
!cjb comments [re-introduced] from GRASP-2013
!   iatjpo - Output. (2j+1) of the dominant coefficients
!   iaspar - Output, Parity (1 or -1) of the dominant coefficients
!   nelmnt - Number of non-zero matrix elements on the current process *
!            stored in the common block                                *
!   nelmnt_a - The maximum number of non-zero matrix elements for all  *
!             nodes                                                    *
!                                                                      *
!     value  meaning  of  method used IV                               *
!       1    LAPACK                  =>  ncf < IOLPCK                  *
!       2    Davidson, Dense-Memory  =>  dense requires less memory    *
!       3    Davidson, Sparse-Memory =>  sparse requires less memory   *
!       4    Davidson, Sparse-Disk   =>  memory requirement too large  *
!                                        for at least one  process     *
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
      use mpi_c
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
      use spicmvmpi_i
      use suite_linear_algebra, only: dcopy, dspevx, dswap
      use iniestmpi_i, only: iniestmpi
      use mpiu_i, only: gdsummpi
      use mpiu_i, only: stopmpi
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
      integer(long)        :: nelmnt_a
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter      :: iolpck = 2000
! GG      REAL(DOUBLE), PARAMETER :: ABSTOL = 1.0D-10
!cjb  NINCOR
!cjb  INTEGER, PARAMETER      :: NINCOR = 1         ! To enforce DISK
      integer, parameter      :: nincor = 268435456 ! = 2 GB or  memory
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!cff   ... variable for deciding between sparse or dense
!     Memory need:
!     sparse : nstore_s = nelmnt + (isize*(nelmnt + ncf+1))/8
!     dense  : nstore_d = (ncf*ncf+1))/(2*nprocs)    (on average)
!     Use Davidson (disk)  if nstore_s > NINCOR
!     Use Davidson (memory, sparse) if nstore_s < nstore_d
!     Use Davidson (memory, dense) if nstore_d < nstore_s
!-----------------------------------------------
      integer(long) :: nstore_s, nstore_d
      integer :: nrows, i, ndense, ncfdum, iccutdum, myiddum, &
         nprocsdum, iofset, nelc, ir, nvecmn, nvex, m, info, loc, nbrkev,   &
         imv, ndense_l, lim, lwork, liwork, maxitr, mblock, nend, &
         ilow, ihigh, niv, nloops, nmv, j, ia,   &
         idummy
      real(double) :: elsto,  dummy, &
          diatmp, crite, critc, critr, ortho, dmungo, amax, wa, abstol
! GG          DIATMP, CRITE, CRITC, CRITR, ORTHO, DMUNGO, AMAX, WA
      logical :: hiend, ldisc, sparse
      character(len=8) :: cnum
      real(double), dimension(:), pointer :: w, z, work, diag
      integer, dimension(:), pointer :: iwork, ifail, jwork
!-----------------------------------------------
!
!-----------------------------------------------------------------------
      abstol = 2*dlamch('S')
      if (myid == 0) write (6, *) 'Calling maneig...'

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
         if (myid .eq. 0) write (24, *) 'Trivial eigenvalue problem.'

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
         if(myid==0) write(*,*) "ncf=",ncf, " iolpck=", iolpck
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

            ndense = (ncf*(ncf + 1))/2
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

! Let each node have a complete copy of EMT

            call gdsummpi (emt, ndense)

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
            if (info /= 0)                                             &
                call stopmpi('maneig: Failure in DSPEVX [LAPACK]',myid)
            call dalloc (work, 'WORK', 'MANEIG')
            call dalloc (iwork, 'IWORK', 'MANEIG')
            call dalloc (ifail, 'IFAIL', 'MANEIG')
            call dalloc (emt, 'EMT', 'MANEIG')

!   Store the eigenpairs in their proper positions EVAL() and EVEC()

            call alloc (eval, nvec,'EVAL', 'MANEIG')
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
            if (myid == 0) write (24,*)'DVDSON routine selected' &
                                      //' for eigenvalue problem;'
!--------------------------------------------------------------
            if (myid == 0) print *,'zou: nelmnt',nelmnt
!            print *, 'myid, nprocs,nincor', myid, nprocs, nincor

!           If on any process more memory is needed than NINCOR, use disc.
!           If not, should the matrix be sparse or dense (distributed)?
!           sparse requires = nelmnt + (isize*(nelmnt+ncf+1))/8
!           dense  storage = (ncf*(ncf+1)/2/nproc - on average
!           (because of integer arithmetic nproc-1 is added to numerator

!           This is a comparion of matrix elements only.
!
            nstore_s = nelmnt_a
            nstore_d = ((ncf+nprocs-1)/nprocs)
            nstore_d = nstore_d*((ncf+1)/2)
            if (myid .eq. 0) write(*,*) "nstore_s=", nstore_s,   &
                             "nstore_d=", nstore_d, 'nincor=', nincor
            ldisc = .false.
            if ( nstore_s > nincor) then
!              .. we need to use disk
               ldisc =.true.
               sparse=.true.
            else

               if (nstore_s .lt. nstore_d ) then
                  sparse = .true.
               else
                  sparse = .false.
               endif
            endif

            call alloc (diag,ncf,'DIAG', 'MANEIG')
            call dinit (ncf, 0.d0, diag, 1)
            if (ldisc) then
!-----------------------------------------------------------------------
!
! (2.2.1) - DVDSON --- Disk, load diagonal
!
!-------------------------------------------------------
               if (myid == 0)  then
                  write (*, *) ' matrix stored on disc;'
                  write (24, *) ' matrix stored on disc;'
               end if

!   Disk storage; necessarily sparse; one column of the matrix in
!   memory

               imv = 1

!   Load diagonal - Each node will have the same, complete copy
!   after this if block
               print *, myid, 'Reading diagonal, unit imcdf =', imcdf
               read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
               if (ncf/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) stop &
                  'maneig:2'

               do i = myid + 1, ncf, nprocs
                  read (imcdf) nelc, elsto, (dummy,ir=2,nelc), diatmp, &
                                            (idummy, ir=1,nelc)
                  diag(i) = diatmp - eav
               end do
               print *, ' finished'
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
                  if (myid == 0) then
                     write (24, *) ' sparse matrix stored in memory'
                     write ( *, *) ' sparse matrix stored in memory'
                  end if

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
                  end do
               else
!-----------------------------------------------------------------------
!
! (2.2.2.2) - DVDSON --- Memory, load all, dense
!
!-------------------------------------------------------
                  if (myid .eq. 0) then
                     write (24,*) ' full matrix stored in in memory'
                     write ( *,*) ' full matrix stored in in memory'
                  end if

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

! Make diagonals global, no matter it is disk or memory mode

            call gdsummpi (diag, ncf)
!
!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECMX
!   plus a smaller number if storage is severely constrained
!
            lim = min(ncf,2*nvecmx + 60)
!           lwork = 2*ncf*lim + lim*lim + (nvecmx+10)*lim + nvecmx
            lwork = 2*ncf*lim + lim*lim*2 + 11*lim + nvecmx
            call alloc (work, lwork, 'WORK', 'MANEIG')
            work(1:lwork) = 0.0d0
            liwork = 6*lim + nvecmx
            call alloc (iwork, liwork, 'IWORK', 'MANEIG')
!*changed by Misha 02/12/97
            crite = 1.0d-17
!GG            CRITC = 1.0D-08
!GG            CRITR = 1.0D-08
!GG            ORTHO = MAX(1D-8,CRITR)
            critc = 1.0d-09
            critr = 1.0d-09
            ortho = max(1d-9,critr)
! end of changes

!            maxitr = MAX (nvecmx*100, ncf/10)
            maxitr = max(nvecmx*100,ncf/10)
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
              if (myid .eq. 0) print *, ' Sparse - Disk, iniestsd'
               call posfile (0, imcdf, nposition)! was within iniestsd before
               call iniestsd (iolpck, ncf, myid, nprocs, niv, work, imcdf, eav)
              print *, 'Returned from iniestsd '
              if (ncf.gt. iolpck) then
                print *, 'Calling GDVD'
                call gdvd (spodmv,ncf,lim,diag,ilow,ihigh,            &
                  jwork,niv,mblock,crite,critc, critr,ortho,maxitr,   &
                  work,lwork,iwork,liwork,hiend,nloops,               &
                  nmv,ierr)
              end if

            case (2)
!******************** sparse and matrix in memory ********************
               if (myid .eq. 0) print *, ' Sparse - Memory, iniestmpi'
               call iniestmpi (iolpck, ncf,niv,work,emt,iendc,irow)
               if(ncf.gt. iolpck) then
                 call gdvd (spicmvmpi,ncf,lim,diag,ilow,ihigh,        &
                  jwork,niv,mblock,crite,critc, critr,ortho,maxitr,   &
                  work,lwork,iwork,liwork,hiend,nloops,               &
                  nmv,ierr)
               end if

               call dalloc (emt, 'EMT', 'MANEIG')
               call dalloc (irow, 'IROW', 'MANEIG')
               call dalloc (iendc, 'IENDC', 'MANEIG')

            case (3)
!*************************** dense and in memory **********************
               if (myid .eq. 0) print *, ' Dense - Memory, iniestdm'
              call iniestdm (iolpck,ncf,niv,work,emt)
              if (ncf.gt. iolpck) then
                call gdvd (dnicmv,ncf,lim,diag,ilow,ihigh,            &
                    jwork,niv,mblock,crite,critc, critr,ortho,maxitr, &
                    work,lwork,iwork,liwork,hiend,nloops,             &
                    nmv,ierr)
              end if
               call dalloc (emt, 'EMT', 'MANEIG')
            end select
!************************************************************************
            call dalloc (diag, 'DIAG', 'MANEIG')
            call dalloc (iwork, 'IWORK', 'MANEIG')
            call dalloc (jwork, 'JWORK', 'MANEIG')
            if (myid .eq. 0) then
               write (24, *) ' ', nloops, ' iterations;'
               write (24, *) ' ', nmv, ' matrix-vector multiplies.'
            endif
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
