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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use memory_man
      use eigv_c
      use fposition_c
      use hmat_c
      use orb_c,           only: csf_count, orbital_count, csf_subshell_occupations
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
      implicit none (type, external)
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
      integer(int64) :: nstore
      integer :: myid, nprocs, nrows, i, ndense, ncfdum, iccutdum, myiddum, &
         nprocsdum, iofset, nelc, ir, nvecmn, nvex, m, info, loc, nbrkev,   &
         imv, ndense_l, lim, lwork, liwork, maxitr, mblock, nend, &
         ilow, ihigh, niv, nloops, nmv, ierr, j, ia
      real(real64) :: elsto,  dummy, &
          diatmp, crite, critc, critr, ortho, dmungo, amax, wa, abstol
! GG          DIATMP, CRITE, CRITC, CRITR, ORTHO, DMUNGO, AMAX, WA
      logical :: hiend, ldisc, sparse
      character(len=8) :: cnum
      real(real64), pointer :: w(:), z(:), work(:), diag(:)
      integer, pointer :: iwork(:), ifail(:), jwork(:)
!-----------------------------------------------------------------------
      abstol = 2*dlamch('S')
      myid = 0
      nprocs = 1
      !IF (MYID == 0) WRITE (6, *) 'Calling maneig...'

!  (nrows+1) is the number of records of the present block's .res file

      nrows = (csf_count - myid - 1 + nprocs)/nprocs
      if (csf_count < nprocs) nrows = csf_count/(myid + 1)
      !CALL posfile (1, imcdf, nrows+1)
      call posfile (0, imcdf, nposition)

      if (csf_count == 1) then
!-----------------------------------------------------------------------
!
! (1) - Trivial    ncf = 1
!
!-------------------------------------------------------
         write (24, *) 'Trivial eigenvalue problem.'

!   Matrix of order 1: the trivial case; we assume that the value
!   of EAV is available

         call alloc (shifted_eigenvalues, 1,'EVAL', 'MANEIG' )
         call alloc (csf_mixing_coefficients, 1, 'EVECO', 'MANEIG')
         shifted_eigenvalues(1) = 0._real64
         csf_mixing_coefficients(1) = 1._real64

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
         ndense = (csf_count*(csf_count + 1))/2

         if (csf_count <= iolpck) then
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

            call alloc (hamiltonian_values, ndense, 'EMT', 'MANEIG')
            call dinit (ndense, 0.0_real64, hamiltonian_values, 1)

!   Read the matrix into position from the disc file; it's already
!   been properly positioned.

            call alloc (work, csf_count,'WORK', 'MANEIG' )
            call alloc (matrix_partner_indices, csf_count, 'IROW', 'MANEIG')
            read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
            if (csf_count/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) stop &
               'maneig:1'

            do i = myid + 1, csf_count, nprocs
               iofset = (i*(i - 1))/2
               read (imcdf) nelc, elsto, (work(ir),ir=1,nelc), (matrix_partner_indices(ir),ir=1,&
                  nelc)
               ! In the row-mode of the lower triangle,
               ! diagonal is the last one
               do ir = 1, nelc - 1
                  hamiltonian_values(iofset+matrix_partner_indices(ir)) = work(ir)
               end do
               hamiltonian_values(iofset+matrix_partner_indices(nelc)) = work(nelc) - hamiltonian_energy_shift

            end do

            call dalloc (work, 'WORK', 'MANEIG')
            call dalloc (matrix_partner_indices, 'IROW', 'MANEIG')

!   Find the eigenpairs
!
!    ivec() - serial numbers of eigenstates of the current block
!    iccmin() - serial numbers of eigenstates of all blocks.
!    nvecmn - minimum serial number of the eigenstates of the block
!    nvecmx - maximum .............
!    nvex - clear from def: NVECMX-NVECMN+1

            nvecmn = csf_count
            do i = 1, nvec
               nvecmn = min(nvecmn,ivec(i))
            end do
            nvex = nvecmx - nvecmn + 1
            call alloc (w, nvex, 'W', 'MANEIG')
            call alloc (z, csf_count*nvex,'Z', 'MANEIG' )
            call alloc (work, csf_count*8,'WORK', 'MANEIG' )
            call alloc (iwork, csf_count*5,'IWORK', 'MANEIG' )
! GG            CALL ALLOC (IFAIL, NVEX, 'IFAIL', 'MANEIG')
            call alloc (ifail, csf_count, 'IFAIL', 'MANEIG')
            call dspevx ('V', 'I', 'U', csf_count, hamiltonian_values, dummy, dummy, nvecmn, nvecmx&
               , abstol, m, w, z, csf_count, work, iwork, ifail, info)
            if (info /= 0) stop 'maneig: Failure in DSPEVX [LAPACK]'
            call dalloc (work, 'WORK', 'MANEIG')
            call dalloc (iwork, 'IWORK', 'MANEIG')
            call dalloc (ifail, 'IFAIL', 'MANEIG')
            call dalloc (hamiltonian_values, 'EMT', 'MANEIG')

!   Store the eigenpairs in their proper positions EVAL() and EVEC()

            call alloc (shifted_eigenvalues, nvec,'EVAL', 'MANEIG' )
            call alloc (csf_mixing_coefficients, csf_count*nvec, 'EVEC', 'MANEIG')

            do i = 1, nvec
               loc = ivec(i)
               shifted_eigenvalues(i) = w(loc - nvecmn + 1)
               iofset = csf_count*(i - 1)
               loc = csf_count*(loc - nvecmn)
               call dcopy (csf_count, z(loc + 1:loc+(csf_count)), 1, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)
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
               nstore = stored_matrix_element_count + stored_matrix_element_count/2 + (csf_count + 1)/2

            call alloc (diag, csf_count, 'DIAG', 'MANEIG')

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
               if (csf_count/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) stop &
                  'maneig:2'

               do i = myid + 1, csf_count, nprocs
                  read (imcdf) nelc, elsto, (dummy,ir=2,nelc), diatmp
                  diag(i) = diatmp - hamiltonian_energy_shift
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
                  write (6, *) 'nelmnt = ', stored_matrix_element_count
                  call alloc (hamiltonian_values, stored_matrix_element_count, 'EMT', 'MANEIG')
                  call alloc (matrix_partner_indices, stored_matrix_element_count, 'IROW', 'MANEIG')
                  call alloc (stored_segment_end_offsets, 0, csf_count , 'IENDC', 'MANEIG')
                  iofset = 0
                  stored_segment_end_offsets(0) = 0

                  read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
                  if (csf_count/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) &
                     stop 'maneig:3'
                  do i = myid + 1, csf_count, nprocs
                     read (imcdf) nelc, elsto, (hamiltonian_values(ir+iofset),ir=1,nelc), (&
                        matrix_partner_indices(ir + iofset),ir=1,nelc)
                     hamiltonian_values(nelc+iofset) = hamiltonian_values(nelc+iofset) - hamiltonian_energy_shift
                     diag(i) = hamiltonian_values(nelc+iofset)
                     iofset = iofset + nelc
                     stored_segment_end_offsets(i) = iofset
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
                  do i = myid + 1, csf_count, nprocs
                     ndense_l = ndense_l + i
                  end do

                  call alloc (hamiltonian_values, ndense_l, 'EMT', 'MANEIG')
                  call dinit (ndense_l, 0.0_real64, hamiltonian_values, 1)
                  call alloc (work, csf_count, 'WORK', 'MANEIG')
                  call alloc (matrix_partner_indices, csf_count, 'IROW', 'MANEIG')

                  read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
                  if (csf_count/=ncfdum .or. myid/=myiddum .or. nprocsdum/=nprocs) &
                     stop 'maneig:4'

                  iofset = 0
                  do i = myid + 1, csf_count, nprocs
                     read (imcdf) nelc, elsto, (work(ir),ir=1,nelc), (matrix_partner_indices(ir),&
                        ir=1,nelc)
                     work(nelc) = work(nelc) - hamiltonian_energy_shift
                     diag(i) = work(nelc)
                     do ir = 1, nelc
                        hamiltonian_values(iofset+matrix_partner_indices(ir)) = work(ir)
                     end do
                     iofset = iofset + i
                  end do
                  call dalloc (work, 'WORK', 'MANEIG')
                  call dalloc (matrix_partner_indices, 'IROW', 'MANEIG')

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
            lim = min(csf_count,2*nvecmx + 60)
!           lwork = 2*ncf*lim + lim*lim + (nvecmx+10)*lim + nvecmx
            lwork = 2*csf_count*lim + lim*lim*2 + 11*lim + nvecmx
            call alloc (work, lwork, 'WORK', 'MANEIG')
            liwork = 6*lim + nvecmx
            call alloc (iwork, liwork, 'IWORK', 'MANEIG')
!*changed by Misha 02/12/97
            crite = 1.0e-17_real64
            critc = 1.0e-9_real64
            critr = 1.0e-9_real64
            ortho = max(1e-9_real64,critr)
! end of changes

!            maxitr = MAX (nvecmx*100, ncf/10)
            maxitr = max(nvecmx*200,csf_count/10)
            !maxitr = MIN (nvect*100, ncf)  ! FROM RSCFVU !!!
            call alloc (jwork, lim,'JWORK', 'MANEIG' )

            call alloc (shifted_eigenvalues, nvecmx, 'EVAL', 'MANEIG')
            call alloc (csf_mixing_coefficients, csf_count*nvecmx, 'EVEC', 'MANEIG')

            dmungo = 10.e99_real64
            call dinit (nvecmx, dmungo, shifted_eigenvalues, 1)

!   Compute the eigenpairs in each block

            nvex = nvecmx
            if (ldisc) then
               mblock = nvex
            else
               mblock = 1
            endif
            nend = csf_count*nvex

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
               call iniestsd (2000, csf_count, myid, nprocs, niv, work, imcdf, hamiltonian_energy_shift)

               !NIV = 0   ! Why equal 0 ???
               !WRITE (6, *) ' Calling gdvd(spodmv,...'
               call gdvd (spodmv, csf_count, lim, diag, ilow, ihigh, jwork, niv, &
                  mblock, crite, critc, critr, ortho, maxitr, work, lwork, &
                  iwork, liwork, hiend, nloops, nmv, ierr)

            case (2)
!******************** sparse and matrix in memory ********************
               write (6, *) ' Sparse - Memory,  iniest2'
!            CALL INIEST2 (1000, NCF,NIV,WORK,EMT,IENDC,IROW)
               call iniest2 (2000, csf_count, niv, work, hamiltonian_values, stored_segment_end_offsets, matrix_partner_indices)
               !WRITE (*, *) NCF, NIV, (WORK(I),I=NCF*NIV + 1,NCF*NIV + NIV)
               !WRITE (*, *) LIM, ILOW, IHIGH, MBLOCK, MAXITR, LWORK, LIWORK
               !WRITE (*, *) IERR
               call gdvd (spicmv2, csf_count, lim, diag, ilow, ihigh, jwork, niv, &
                  mblock, crite, critc, critr, ortho, maxitr, work, lwork, &
                  iwork, liwork, hiend, nloops, nmv, ierr)
               !WRITE (*, *) 'after gdvd...'
               !WRITE (*, *) NCF, NIV, (WORK(I),I=NCF*NIV + 1,NCF*NIV + NIV)
               !WRITE (*, *) LIM, ILOW, IHIGH, MBLOCK, MAXITR, LWORK, LIWORK
               !WRITE (*, *) HIEND, NLOOPS, NMV, IERR

               call dalloc (hamiltonian_values, 'EMT', 'MANEIG')
               call dalloc (matrix_partner_indices, 'IROW', 'MANEIG')
               call dalloc (stored_segment_end_offsets, 'IENDC', 'MANEIG')

            case (3)
!*************************** dense and in memory **********************
               write (6, *) ' Dense - Memory,  iniestdm'
!            CALL INIESTDM (1000,NCF,NIV,WORK,EMT)
               call iniestdm (2000, csf_count, niv, work, hamiltonian_values)
               call gdvd (dnicmv, csf_count, lim, diag, ilow, ihigh, jwork, niv, &
                  mblock, crite, critc, critr, ortho, maxitr, work, lwork, &
                  iwork, liwork, hiend, nloops, nmv, ierr)
               call dalloc (hamiltonian_values, 'EMT', 'MANEIG')
            end select
!************************************************************************
            call dalloc (diag, 'DIAG', 'MANEIG')
            call dalloc (iwork, 'IWORK', 'MANEIG')
            call dalloc (jwork, 'JWORK', 'MANEIG')

            write (24, *) ' ', nloops, ' iterations;'
            write (24, *) ' ', nmv, ' matrix-vector multiplies.'

            if (ierr /= 0) then
               write (stderr_unit, *) 'MANEIG: Returned from DVDSON with'
               write (stderr_unit, *) ' IERR = ', ierr, '.'
               stop 'maneig: DVDSON wrong'
            endif

!   Put the eigenpairs in order, overwriting as necessary

            call dcopy (nvex, work(nend+1:nend+(nvex)), 1, shifted_eigenvalues, 1)
            call dcopy (csf_count*nvex, work(1:csf_count*nvex), 1, csf_mixing_coefficients, 1)
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
                  call dswap (1, shifted_eigenvalues(iofset:(iofset)+((1)-1)), 1, shifted_eigenvalues(i:(i)+((1)-1)), 1)
                  iwork(i) = iwork(iofset)
                  iwork(iofset) = loc
                  iofset = csf_count*(iofset - 1)
                  loc = csf_count*(i - 1)
                  call dswap (csf_count, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1, &
                     & csf_mixing_coefficients(loc+1:loc+(csf_count)), 1)
               end do
               call dalloc (iwork, 'IWORK', 'MANEIG')
               call ralloc (shifted_eigenvalues,  nvec, 'EVAL', 'MANEIG')
               call ralloc (csf_mixing_coefficients, csf_count*nvec, 'EVEC', 'MANEIG' )

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

         iofset = (j - 1)*csf_count

         amax = 0._real64
         do i = 1, csf_count
            wa = abs(csf_mixing_coefficients(i+iofset))
            if (wa <= amax) cycle
            amax = wa
            ia = i
         end do

!          Find the angular momentum and parity of the dominant component

         iatjpo = itjpo(ia)
         iaspar = ispar(ia)

!          Change sign of eigenvactor if dominant component is negative

         if (csf_mixing_coefficients(ia+iofset) >= 0._real64) cycle
         csf_mixing_coefficients(1+iofset:csf_count+iofset) = -csf_mixing_coefficients(1+iofset:csf_count+iofset)
      end do

      return
      end subroutine maneig
