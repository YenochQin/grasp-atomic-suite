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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use memory_man
      use eigv_c
      use fposition_c, only: nposition
      use hmat_c
      use orb_c,           only: csf_count, orbital_count, csf_subshell_occupations
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use where_c, only: double, imcdf
      use wchblk_c, only: double, jblock
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dnicmv_i, only: dnicmv
      use spicmv2_i, only: spicmv2
      use spodmv_i, only: spodmv
      use posfile_i, only: posfile
      use dinit_i, only: dinit
      use dspevx_i
      use iniestsd_i, only: iniestsd
      use gdvd_i
      use iniest2_i, only: iniest2
      use iniestdm_i, only: iniestdm
      use itjpo_i, only: itjpo
      use ispar_i, only: ispar
      use spicmvmpi_i, only: spicmvmpi
      use suite_linear_algebra, only: dcopy, dspevx, dswap
      use iniestmpi_i, only: iniestmpi
      use mpiu_i, only: gdsummpi
      use mpiu_i, only: stopmpi
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
      integer(int64)        :: nelmnt_a
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
      integer(int64) :: nstore_s, nstore_d
      integer :: nrows, i, ndense, ncfdum, iccutdum, myiddum, &
         nprocsdum, iofset, nelc, ir, nvecmn, nvex, m, info, loc, nbrkev,   &
         imv, ndense_l, lim, lwork, liwork, maxitr, mblock, nend, &
         ilow, ihigh, niv, nloops, nmv, j, ia,   &
         idummy
      real(real64) :: elsto,  dummy, &
          diatmp, crite, critc, critr, ortho, dmungo, amax, wa, abstol
! GG          DIATMP, CRITE, CRITC, CRITR, ORTHO, DMUNGO, AMAX, WA
      logical :: hiend, ldisc, sparse
      character(len=8) :: cnum
      real(real64), pointer :: w(:), z(:), work(:), diag(:)
      integer, pointer :: iwork(:), ifail(:), jwork(:)
!-----------------------------------------------
!
!-----------------------------------------------------------------------
      abstol = 2*dlamch("S")
      if (mpi_rank == 0) write (6, *) "Calling maneig..."

!  (nrows+1) is the number of records of the present block's .res file

      nrows = (csf_count - mpi_rank - 1 + mpi_process_count)/mpi_process_count
      if (csf_count < mpi_process_count) nrows = csf_count/(mpi_rank + 1)
      !CALL posfile (1, imcdf, nrows+1)
      call posfile (0, imcdf, nposition)

      if (csf_count == 1) then
!-----------------------------------------------------------------------
!
! (1) - Trivial    ncf = 1
!
!-------------------------------------------------------
         if (mpi_rank == 0) write (24, *) "Trivial eigenvalue problem."

!   Matrix of order 1: the trivial case; we assume that the value
!   of EAV is available

         call alloc (shifted_eigenvalues, 1,"EVAL", "MANEIG" )
         call alloc (csf_mixing_coefficients, 1, "EVECO", "MANEIG")
         shifted_eigenvalues(1) = 0.0_real64
         csf_mixing_coefficients(1) = 1.0_real64

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
         if(mpi_rank==0) write(*,*) "ncf=",csf_count, " iolpck=", iolpck
         if (csf_count <= iolpck) then
!-----------------------------------------------------------------------
!
! (2.1) - LAPACK    Dense, Memory,
!
!-------------------------------------------------------
            if (mpi_rank == 0) then
               write (6, *) &
                  "LAPACK routine DSPEVX selected for eigenvalue problem."
               write (24, *) &
                  "LAPACK routine DSPEVX selected for eigenvalue problem."
            end if

!   Allocate storage for the dense representation of the matrix
!   and initialize emt

            ndense = (csf_count*(csf_count + 1))/2
            call alloc (hamiltonian_values, ndense, "EMT", "MANEIG")
            call dinit (ndense, 0.0_real64, hamiltonian_values, 1)

!   Read the matrix into position from the disc file; it's already
!   been properly positioned.

            call alloc (work, csf_count,"WORK", "MANEIG" )
            call alloc (matrix_partner_indices, csf_count, "IROW", "MANEIG")
            read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
            if (csf_count/=ncfdum .or. mpi_rank/=myiddum .or. nprocsdum/=mpi_process_count) stop &
               "maneig:1"

            do i = mpi_rank + 1, csf_count, mpi_process_count
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

            call dalloc (work, "WORK", "MANEIG")
            call dalloc (matrix_partner_indices, "IROW", "MANEIG")

! Let each node have a complete copy of EMT

            call gdsummpi (hamiltonian_values, ndense)

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
            call alloc (w, nvex, "W", "MANEIG")
            call alloc (z, csf_count*nvex,"Z", "MANEIG" )
            call alloc (work, csf_count*8,"WORK", "MANEIG" )
            call alloc (iwork, csf_count*5,"IWORK", "MANEIG" )
! GG            CALL ALLOC (IFAIL, NVEX, 'IFAIL', 'MANEIG')
            call alloc (ifail, csf_count, "IFAIL", "MANEIG")
            call dspevx ("V", "I", "U", csf_count, hamiltonian_values, dummy, dummy, nvecmn, nvecmx&
               , abstol, m, w, z, csf_count, work, iwork, ifail, info)
            if (info /= 0)                                             &
                call stopmpi("maneig: Failure in DSPEVX [LAPACK]",mpi_rank)
            call dalloc (work, "WORK", "MANEIG")
            call dalloc (iwork, "IWORK", "MANEIG")
            call dalloc (ifail, "IFAIL", "MANEIG")
            call dalloc (hamiltonian_values, "EMT", "MANEIG")

!   Store the eigenpairs in their proper positions EVAL() and EVEC()

            call alloc (shifted_eigenvalues, nvec,"EVAL", "MANEIG")
            call alloc (csf_mixing_coefficients, csf_count*nvec, "EVEC", "MANEIG")

            do i = 1, nvec
               loc = ivec(i)
               shifted_eigenvalues(i) = w(loc - nvecmn + 1)
               iofset = csf_count*(i - 1)
               loc = csf_count*(loc - nvecmn)
               call dcopy (csf_count, z(loc + 1:loc+(csf_count)), 1, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)
            end do
            call dalloc (w, "W", "MANEIG")
            call dalloc (z, "Z", "MANEIG")

         else
!-----------------------------------------------------------------------

! (2.2) - DVDSON --- preparation work
!
!-------------------------------------------------------
            if (mpi_rank == 0) write (24,*)"DVDSON routine selected" &
                                      //" for eigenvalue problem;"
!--------------------------------------------------------------
            if (mpi_rank == 0) print *,"zou: nelmnt",stored_matrix_element_count
!            print *, 'myid, nprocs,nincor', myid, nprocs, nincor

!           If on any process more memory is needed than NINCOR, use disc.
!           If not, should the matrix be sparse or dense (distributed)?
!           sparse requires = nelmnt + (isize*(nelmnt+ncf+1))/8
!           dense  storage = (ncf*(ncf+1)/2/nproc - on average
!           (because of integer arithmetic nproc-1 is added to numerator

!           This is a comparion of matrix elements only.
!
            nstore_s = nelmnt_a
            nstore_d = ((csf_count+mpi_process_count-1)/mpi_process_count)
            nstore_d = nstore_d*((csf_count+1)/2)
            if (mpi_rank == 0) write(*,*) "nstore_s=", nstore_s,   &
                             "nstore_d=", nstore_d, "nincor=", nincor
            ldisc = .false.
            if ( nstore_s > nincor) then
!              .. we need to use disk
               ldisc =.true.
               sparse=.true.
            else

               if (nstore_s < nstore_d ) then
                  sparse = .true.
               else
                  sparse = .false.
               end if
            end if

            call alloc (diag,csf_count,"DIAG", "MANEIG")
            call dinit (csf_count, 0.0_real64, diag, 1)
            if (ldisc) then
!-----------------------------------------------------------------------
!
! (2.2.1) - DVDSON --- Disk, load diagonal
!
!-------------------------------------------------------
               if (mpi_rank == 0)  then
                  write (*, *) " matrix stored on disc;"
                  write (24, *) " matrix stored on disc;"
               end if

!   Disk storage; necessarily sparse; one column of the matrix in
!   memory

               imv = 1

!   Load diagonal - Each node will have the same, complete copy
!   after this if block
               print *, mpi_rank, "Reading diagonal, unit imcdf =", imcdf
               read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
               if (csf_count/=ncfdum .or. mpi_rank/=myiddum .or. nprocsdum/=mpi_process_count) stop &
                  "maneig:2"

               do i = mpi_rank + 1, csf_count, mpi_process_count
                  read (imcdf) nelc, elsto, (dummy,ir=2,nelc), diatmp, &
                                            (idummy, ir=1,nelc)
                  diag(i) = diatmp - hamiltonian_energy_shift
               end do
               print *, " finished"
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
                  if (mpi_rank == 0) then
                     write (24, *) " sparse matrix stored in memory"
                     write ( *, *) " sparse matrix stored in memory"
                  end if

                  imv = 2
                  write (6, *) "nelmnt = ", stored_matrix_element_count
                  call alloc (hamiltonian_values, stored_matrix_element_count, "EMT", "MANEIG")
                  call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "MANEIG")
                  call alloc (stored_segment_end_offsets, 0, csf_count , "IENDC", "MANEIG")
                  iofset = 0
                  stored_segment_end_offsets(0) = 0

                  read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
                  if (csf_count/=ncfdum .or. mpi_rank/=myiddum .or. nprocsdum/=mpi_process_count) &
                     stop "maneig:3"
                  do i = mpi_rank + 1, csf_count, mpi_process_count
                     read (imcdf) nelc, elsto, (hamiltonian_values(ir+iofset),ir=1,nelc), (&
                        matrix_partner_indices(ir + iofset),ir=1,nelc)
                     hamiltonian_values(nelc+iofset) = hamiltonian_values(nelc+iofset) - hamiltonian_energy_shift
                     diag(i) = hamiltonian_values(nelc+iofset)
                     iofset = iofset + nelc
                     stored_segment_end_offsets(i) = iofset
                  end do
               else
!-----------------------------------------------------------------------
!
! (2.2.2.2) - DVDSON --- Memory, load all, dense
!
!-------------------------------------------------------
                  if (mpi_rank == 0) then
                     write (24,*) " full matrix stored in in memory"
                     write ( *,*) " full matrix stored in in memory"
                  end if

                  imv = 3

! Find NDENSE_L, the number of elements on the node (dense form)

                  ndense_l = 0
                  do i = mpi_rank + 1, csf_count, mpi_process_count
                     ndense_l = ndense_l + i
                  end do

                  call alloc (hamiltonian_values, ndense_l, "EMT", "MANEIG")
                  call dinit (ndense_l, 0.0_real64, hamiltonian_values, 1)
                  call alloc (work, csf_count, "WORK", "MANEIG")
                  call alloc (matrix_partner_indices, csf_count, "IROW", "MANEIG")

                  read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
                  if (csf_count/=ncfdum .or. mpi_rank/=myiddum .or. nprocsdum/=mpi_process_count) &
                     stop "maneig:4"

                  iofset = 0
                  do i = mpi_rank + 1, csf_count, mpi_process_count
                     read (imcdf) nelc, elsto, (work(ir),ir=1,nelc), (matrix_partner_indices(ir),&
                        ir=1,nelc)
                     work(nelc) = work(nelc) - hamiltonian_energy_shift
                     diag(i) = work(nelc)
                     do ir = 1, nelc
                        hamiltonian_values(iofset+matrix_partner_indices(ir)) = work(ir)
                     end do
                     iofset = iofset + i
                  end do
                  call dalloc (work, "WORK", "MANEIG")
                  call dalloc (matrix_partner_indices, "IROW", "MANEIG")

               end if
!               ...Memory mode - sparse or dense
!-----------------------------------------------------------------------
!  (2.2.2.3e)          *** E n d   o f   D V D S O N   m e m o r y
!-----------------------------------------------------------------------
            end if
!            ...Disk or Memory
!-----------------------------------------------------------------------
!  (2.2.3e)      *** E n d   o f   D V D S O N
!-----------------------------------------------------------------------

! Make diagonals global, no matter it is disk or memory mode

            call gdsummpi (diag, csf_count)
!
!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECMX
!   plus a smaller number if storage is severely constrained
!
            lim = min(csf_count,2*nvecmx + 60)
!           lwork = 2*ncf*lim + lim*lim + (nvecmx+10)*lim + nvecmx
            lwork = 2*csf_count*lim + lim*lim*2 + 11*lim + nvecmx
            call alloc (work, lwork, "WORK", "MANEIG")
            work(1:lwork) = 0.0_real64
            liwork = 6*lim + nvecmx
            call alloc (iwork, liwork, "IWORK", "MANEIG")
!*changed by Misha 02/12/97
            crite = 1.0e-17_real64
!GG            CRITC = 1.0D-08
!GG            CRITR = 1.0D-08
!GG            ORTHO = MAX(1D-8,CRITR)
            critc = 1.0e-9_real64
            critr = 1.0e-9_real64
            ortho = max(1e-9_real64,critr)
! end of changes

!            maxitr = MAX (nvecmx*100, ncf/10)
            maxitr = max(nvecmx*100,csf_count/10)
            !maxitr = MIN (nvect*100, ncf)  ! FROM RSCFVU !!!
            call alloc (jwork, lim,"JWORK", "MANEIG" )

            call alloc (shifted_eigenvalues, nvecmx, "EVAL", "MANEIG")
            call alloc (csf_mixing_coefficients, csf_count*nvecmx, "EVEC", "MANEIG")

            dmungo = 10.0e99_real64
            call dinit (nvecmx, dmungo, shifted_eigenvalues, 1)

!   Compute the eigenpairs in each block

            nvex = nvecmx
            if (ldisc) then
               mblock = nvex
            else
               mblock = 1
            end if
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
              if (mpi_rank == 0) print *, " Sparse - Disk, iniestsd"
               call posfile (0, imcdf, nposition)! was within iniestsd before
               call iniestsd (iolpck, csf_count, mpi_rank, mpi_process_count, niv, work, imcdf, hamiltonian_energy_shift)
              print *, "Returned from iniestsd "
              if (csf_count> iolpck) then
                print *, "Calling GDVD"
                call gdvd (spodmv,csf_count,lim,diag,ilow,ihigh,            &
                  jwork,niv,mblock,crite,critc, critr,ortho,maxitr,   &
                  work,lwork,iwork,liwork,hiend,nloops,               &
                  nmv,mpi_error_code)
              end if

            case (2)
!******************** sparse and matrix in memory ********************
               if (mpi_rank == 0) print *, " Sparse - Memory, iniestmpi"
               call iniestmpi (iolpck, csf_count,niv,work,hamiltonian_values,stored_segment_end_offsets,matrix_partner_indices)
               if(csf_count> iolpck) then
                 call gdvd (spicmvmpi,csf_count,lim,diag,ilow,ihigh,        &
                  jwork,niv,mblock,crite,critc, critr,ortho,maxitr,   &
                  work,lwork,iwork,liwork,hiend,nloops,               &
                  nmv,mpi_error_code)
               end if

               call dalloc (hamiltonian_values, "EMT", "MANEIG")
               call dalloc (matrix_partner_indices, "IROW", "MANEIG")
               call dalloc (stored_segment_end_offsets, "IENDC", "MANEIG")

            case (3)
!*************************** dense and in memory **********************
               if (mpi_rank == 0) print *, " Dense - Memory, iniestdm"
              call iniestdm (iolpck,csf_count,niv,work,hamiltonian_values)
              if (csf_count> iolpck) then
                call gdvd (dnicmv,csf_count,lim,diag,ilow,ihigh,            &
                    jwork,niv,mblock,crite,critc, critr,ortho,maxitr, &
                    work,lwork,iwork,liwork,hiend,nloops,             &
                    nmv,mpi_error_code)
              end if
               call dalloc (hamiltonian_values, "EMT", "MANEIG")
            end select
!************************************************************************
            call dalloc (diag, "DIAG", "MANEIG")
            call dalloc (iwork, "IWORK", "MANEIG")
            call dalloc (jwork, "JWORK", "MANEIG")
            if (mpi_rank == 0) then
               write (24, *) " ", nloops, " iterations;"
               write (24, *) " ", nmv, " matrix-vector multiplies."
            end if
            if (mpi_error_code /= 0) then
               write (stderr_unit, *) "MANEIG: Returned from DVDSON with"
               write (stderr_unit, *) " IERR = ", mpi_error_code, "."
               stop "maneig: DVDSON wrong"
            end if

!   Put the eigenpairs in order, overwriting as necessary

            call dcopy (nvex, work(nend+1:nend+(nvex)), 1, shifted_eigenvalues, 1)
            call dcopy (csf_count*nvex, work(1:csf_count*nvex), 1, csf_mixing_coefficients, 1)
            call dalloc (work, "WORK", "MANEIG")

!   Rearrange and reallocate storage for the eigenpairs
!   as necessary

            if (nvec < nvecmx) then
               call alloc (iwork, nvecmx, "IWORK", "MANEIG")
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
               call dalloc (iwork, "IWORK", "MANEIG")
               call ralloc (shifted_eigenvalues,  nvec, "EVAL", "MANEIG")
               call ralloc (csf_mixing_coefficients, csf_count*nvec, "EVEC", "MANEIG" )

            end if

         end if
! (2.3e)              *** E N D   O F    N O N - T R I V I A L   C A S E

      end if
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

         amax = 0.0_real64
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

         if (csf_mixing_coefficients(ia+iofset) >= 0.0_real64) cycle
         csf_mixing_coefficients(1+iofset:csf_count+iofset) = -csf_mixing_coefficients(1+iofset:csf_count+iofset)
      end do

      return
      end subroutine maneig
