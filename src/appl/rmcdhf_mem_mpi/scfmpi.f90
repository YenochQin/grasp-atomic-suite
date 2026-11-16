!***********************************************************************
!                                                                      *
      subroutine scfmpi(eol, rwffile2)
!                                                                      *
!   This  subroutine  performs  the SCF iterations. The procedure is   *
!   essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3   *
!   (1986) 290.                                                        *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
!               [RSCF92]: improvmpi, matrixmpi, MAXARR, newcompi,      *
!                         ORBOUT, ORTHSC, setlagMpi.                   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
!   ifort -i8 version by Alexander Kramida (AK) Last rev. 22 Mar 2016  *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:16:00   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use blkidx_c
      use default_c
      use def_c
      use debug_c
      use eigv_c
      use fixd_c
      use hblock_c
      use iounit_c
      use lagr_c
      use mcpa_c
      use mpi_c
      use pos_c
      use peav_c
      use orb_c
      use orba_c
      use scf_c
      use syma_c
      use stat_c
      use orthct_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use matrixmpi_i
      use newcompi_i
      use setlagmpi_i
      use improvmpi_i
      use maxarr_i
      use prwf_i
      use orthsc_i
      use orbout_i
      use endsum_i
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical  :: eol
      character  :: rwffile2*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, nit, jseq, kount, k, l_mpi
      real(double) :: wtaev, wtaev0, dampmx
      logical :: convg, lsort, dvdfirst
!-----------------------------------------------
!CFF   .. set the logical variable dvdfirst
      dvdfirst = .true.
!-----------------------------------------------------------------------
!AK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
!      ISIZE = sizeof(NCF)
      l_mpi = mpi_logical
!      if (ISIZE.EQ.8) L_MPI = MPI_LOGICAL8

      ncftot = ncf
      !IF (myid .EQ. 0) PRINT *, '===SCF==='

!=======================================================================
!   Determine Orthonomalization order --- lsort
!=======================================================================

      if (ndef == 0) then
         lsort = .false.
      else
         if (myid .eq. 0) then
  123       continue
            write (istde, *) 'Orthonomalization order? '
            write (istde, *) '     1 -- Update order'
            write (istde, *) '     2 -- Self consistency connected'
            read (istdi, *) j
            if (j == 1) then
               lsort = .false.
            else if (j == 2) then
               lsort = .true.
            else
               write (istde, *) 'Input is wrong, redo...'
               go to 123
            endif
         endif
         call mpi_bcast (lsort, 1, l_mpi, 0, mpi_comm_world, ierr)
      endif

!=======================================================================
!   Deallocate storage that will no longer be used
!=======================================================================

!GG      CALL DALLOC (JQSA, 'JQSA', 'SCFmpi')

!=======================================================================
!   Allocate and fill in auxiliary arrays
!=======================================================================

      call alloc (ncfpast, nblock, 'NCFPAST', 'SCFmpi')
      call alloc (ncminpast, nblock, 'NCMINPAST', 'SCFmpi')
      call alloc (nevecpast, nblock, 'NEVECPAST', 'SCFmpi')
      call alloc (eavblk, nblock, 'EAVBLK', 'SCFmpi')

      ncfpast(1) = 0
      ncminpast(1) = 0
      nevecpast(1) = 0
      do i = 2, nblock
         ncfpast(i) = ncfpast(i-1) + ncfblk(i - 1)
         ncminpast(i) = ncminpast(i-1) + nevblk(i - 1)
         nevecpast(i) = nevecpast(i-1) + nevblk(i - 1)*ncfblk(i - 1)
      end do

      !*** Size of the eigenvector array for all blocks
      nvecsiz = nevecpast(nblock) + nevblk(nblock)*ncfblk(nblock)

      if (eol) then
         call alloc (eval, ncmin, 'EVAL', 'SCFmpi')
         call alloc (evec, nvecsiz, 'EVEC', 'SCFmpi')
         call alloc (iatjpo, ncmin, 'IATJPO', 'SCFmpi')
         call alloc (iaspar, ncmin, 'IASPAR', 'SCFmpi')
      endif

!=======================================================================
!
!=======================================================================
      nddim = 64
      nxdim = 64
      nydim = 64
!      NDDIM = 0
!      NYDIM = 0
!      NYDIM = 0
!AK    Must initialize DA,XA,YA,... here since MPI routines rely on their existence in all processes
      call alloc (nda, nddim, 'NDA','SCFmpi')
      call alloc (da, nddim, 'DA','SCFmpi')
      call alloc (nxa, nxdim, 'NXA','SCFmpi')
      call alloc (xa, nxdim, 'XA','SCFmpi')
      call alloc (nya, nydim, 'NYA','SCFmpi')
      call alloc (ya, nydim, 'YA','SCFmpi')

!     This call should only be made AFTER the call to newco
!     CALL setlag (EOL)

!   For (E)OL calculations, determine the level energies and
!   mixing coefficients

!CFF   .. set the logical variable dvdfirst
      if (eol) then
         call matrixmpi (dvdfirst)
         call newcompi (wtaev)
      endif
      wtaev0 = 0.0
      dvdfirst = .false.
      do nit = 1, nscf
         if (myid == 0) write (*, 301) nit

!   For all pairs constrained through a Lagrange multiplier, compute
!   the Lagrange multiplier

         call setlagmpi (eol)

!   Improve all orbitals in turn

         dampmx = 0.0
         if (myid == 0) write (*, 302)
         do j = 1, nw
            jseq = iorder(j)
            if (lfix(jseq)) cycle
            call improvmpi (eol, jseq, lsort, dampmx)
         end do
!
!   For KOUNT = 1 to NSIC: find the least self-consistent orbital;
!   improve it
!
!         write(istde,*) 'nsic=',nsic
         do kount = 1, nsic
            call maxarr (k)
            if (k == 0) then
               convg = .true.
               go to 3
            else
               if (scnsty(k) <= accy) then
                  convg = .true.
                  go to 3
               endif
            endif
            call improvmpi (eol, k, lsort, dampmx)
         end do

         call maxarr (k)

         if (k == 0) then
            convg = .true.
         else
            if (scnsty(k) <= accy) then
               convg = .true.
            else
               convg = .false.
            endif
         endif

    3    continue
         if (ldbpr(24) .and. myid==0) call prwf (0)

!   Perform Gram-Schmidt process
!   For OL calculation, orthst is true and orbitals are orthonormalized
!   in subroutine improv. For AL calculation, orthst is false.
         if (.not.orthst) call orthsc

!   Write the subshell radial wavefunctions to the .rwf file

         if (myid == 0) call orbout (rwffile2)

         if (eol) then
            call matrixmpi(dvdfirst)
            call newcompi (wtaev)
         endif
!        Make this a relative convergence test
!        IF(ABS(WTAEV-WTAEV0).LT.1.0D-9.and.
!    &   DAMPMX.LT.1.0D-4) CONVG=.true.
!        PRINT *, 'WTAEV, WTAEV0', WTAEV, WTAEV0,
!    &             ABS((WTAEV-WTAEV0)/WTAEV)
!GG         IF (ABS((WTAEV - WTAEV0)/WTAEV) < 0.001*ACCY) CONVG = .TRUE.
!cjb unified convergence criteria in RMCDHF and RMCDHF_MPI
!cjb     IF(DABS(WTAEV-WTAEV0).LT.1.0D-8.and.              &
!cjb                   DAMPMX.LT.1.0D-2) CONVG=.true.
         if (abs((wtaev - wtaev0)/wtaev) < 0.001*accy) convg = .true.
         wtaev0 = wtaev
         if (.not.convg) cycle
         if (ldbpr(25) .and. .not.ldbpr(24) .and. myid==0) call prwf (0)
            !IF (EOL) CALL matrixmpi (dvdfirst)
         go to 5

      end do

      if (myid==0) write (istde,*)' Maximum iterations in SCF Exceeded.'

    5 continue
      do i = 31, 32 + kmaxf
         close(i)                                ! The MCP coefficient files
      end do

      if (myid == 0) then
         !CLOSE (23)     ! The .rwf file
         close(25)                               ! The .mix file
      endif
!
!   Complete the summary - moved from rscf92 for easier alloc/dalloc
!
      if (myid .eq. 0) call endsum
!
!   Deallocate storage
!
      call dalloc (wt, 'WT', 'SCFmpi')                       !Either getold or getald

      if (nec > 0) then
         call dalloc (iecc, 'IECC', 'SCFmpi')
         call dalloc (ecv, 'ECV', 'SCFmpi')
         call dalloc (iqa, 'IQA', 'SCFmpi')
      endif

      if (nddim > 0) then
         call dalloc (da, 'DA', 'SCFmpi')
         call dalloc (nda, 'NDA', 'SCFmpi')
         nddim = 0
      endif

      if (nxdim > 0) then
         call dalloc (xa, 'XA', 'SCFmpi')
         call dalloc (nxa, 'NXA', 'SCFmpi')
         nxdim = 0
      endif

      if (nydim > 0) then
         call dalloc (ya, 'YA', 'SCFmpi')
         call dalloc (nya, 'NYA', 'SCFmpi')
         nydim = 0
      endif

      if (eol) then
         call dalloc (eval, 'EVAL', 'SCFmpi')
         call dalloc (evec, 'EvEC', 'SCFmpi')
         call dalloc (iatjpo, 'IATJPO', 'SCFmpi')
         call dalloc (iaspar, 'IASPAR', 'SCFmpi')
         call dalloc (ncmaxblk, 'NCMAXBLK', 'SCFmpi') ! getold.f
         call dalloc (eavblk, 'EAVBLK', 'SCFmpi')     ! getold.f
         call dalloc (idxblk, 'IDXBLK', 'SCFmpi')     ! Allocated in getold.f
         call dalloc (iccmin, 'ICCMIN', 'SCFmpi')       ! Allocated in items.f<-getold.f
      endif
!
      call dalloc (ncfpast, 'NCFPAST', 'SCFmpi')
      call dalloc (ncminpast, 'NCMINPAST', 'SCFmpi')
      call dalloc (nevecpast, 'NEVECPAST', 'SCFmpi')

  301 format(/,' Iteration number ',1i3,/,' --------------------')
  302 format(41x,'Self-            Damping'/,&
         'Subshell    Energy    Method   P0    ',&
         'consistency  Norm-1  factor  JP',' MTP INV NNP'/)

      return
      end subroutine scfmpi
