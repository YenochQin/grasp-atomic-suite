!***********************************************************************
!                                                                      *
      subroutine matrixmpi(dvdfirst)
!                                                                      *
!   Calls routines to form the Hamiltonian matrix and to diagonalise   *
!   it. The total angular momenta and parity of each  ASF  is found;   *
!   the eigenvectors are normalised so  that the sign of the largest   *
!   element is positive.                                               *
!                                                                      *
!   Call(s) to: [LIB92] ALLOC, DALLOC.                                 *
!               [RSCF92]: SETHAM, HMOUT.                               *
!               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL                       *
!                                                                      *
!                                         Last revision: 24 Dec 1992   *
! Block version by Xinghong He                           07 Aug 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:52:04   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double, long
      use memory_man
      use damp_c
      use def_c, only: iccmin, ncmin, ncmax
      use debug_c
      use eigv_c
      use hblock_c
      use hmat_c
      use iounit_c
      use mcpa_c
      use mpi_c
      use orb_c
      use pos_c
      use peav_c
      use syma_c,          only: jpgg
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setham_i
      use maneigmpi_i
      use suite_linear_algebra, only: dcopy, ddot, dscal
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) ::  dvdfirst
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nfile, nelec, ntmp, jblock, jblockt, i, ncfpat, ncminpat, &
         nevecpat, ncft, ncoeff, lab, ncontr, itmp, ir, idiag, j, iofset, &
         jother, ia, jofset, iattmp, iastmp, nelmntgg
      real(double), dimension(:), pointer :: cmvl
      real(double) :: amax, cdampj, dnfac, evecij, sum,  tmp, omcdaj, ovrlap, wa
      logical :: first
      character :: mcplab*3

      save first
!
      data first/ .true./
!
!     POINTER (cmvl(1))
!
!-----------------------------------------------------------------------

      if (myid == 0) write (6, *)

!   Allocate memory for CMVL once (the maximum size)
!   Save previous estimate of eigenvectors

      if (.not.first) then
         call alloc (cmvl, nvecsiz, 'CMVL', 'MATRIXmpi')
         call dcopy (nvecsiz, evec, 1, cmvl, 1)
      endif

!=======================================================================
!   Position the files - MCP files (unit NFILE) for reading
!   and mixing coefficients file (unit 25) for writing
!=======================================================================

      do nfile = 30, 32 + kmaxf
         rewind (nfile)
         if (nfile == 30) then
            read (nfile)
            read (nfile)
            read (nfile)
         endif
         read (nfile)
         read (nfile)
         read (nfile)
      end do

      ! To put in ncmin and nvecsiz. Values read here are the same
      ! as those from elsewhere (shch as common blocks)
      if (myid == 0) then
         rewind (25)
         read (25)                               ! 'G92MIX'
         read (25) nelec, ncftot, nw, ntmp, ntmp, nblock
         backspace (25)
         write (25) nelec, ncftot, nw, ncmin, nvecsiz, nblock
      endif

!=======================================================================
!   Do the job block by block
!=======================================================================

               !------------------------------------------------
      do jblock = 1, nblock                      ! block do-loop
               !------------------------------------------------

!=======================================================================
!   Read indeces of non-zero elements from mcp.30 file. Note the
!   format has been changed to lower-triangle-by-rows.
! Length of iendc can be reduced
!=======================================================================

         read (30) mcplab, jblockt, ncf
         if (jblockt/=jblock .or. ncf/=ncfblk(jblock)) error stop &
            'MATRIXMPI: jblockt .NE. jblock .OR. ncf1 .NE. ncf2'
         read (30) nelmntgg
         nelmnt = int(nelmntgg,kind=long)
         call alloc (irow, nelmnt, 'IROW', 'MATRIXmpi')
         call alloc (emt, nelmnt, 'EMT', 'MATRIXmpi')
        call alloc (iendc, 0, ncf, 'IENDC', 'MATRIXMPI' )

!                                                ! may not be necessary if iendc is ALWAYS used
!GG         IENDC(:NCF) = 0                        ! the way it is assigned here.
         iendc(0:ncf) = 0

      !...EMT will be accumulated in setham
         emt(:nelmnt) = 0.d0

         read (30) (iendc(i),i=myid+1,ncf,nprocs), (irow(i),i=1,nelmnt)

         ncfpat = ncfpast(jblock)
         ncminpat = ncminpast(jblock)
         nevecpat = nevecpast(jblock)

!=======================================================================
!   Skip current block if no eigenlaue is required
!=======================================================================

         if (nevblk(jblock) == 0) then
            do nfile = 31, 32 + kmaxf
               read (nfile) mcplab, jblockt, ncft, ncoeff
               if (jblockt /= jblock)                  &
                          error stop 'MATRIXmpi: jblockt .NE. jblock'
               if (ncft /= ncf) error stop 'MATRIXmpi: ncft .NE. ncf'

               read (nfile) lab, ncontr
               do while(lab/=0 .or. ncontr/=0)
                  read (nfile) (itmp,itmp,tmp,i=1,ncontr)
                  read (nfile) lab, ncontr
               end do
            end do

            call dalloc (iendc, 'IENDC', 'MATRIXmpi')
            call dalloc (irow, 'IROW', 'MATRIXmpi')

            cycle

         endif

!=======================================================================
!   Generate the Hamiltonian matrix - average energy is removed here
!=======================================================================

         call setham (jblock, myid, nprocs)
!
!   Determine average energy
!
         eav = 0.d0
         do ir = myid + 1, ncf, nprocs
            eav = eav + emt(iendc(ir))
         end do
         call mpi_barrier(mpi_comm_world,ierr)
         call mpi_allreduce (eav,tmp,1,mpi_double_precision,mpi_sum,  &
                             mpi_comm_world,ierr)
         eav = tmp / ncf
!GG         EAV = EAV/NCF
         eavblk(jblock) = eav

      ! Print Hamiltonian matrix and average energy
      ! hmout is not general
      !call hmout (0, 1, ncf)

         if (myid == 0) write (*, 302) eav

      ! Subtract the average energy from the diagonal elements
      ! to reduce the condition number of the matrix
!         DO I = 1, (NCF - (MYID + 1) + NPROCS)/NPROCS
!            EMT(IENDC(NPROCS*(I-1)+MYID+1)) = EMT(IENDC(NPROCS*(I-1)+MYID+1))&
!                - EAV
         do i = myid + 1, ncf, nprocs
            idiag = iendc(i)     ! new mode: each row ends in diagonal
            emt(idiag) = emt(idiag) - eav
         end do

!=======================================================================
!   Compute and store eigenpairs
!=======================================================================
         call maneigmpi (dvdfirst, ldbpg(3),                      &
                   jblock, ncfpat, ncminpat, nevecpat, ncftot)

!=======================================================================
!   Damp and Schmidt orthogonalise eigenvectors for OL calculations
!=======================================================================

         if (.not.first) then

            do j = 1, nevblk(jblock)

               iofset = (j - 1)*ncf + nevecpat
               jother = j

            ! cdamp has the original non-block feature
               cdampj = cdamp(j + ncminpat)
               if (cdampj == 0.d0) cycle         ! So SURE ???

               omcdaj = 1.d0 - cdampj

           !...Damp eigenvector and determine the new dominant component
  123          continue
               amax = 0.d0
               do i = 1, ncf
                  evecij = omcdaj*evec(i+iofset) + cdampj*cmvl(i + iofset)
                  evec(i+iofset) = evecij
                  wa = abs(evecij)
                  if (wa <= amax) cycle
                  amax = wa
                  ia = i
               end do

            !...compute the normalization factor
               sum = 0.d0
               do i = 1, ncf
                  sum = sum + evec(i+iofset)**2
               end do
               dnfac = 1.d0/sqrt(sum)

            !...Renormalize and invert as necessary
               if (evec(ia+iofset) < 0.d0) dnfac = -dnfac
               call dscal (ncf, dnfac, evec(iofset+1:iofset+(ncf)), 1)

            !...Schmidt orthogonalise
  234          continue
               jother = jother - 1
               if (jother < 1) cycle
               jofset = (jother - 1)*ncf + nevecpat
               ovrlap = ddot(ncf,evec(iofset+1:iofset+(ncf)),1,evec(jofset+1:jofset+(ncf)),1)
               if (ovrlap /= 0.d0) then          ! So SURE ???
                  omcdaj = 1.d0
                  cdampj = -ovrlap
                  call dcopy (ncf, evec(jofset+1:jofset+(ncf)), 1, cmvl(iofset + 1:iofset+(ncf)), 1)
                  go to 123
               else
                  go to 234
               endif
            end do
         endif

!   Write out the eigenpair information: ASF symmetries, eigenvalues,
!   and eigenvectors to GRASP92 mixing coefficients File

         if (nevblk(jblock) == 0) then
            iattmp = 999
            iastmp = 999
         else
!GGGG
         iattmp = iabs(jpgg(jblock))
         if(jpgg(jblock) .ge. 0) then
            iastmp = 1
         else
            iastmp = -1
         end if
!GG            IATTMP = IATJPO(NCMINPAT + 1)
!GG            IASTMP = IASPAR(NCMINPAT + 1)
         endif

         if (myid == 0) then
            write (25) jblock, ncf, nevblk(jblock), iattmp, iastmp
            write (25) (iccmin(i + ncminpat),i=1,nevblk(jblock))
            write (25) eav, (eval(i + ncminpat),i=1,nevblk(jblock))
            write (25) ((evec(i+(j-1)*ncf+nevecpat),i=1,ncf),j=1,nevblk(jblock)&
               )
         endif

         call dalloc (emt, 'EMT', 'MATRIXmpi')
         call dalloc (iendc, 'IENDC', 'MATRIXmpi')
         call dalloc (irow, 'IROW', 'MATRIXmpi')

!----------------------
      end do
!----------------------
!
!   Deallocate the temporary storage
!
      if (.not.first) then
         call dalloc (cmvl, 'CMVL', 'MATRIXmpi')
      else
         first = .false.
      endif

  302 format(' Average energy = ',1p,d18.10,' Hartrees')

      return
      end subroutine matrixmpi
