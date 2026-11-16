!***********************************************************************
!                                                                      *
      subroutine maneigmpi(dvdfirst, lprint, jblock, &
                        ncfpat, ncminpat, nevecpat, ncftot)
!                                                                      *
!   This module  manages the  operation of the  eigensolvers and the   *
!   storage of the eigenpairs.  There are two principal branches:      *
!                                                                      *
!      (1) Matrix of order 1: the trivial case                         *
!      (2) Matrix of order greater than 1: eigenpairs are found        *
!          using DVDSON; this involves up to three steps:              *
!                    (a) The matrix is analysed to determine its       *
!                        block structure (only irreducibe matrices     *
!                        are correctly treated by DVDSON)              *
!                    (b) Eigenpairs are extracted for each block       *
!                    (c) The appropriate eigenpairs are selected and   *
!                        stored                                        *
!                                                                      *
!   We  assume that  the sparse representation  of the matrix  is in   *
!   core.                                                              *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, RALLOC.          *
!               [RSCF92]: POSNFL, SPICMV2.
!               [DVDSON]: DVDSON                                       *
!               [AUXBLAS]: DINIT/SINIT                                 *
!               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL, DSWAP/SSWAP          *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
!   Modified by Xinghong He               Last revision: 17 Aug 1998   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use def_c
      use eigv_c
      use hblock_c
      use hmat_c
      use mpi_c
      use orb_c
      use syma_c,          only: jpgg
      use wchblk_c, jblockk=>jblock
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use spicmvmpi_i
      use iniestmpi_i
      use gdvd_i
      use itjpo_i
      use ispar_i
      use suite_linear_algebra, only: dcopy, dscal
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(in) ::  dvdfirst
      integer             :: jblock
      integer, intent(in) :: ncfpat
      integer, intent(in) :: ncminpat
      integer, intent(in) :: nevecpat
      integer             :: ncftot
      logical             :: lprint
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nvect, lim, lwork, nvex, niv, maxitr, n1000, mblock, &
         ilow, ihigh, liwork, ic, nloops, nmv, nend, j, jstate, &
         iofset, i, ia
      integer, dimension(:), pointer :: iwork, jwork
      real(double) :: pnwork, crite, critc, critr, ortho, amax, wa, dnfac
      real(double), dimension(:), pointer :: work, diag, atmp
      logical :: hiend
!-----------------------------------------------

      !PRINT *, 'maneig ...'

!      ...spicmvmpi needs this COMMON /WCHBLK/JBLOCKK
      jblockk = jblock
!
!=======================================================================
!   Trivial case
!=======================================================================
      if (csf_count == 1) then
         shifted_eigenvalues(ncminpat+1) = 0.d0
         csf_mixing_coefficients(nevecpat+1) = 1.d0
         go to 123                               ! Don't like big ELSE
      endif
!
!=======================================================================
!   Non-trivial case - Use Davidson eigensolver
!=======================================================================
!
!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECT
!   plus a smaller number if storage is severely constrained
!
      nvect = ncmaxblk(jblock)
!CFF   ... make this more like rscf
      lim   = min (csf_count, 2*nvect + 80)
!GG      lim   = MIN (ncf, 2*nvect + 40)
!      LIM = MIN(NCF,2*NVECT + 20)
      lwork = 2*csf_count*lim + lim*lim*2 + 11*lim + nvect
      call alloc (work, lwork, 'WORK', 'MANEIGmpi')

      !...At most 14 ? restriction removed xhh 98-05-19
      !nvex = MIN (nvect,ncfblk(jblock),14)
      nvex = min(nvect,csf_counts_by_block(jblock))
      niv = nvex
      maxitr = min(nvect*200,csf_count)
!      N1000 = 2000
      n1000 = 4000
!
!   Initial estimates for eigenvectors
!
!CFF
      if (dvdfirst .or. (csf_count .le. n1000) ) then
!GG         CALL INIESTmpi (N1000, NCF, NIV, WORK, EMT, IENDC, IROW)
         call iniestmpi (n1000, csf_count, niv, work, hamiltonian_values, stored_segment_end_offsets(0:csf_count), &
            & matrix_partner_indices)
      else
!CFF        .. use current estimates
      nend = csf_count * nvex
         do j = 1, eigenvector_counts_by_block(jblock)
            work( nend + iccmin(j+ncminpat) ) = shifted_eigenvalues(ncminpat+j)
            call dcopy (csf_count, csf_mixing_coefficients( nevecpat + csf_count*(j-1) + 1:nevecpat + csf_count*(j-1)+ &
               & (csf_count)),    1,    &
                          work( csf_count*(iccmin(j+ncminpat)-1) + 1 :( csf_count*(iccmin(j+ncminpat)-1) + 1 )+((csf_count)-1)), 1)
         enddo
      endif

! iniest looks for eigenvectors of n1000*n1000 matrix so there
! is no need to call dvdson if block size <= n1000

      if (csf_count > n1000) then
         if (mpi_rank == 0) write (*,*) 'Calling dvdson!!!', maxitr,nvect

         ! Call Davidson eigensolver

         mblock = 1
         ilow = 1
         ihigh = nvex
         liwork = 6*lim + nvect
         crite = 1.0d-17
!         CRITC = 1.0D-08
!         CRITR = 1.0D-08
!         ORTHO = MAX(1D-8,CRITR)
         critc = 1.0d-09
         critr = 1.0d-09
         ortho = max (1d-9, critr)
!
!   Store the diagonals in a separate array and make it global
!
         call alloc (diag, csf_count, 'DIAG', 'MANEIGmpi')
         call alloc (atmp, csf_count, 'ATMP', 'MANEIGmpi')
         do i = 1, csf_count
            atmp(i) = 0.d0
            diag(i) = 0.d0    ! this one may not be necessary
         enddo

         do ic = mpi_rank + 1, csf_count, mpi_process_count
!GG            DIAG(IC) = EMT(IENDC(IC))
            atmp(ic) = hamiltonian_values(stored_segment_end_offsets(ic))
         end do
         call mpi_allreduce (atmp, diag, csf_count, mpi_double_precision,  &
                             mpi_sum, mpi_comm_world, mpi_error_code)
         call dalloc(atmp,'ATMP', 'MANEIGmpi')

         call alloc (iwork, liwork, 'IWORK', 'MANEIGmpi')
         call alloc (jwork, lim, 'JWORK', 'MANEIGmpi')
         if (csf_count.gt.1000) then
            call gdvd (spicmvmpi,csf_count,lim,diag,ilow,ihigh,jwork,niv,mblock, &
            crite, critc, critr, ortho, maxitr, work, lwork, iwork, liwork, &
            hiend, nloops, nmv, mpi_error_code)
         end if
         call dalloc (diag, 'DIAG', 'MANEIGmpi')
         call dalloc (iwork, 'IWORK', 'MANEIGmpi')
         call dalloc (jwork, 'JWORK', 'MANEIGmpi')

         if (mpi_rank .eq. 0) then
            write (*,301) nloops, nmv
            if (mpi_error_code .ne. 0) then
               write (*,302) mpi_error_code
            endif
         endif
      endif
!
!   Pick up the eigen pairs and store in EVAL and EVEC
!
      nend = csf_count*nvex
      do j = 1, eigenvector_counts_by_block(jblock)
         shifted_eigenvalues(ncminpat+j) = work(nend + iccmin(j + ncminpat))
         call dcopy (csf_count, work(csf_count*(iccmin(j + ncminpat) - 1) + 1:csf_count*(iccmin(j + ncminpat) - 1)+ &
            & (csf_count)), 1, &
      & csf_mixing_coefficients( nevecpat+csf_count*(j-1)+1:nevecpat+csf_count*(j-1)+(csf_count)), 1)
      end do
!     print *, ncminpat,(eval(ncminpat+j),j=1,nevblk(jblock)),
!    1 'zou,from maneig'
!
!   Deallocate storage
!
      call dalloc (work, 'WORK', 'MANEIGmpi')

  123 continue
      do jstate = 1, eigenvector_counts_by_block(jblock)
!
!   Find the dominant component of each eigenvector
!
         iofset = nevecpat + csf_count*(jstate - 1)

         amax = 0.d0
         do i = 1, csf_count
            wa = abs(csf_mixing_coefficients(i+iofset))
            if (wa <= amax) cycle
            amax = wa
            ia = i
         end do
!
!   Find the angular momentum and parity of the dominant component
!
!GG         IATJPO(JSTATE+NCMINPAT) = ITJPO(IA + NCFPAT)
!GG         IASPAR(JSTATE+NCMINPAT) = ISPAR(IA + NCFPAT)
!
!   Redefine eigenvectors so that the dominant component
!   is positive
!
         if (csf_mixing_coefficients(ia+iofset) >= 0.d0) cycle
         dnfac = -1.d0
         call dscal (csf_count, dnfac, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)
!===============================================================

      end do

  301 format('DVDSON: ',1i3,' loops; ',1i3,' matrix-vector multiplies.')
  302 format(' Returned from DVDSON with IERR = ',1i4)
  303 format(/,' ***** WARNING *****'/,/,&
         ' The angular momentum and parity of level ',1i2,' have changed:'/,&
         ' Last iteration: (2J+1) = ',1i2,', parity = ',1i2,';'/,&
         ' this iteration: (2J+1) = ',1i2,', parity = ',1i2,'.')

      return
      end subroutine maneigmpi
