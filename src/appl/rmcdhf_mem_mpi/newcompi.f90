!***********************************************************************
!                                                                      *
      subroutine newcompi(sum)
!                                                                      *
!   This routine computes the level weights, the generalized occupa-   *
!   tion numbers,  and average energy for  EOL calculations;  this     *
!   information and the eigenvectors are then printed out.             *
!                                                                      *
!   Call(s) to: [LIB92]: IQ.                                           *
!               [RSCF92]: CSFWGT, DSUBRS.                              *
!                                         Last revision: 24 Dec 1992   *
!   Block version by Xinghong He          Last revision: 05 Aug 1998   *
!                                                                      *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:53:04   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use debug_c
      use def_c
      use eigv_c
      use orb_c
      use scf_c
      use syma_c,          only: jpgg
      use iounit_c
      use hblock_c
      use pos_c
      use peav_c
      use blkidx_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dsubrs_i
      use iq_i
      use csfwgt_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(out) :: sum
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, nb, i, noff, jall, jblock
      real(double) :: weitj, ee
      logical :: eol
!-----------------------------------------------
!   Compute weighting factors
!
      sum = 0.d0
      do j = 1, ncmin
         weitj = weight(j)
         if (weitj == (-2.d0)) then
            wt(j) = 1.d0
         else if (weitj == (-1.d0)) then
!GG            WT(J) = IATJPO(J)
            jblock = idxblk(j)     ! Block number of this state
            wt(j) = abs(jpgg(jblock))
         else
            wt(j) = weight(j)
         endif
         sum = sum + wt(j)
      end do

      wt(:ncmin) = wt(:ncmin)/sum
!
!   Compute generalised occupation numbers
!   <----- Distributed ----->
!
      eol = .true.

      do j = 1, orbital_count
         sum = 0.d0
         do nb = 1, symmetry_block_count
            do i = mpi_rank + 1, csf_counts_by_block(nb), mpi_process_count
               sum = sum + dsubrs(eol,i,i,nb)*iq(j,i + csf_offsets_by_block(nb))
            end do
         end do
!         UCF(J) = SUM
         call mpi_allreduce (sum, ucf(j), 1, mpi_double_precision, &
                             mpi_sum, mpi_comm_world, mpi_error_code)
      end do
!
!   Write out level energies and weights
!
      if(mpi_rank == 0) write (*, 300)
      sum = 0.d0
      noff = 0
      do jall = 1, ncmin
         nb = idxblk(jall)                       ! Block number of this state
         ee = eavblk(nb) + shifted_eigenvalues(jall)
         if (mpi_rank .eq. 0) then
            write (*, 301) iccmin(jall), ee, wt(jall)
            if (ldbpg(5)) then
               write (99, *) jall, nb, csf_counts_by_block(nb), mixing_offsets_by_block(nb)
               write (99, 302)
               write (99, 303) (csf_mixing_coefficients(i + noff),i=1,csf_counts_by_block(nb))
               noff = noff + csf_counts_by_block(nb)
            endif
         endif
         sum = sum + wt(jall)*ee
      end do
            !~~~~~~~~~
            if (mpi_rank .eq. 0) then
            !~~~~~~~~~

      call csfwgt (.true.)
!
!   Write out average energy
!
      if (ncmin > 1) write (*, 304) sum
      write (0,'(A25,1PD18.10)') 'Weighted average energy: ', sum
!
!   Write out generalized occupation numbers
!
      write (*, 305)
      write (*, 303) (ucf(i),i=1,orbital_count)

            !~~~~~~~~~
            endif
            !~~~~~~~~~

  300 format(/,'Optimise on the following level(s):'/)
  301 format('Level ',1i2,4x,'Energy = ',1p,1d19.12,4x,'Weight = ',1d12.5)
  302 format(/,'Configuration mixing coefficients:')
  303 format(1x,1p,6d12.4)
  304 format(/,'Weighted average energy of these levels = ',1p,d18.10)
  305 format(/,'Generalised occupation numbers:'/)

      return
      end subroutine newcompi
