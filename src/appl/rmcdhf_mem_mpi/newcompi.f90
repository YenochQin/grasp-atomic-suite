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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: iccmin, ncmin, weight, wt
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use orbital_csf_basis, only: max_orbitals, orbital_count
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use level_symmetry_metadata,          only: jpgg
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
          csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
      use block_mean_energies, only: eavblk
      use transition_block_indices, only: idxblk
      use mpi_state, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_error_code, &
         & mpi_process_count, mpi_rank, mpi_sum
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dsubrs_i, only: dsubrs
      use csf_quantum_queries, only: iq
      use csfwgt_i, only: csfwgt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64) , intent(out) :: sum
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, nb, i, noff, jall, jblock
      real(real64) :: weitj, ee
      logical :: eol
!-----------------------------------------------
!   Compute weighting factors
!
      sum = 0.0_real64
      do j = 1, ncmin
         weitj = weight(j)
         if (weitj == (-2.0_real64)) then
            wt(j) = 1.0_real64
         else if (weitj == (-1.0_real64)) then
!GG            WT(J) = IATJPO(J)
            jblock = idxblk(j)     ! Block number of this state
            wt(j) = abs(jpgg(jblock))
         else
            wt(j) = weight(j)
         end if
         sum = sum + wt(j)
      end do

      wt(:ncmin) = wt(:ncmin)/sum
!
!   Compute generalised occupation numbers
!   <----- Distributed ----->
!
      eol = .true.

      do j = 1, orbital_count
         sum = 0.0_real64
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
      if(mpi_rank == 0) write (*, '(/,"Optimise on the following level(s):"/)')
      sum = 0.0_real64
      noff = 0
      do jall = 1, ncmin
         nb = idxblk(jall)                       ! Block number of this state
         ee = eavblk(nb) + shifted_eigenvalues(jall)
         if (mpi_rank == 0) then
            write (*, '("Level ",1i2,4x,"Energy = ",1p,1d19.12,4x,"Weight = ",1d12.5)') iccmin(jall), ee, wt(jall)
            if (ldbpg(5)) then
               write (99, *) jall, nb, csf_counts_by_block(nb), mixing_offsets_by_block(nb)
               write (99, '(/,"Configuration mixing coefficients:")')
               write (99, '(1x,1p,6d12.4)') (csf_mixing_coefficients(i + noff),i=1,csf_counts_by_block(nb))
               noff = noff + csf_counts_by_block(nb)
            end if
         end if
         sum = sum + wt(jall)*ee
      end do
            !~~~~~~~~~
            if (mpi_rank == 0) then
            !~~~~~~~~~

      call csfwgt (.true.)
!
!   Write out average energy
!
      if (ncmin > 1) write (*, '(/,"Weighted average energy of these levels = ",1p,d18.10)') sum
      write (0,"(A25,1PD18.10)") "Weighted average energy: ", sum
!
!   Write out generalized occupation numbers
!
      write (*, '(/,"Generalised occupation numbers:"/)')
      write (*, '(1x,1p,6d12.4)') (ucf(i),i=1,orbital_count)

            !~~~~~~~~~
            end if
            !~~~~~~~~~


      return
      end subroutine newcompi
