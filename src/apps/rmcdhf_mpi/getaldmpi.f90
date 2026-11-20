!***********************************************************************
!
      subroutine getaldmpi
!
!   Interactively determines the data governing AL problem.
!
!   Quantities to be determined:
!     ncmin, wt(1:ncf), ucf(1:nw), nscf, nsic, orthst
!
!   Call(s) to: [LIB92]: ALLOC, IQ.
!
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992
!   Block version by Xinghong He          Last revision: 13 Jul 1998
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use array_memory_management, only: alloc

      use atomic_calculation_state, only: ncmin, nscf, nsic, wt
      use fixed_orbital_flags, only: fixed_orbital_count, is_orbital_fixed, max_orbitals
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count
      use scf_orthogonality_control, only: orthst
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getaldwt_i, only: getaldwt
      use csf_quantum_queries, only: iq
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: iqadum, orbital_index, csf_index
      real(real64) :: sum
      logical      :: getyn, yes
!-----------------------------------------------
!
      if (mpi_rank == 0) then
        write (stderr_unit, *) "(E)AL type calculation; H(DC) will not be ", &
                        "diagonalised;"
        write (stderr_unit, *) "getald ..."
        write (stderr_unit, *) "ncf=", csf_count
      end if

      call alloc (wt, csf_count, "WT", "GETALDmpi")

      if (mpi_rank == 0) call getaldwt (csf_count, wt)

      call mpi_bcast (wt(1), csf_count, mpi_double_precision, 0,        &
                      mpi_comm_world, mpi_error_code)

      do orbital_index = 1, orbital_count
         sum = 0.0_real64
         do csf_index = 1, csf_count
            sum = sum + wt(csf_index)*dble(iq(orbital_index,csf_index))
         end do
         ucf(orbital_index) = sum
      end do

      ncmin = 0
      nscf = 12
      nsic = 2 + (orbital_count - fixed_orbital_count)/4
      orthst = .false.

      return
      end subroutine getaldmpi
