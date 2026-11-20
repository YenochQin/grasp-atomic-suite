!***********************************************************************
!                                                                      *
      subroutine genmat2(irestart, nelmnt_a, elsto)
!
!   Get eav and do writings to the summary file .csum
!   The mpi version (genmat2mpi) also gets nelmnt_a and elsto
!
! Xinghong He 98-06-15
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use sparse_hamiltonian_storage, only: cutofftmp, ncftmp, ncoectmp, ncoeitmp, ncoretmp, nctectmp, ncteitmp, nelmnttmp, &
         & nkeitmp, nmcbptmp, ntpitmp, nvintitmp, nvpitmp  !! rather than setham_to_genmat2
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)        :: irestart
      integer(int64), intent(out) :: nelmnt_a
      real(real64), intent(in)   :: elsto
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntpi_a(6)
      real(real64)          :: densty
!-----------------------------------------------
!-----------------------------------------------------------------------
! Transfer parameters around and print out on node-0.
!  .Temporary variables (ending with "_a") are used to store the
!   sums over nodes;
!  .Outputs to stream 24 are done on node-0;
! The following parameters are accumulated in setham which will not
! contain the correct values in restart mode. And thus is skipped.
!-----------------------------------------------------------------------
      if (irestart /= 0) then                    ! non-restart mode
         write (24, '("CUTOFF set to ",1p,d17.10)') cutofftmp
         write (24, '("Dirac-Coulomb one-e radial integrals:",1i8)') ncoeitmp
         write (24, '("One-e angular integrals that exceed CUTOFF: ",1i8)') ncoectmp
         write (24, '("Coulomb two-e radial integrals: ",1i8)') ncteitmp
         write (24, '("Two-e angular integrals that exceed CUTOFF: ",1i8)') nctectmp
         if (ltrans) then
            write (24, '("Transverse two-e radial integrals: "/,6i8)') ntpitmp
            write (24, '("MCBP coefficients that exceed CUTOFF: ",1i8)') nmcbptmp
            write (24, '("Core coefficients that exceed CUTOFF: ",1i8)') ncoretmp
         end if
         if (lvp) write (24, '("Vacuum polarisation integrals: ",1i8)') nvpitmp
         if (lnms) write (24, '("Kinetic energy integrals: ",1i8)') nkeitmp
         if (lsms) write (24, '("Vinti integrals: ",1i8)') nvintitmp
      else
         write (24, *) "Restart mode --- no report on radial integrals"
      end if                                      !(irestart .NE. 0) ! non-restart mode


!-----------------------------------------------------------------------
! ELSTO, EAV are not only for print-out, but also used later.
! density of the Hamiltonian matrix is only for print-out.
! At this place Hamiltonian matrix elements (EMT on each node) do
!   _not_ contain ELSTO. ELSTO will be added to the total energy
!  later with EAV.
!-----------------------------------------------------------------------

      nelmnt_a = nelmnttmp

      densty = dble(nelmnttmp)/dble((ncftmp*(ncftmp + 1))/2)

      hamiltonian_energy_shift = hamiltonian_energy_shift/dble(ncftmp) + elsto

      write (24, '("Elements that exceed CUTOFF in the lower"," triangle of the H matrix: ",1i8)') nelmnttmp
      write (24, '("Density of the H(amiltonian) matrix: ",1p,d22.15)') densty
      write (24, *)
      write (24, '("Average energy = ",1p,d19.12," Hartrees.")') hamiltonian_energy_shift
      write (24, *) hamiltonian_energy_shift, elsto, ncftmp


      return
      end subroutine genmat2
