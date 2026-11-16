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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double, long
      use hmat_c  !! rather than setham_to_genmat2
      use decide_c
      use eigv_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)        :: irestart
      integer(long), intent(out) :: nelmnt_a
      real(double), intent(in)   :: elsto
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(6) :: ntpi_a
      real(double)          :: densty
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
         write (24, 301) cutofftmp
         write (24, 302) ncoeitmp
         write (24, 303) ncoectmp
         write (24, 304) ncteitmp
         write (24, 305) nctectmp
         if (ltrans) then
            write (24, 306) ntpitmp
            write (24, 307) nmcbptmp
            write (24, 308) ncoretmp
         endif
         if (lvp) write (24, 309) nvpitmp
         if (lnms) write (24, 310) nkeitmp
         if (lsms) write (24, 311) nvintitmp
      else
         write (24, *) 'Restart mode --- no report on radial integrals'
      endif                                      !(irestart .NE. 0) ! non-restart mode


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

      write (24, 312) nelmnttmp
      write (24, 313) densty
      write (24, *)
      write (24, 300) hamiltonian_energy_shift
      write (24, *) hamiltonian_energy_shift, elsto, ncftmp

  300 format('Average energy = ',1p,d19.12,' Hartrees.')
  301 format('CUTOFF set to ',1p,d17.10)
  302 format('Dirac-Coulomb one-e radial integrals:',1i8)
  303 format('One-e angular integrals that exceed CUTOFF: ',1i8)
  304 format('Coulomb two-e radial integrals: ',1i8)
  305 format('Two-e angular integrals that exceed CUTOFF: ',1i8)
  306 format('Transverse two-e radial integrals: '/,6i8)
  307 format('MCBP coefficients that exceed CUTOFF: ',1i8)
  308 format('Core coefficients that exceed CUTOFF: ',1i8)
  309 format('Vacuum polarisation integrals: ',1i8)
  310 format('Kinetic energy integrals: ',1i8)
  311 format('Vinti integrals: ',1i8)
  312 format('Elements that exceed CUTOFF in the lower',&
         ' triangle of the H matrix: ',1i8)
  313 format('Density of the H(amiltonian) matrix: ',1p,d22.15)

      return
      end subroutine genmat2
