!***********************************************************************
!                                                                      *
      subroutine genintrk (myid, nprocs, n, j2max)
!                                                                      *
!  Input:
!     myid, nprocs
!  Output:
!   N - Number of integrals
!   j2max - max of 2*J
!
!       Generate the list of Rk Integrals that could arise from a set  *
!       of orbitals.                                                   *
!                                                                      *
!     Written by Per Jonsson                                           *
!   MPI version by Xinghong He            Last revision: 22 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use array_memory_management, only: alloc
      use distributed_slater_integrals, only: indteirk, valteirk, rk_local_count
      use coulomb_rank_offsets, only: kmax, kstart
      use orbital_csf_basis, only: orbital_count, orbital_l, orbital_two_j
      use angular_symmetry_labels, only: jlbl, jlbp, jlbr
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use slater_radial_integral, only: slater
      use triangrk_i, only: triangrk
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: myid, nprocs
      integer, intent(out) :: n, j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: gen
      integer :: key, i, k, ia, ib, ic, id
!-----------------------------------------------------------------------

      key = orbital_count + 1
      kstart(0) = 1
!
!   Find 2*JMAX; it should not be greater than PARAMETER KMAX
!
      j2max = orbital_two_j(1)
      do i = 2, orbital_count
         if (orbital_two_j(i) > j2max) j2max = orbital_two_j(i)
      end do

      if (j2max > kmax) then
         stop "genintrk: KMAX too small"
      end if
!
!   Make the RK integrals: IA <= IB, IA <= IC, IA <= ID, IB <= ID
!   When GEN is false, sweep through to find dimension
!
      gen = .false.

  999 n = 0
      do k = 0, j2max
         do ia = 1, orbital_count
            do ib = ia, orbital_count
               do ic = ia, orbital_count
                  if (triangrk(orbital_l(ia),k,orbital_l(ic))) then
                     do id = ib, orbital_count
                        if (triangrk(orbital_l(ib),k,orbital_l(id))) then
                           n = n + 1
                           if (gen .and. (mod(n,nprocs) == myid)) then
                              indteirk(n) = ((ia*key+ib)*key+ic)*key+id
                              valteirk(n) = slater (ia,ib,ic,id,k)
                           end if
                        end if
                     end do
                  end if
               end do
            end do
         end do
         kstart(k+1) = n + 1
      end do
!
!     Allocate memory for integral book keeping
!
      if (.not. gen) then
         call alloc (indteirk,n,"INDTEIRK", "GENINTRK")
         call alloc (valteirk,n,"VALTEIRK", "GENINTRK")

! Initialization is necessary in the mpi version

         do i = 1, n
            indteirk(i) = 0
            valteirk(i) = 0.0_real64
         end do

         if (myid == 0)                                               &
     &      print *, "Allocating space for ",n," Rk integrals"

         gen = .true.
         goto 999
      end if

      return
      end subroutine genintrk
