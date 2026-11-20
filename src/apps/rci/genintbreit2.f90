!***********************************************************************
!                                                                      *
      subroutine genintbreit2 (mpi_rank, mpi_process_count, nb, max_orbital_two_j)
!                                                                      *
!  Input:
!     myid, nprocs
!  Output:
!   N - Number of integrals
!   j2max - max of 2*J
!
!       Generate the list of Breit Integrals of type 2                 *
!       that could arise from a set of orbitals                        *
!       of orbitals.                                                   *
!                                                                      *
!       This routine is similar to genintrk                            *
!                                                                      *
!     Written by Per Jonsson            October 2014                   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use array_memory_management, only: alloc
      use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
          valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
      use orbital_csf_basis, only: orbital_two_j, orbital_count
      use breit_rank_offsets, only: kmax, kstartbreit1, kstartbreit2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use brintf_i, only: brintf
      use triangbreit2_i, only: triangbreit2
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mpi_rank, mpi_process_count
      integer             :: nb, max_orbital_two_j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: gen
      integer :: i, key, k, orbital_index, ib, ic, id
!-----------------------------------------------
!

      key = orbital_count + 1
!   Type 2 angular coefficients use rank nu+1 (nu >= 0).
!   The empty rank-zero interval and rank one start at the same offset.
      kstartbreit2(0:1) = 1
!
!   Find 2*JMAX; it should not be greater than PARAMETER KMAX
!
      max_orbital_two_j = orbital_two_j(1)
      do i = 2, orbital_count
         if (orbital_two_j(i) > max_orbital_two_j) max_orbital_two_j = orbital_two_j(i)
      end do

      if (max_orbital_two_j > kmax) then
         stop "genintbreit2: KMAX too small"
      end if
!
!   Make the breit integrals. Note that there are no symmetry relations that
!   can be utilized. See paper by Grasnt.
!   When GEN is false, sweep through to find dimension
!
      gen = .false.

  999 nb = 0
      do k = 1, max_orbital_two_j
         do orbital_index = 1, orbital_count
            do ib = 1, orbital_count
               do ic = 1, orbital_count
                  do id = 1, orbital_count
                     if (triangbreit2(orbital_index,ib,ic,id,k)) then
                        nb = nb + 1
                        if (gen .and.(mod(nb,mpi_process_count) == mpi_rank)) then
                           indtp2(nb) = ((orbital_index*key+ib)*key+ic)*key+id
                           valtp2(nb) = brintf (2,orbital_index,ib,ic,id,k)
                        end if
                     end if
                  end do
               end do
            end do
         end do
         kstartbreit2(k+1) = nb + 1
      end do
!
!     Allocate memory for integral book keeping
!
      if (.not. gen) then
         call alloc (indtp2,nb,"INDT2","genintbreit2")
         call alloc (valtp2,nb,"VALT2","genintbreit2")

! Initialization is necessary in the mpi version

         do i = 1, nb
            indtp2(i) = 0
            valtp2(i) = 0.0_real64
         end do

         if (mpi_rank == 0)                                        &
            print *, "Computing",nb," Breit integrals of type 2"

         gen = .true.
         goto 999
      end if

      return
      end subroutine genintbreit2
