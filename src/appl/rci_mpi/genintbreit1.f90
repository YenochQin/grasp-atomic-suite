!***********************************************************************
!                                                                      *
      subroutine genintbreit1 (myid, nprocs, nb, j2max)
!                                                                      *
!  Input:
!     myid, nprocs
!  Output:
!   N - Number of integrals
!   j2max - max of 2*J
!
!       Generate the list of Breit Integrals of type 1                 *
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
      use memory_man
      use bilst_c
      use orb_c, only: orbital_two_j, orbital_count
      use kkstartbreit_c, only: kmax, kstartbreit1, kstartbreit2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use brintf_i, only: brintf
      use triangbreit1_i, only: triangbreit1
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: myid, nprocs
      integer             :: nb, j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: gen
      integer :: i, key, k, ia, ib, ic, id
!-----------------------------------------------
!
      key = orbital_count + 1
      kstartbreit1(0) = 1
!
!   Find 2*JMAX; it should not be greater than PARAMETER KMAX
!
      j2max = orbital_two_j(1)
      do i = 2, orbital_count
         if (orbital_two_j(i) > j2max) j2max = orbital_two_j(i)
      end do

      if (j2max > kmax) then
         stop "genintbreit1: KMAX too small"
      end if
!
!   Make the breit integrals. Note that there are no symmetry relations that
!   can be utilized. See paper by Grasnt.
!   When GEN is false, sweep through to find dimension
!
      gen = .false.

  999 nb = 0
      do k = 0, j2max
         do ia = 1, orbital_count
            do ib = 1, orbital_count
               do ic = 1, orbital_count
                  do id = 1, orbital_count
                     if (triangbreit1(ia,ib,ic,id,k)) then
                        nb = nb + 1
                        if (gen .and.(mod(nb,nprocs) == myid)) then
                           indtp1(nb) = ((ia*key+ib)*key+ic)*key+id
                           valtp1(nb) = brintf (1,ia,ib,ic,id,k)
                        end if
                     end if
                  end do
               end do
            end do
         end do
         kstartbreit1(k+1) = nb + 1
      end do
!
!     Allocate memory for integral book keeping
!
      if (.not. gen) then
         call alloc (indtp1,nb,"INDT1","genintbreit1")
         call alloc (valtp1,nb,"VALT1","genintbreit1")

! Initialization is necessary in the mpi version

         do i = 1, nb
            indtp1(i) = 0
            valtp1(i) = 0.0_real64
         end do

         if (myid == 0)                                         &
            print *, "Computing",nb," Breit integrals of type 1"

         gen = .true.
         goto 999
      end if

      return
      end subroutine genintbreit1
