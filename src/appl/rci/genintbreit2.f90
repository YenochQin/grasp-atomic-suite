!***********************************************************************
!                                                                      *
      subroutine genintbreit2 (myid, nprocs, nb, j2max)
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
      use memory_man
      use bilst_c
      use orb_c, only: nkj, nw
      use kkstartbreit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use brintf_i
      use triangbreit2_i
      implicit none
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

      key = nw + 1
      kstartbreit2(0) = 1
!
!   Find 2*JMAX; it should not be greater than PARAMETER KMAX
!
      j2max = nkj(1)
      do i = 2, nw
         if (nkj(i) .gt. j2max) j2max = nkj(i)
      enddo

      if (j2max .gt. kmax) then
         stop 'genintbreit2: KMAX too small'
      endif
!
!   Make the breit integrals. Note that there are no symmetry relations that
!   can be utilized. See paper by Grasnt.
!   When GEN is false, sweep through to find dimension
!
      gen = .false.

  999 nb = 0
      do k = 0, j2max
         do ia = 1, nw
            do ib = 1, nw
               do ic = 1, nw
                  do id = 1, nw
                     if (triangbreit2(ia,ib,ic,id,k)) then
                        nb = nb + 1
                        if (gen .and.(mod(nb,nprocs) .eq. myid)) then
                           indtp2(nb) = ((ia*key+ib)*key+ic)*key+id
                           valtp2(nb) = brintf (2,ia,ib,ic,id,k)
                        endif
                     endif
                  enddo
               enddo
            enddo
         enddo
         kstartbreit2(k+1) = nb + 1
      enddo
!
!     Allocate memory for integral book keeping
!
      if (.not. gen) then
         call alloc (indtp2,nb,'INDT2','genintbreit2')
         call alloc (valtp2,nb,'VALT2','genintbreit2')

! Initialization is necessary in the mpi version

         do i = 1, nb
            indtp2(i) = 0
            valtp2(i) = 0.d0
         enddo

         if (myid .eq. 0)                                        &
            print *, 'Computing',nb,' Breit integrals of type 2'

         gen = .true.
         goto 999
      endif

      return
      end subroutine genintbreit2
