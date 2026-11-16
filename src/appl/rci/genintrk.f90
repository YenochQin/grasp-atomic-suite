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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use cteilsrk_c
      use kkstart_c
      use orb_c
      use jlabl_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use slater_i
      use triangrk_i, only: triangrk
      implicit none
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

      key = nw + 1
      kstart(0) = 1
!
!   Find 2*JMAX; it should not be greater than PARAMETER KMAX
!
      j2max = nkj(1)
      do i = 2, nw
         if (nkj(i) .gt. j2max) j2max = nkj(i)
      enddo

      if (j2max .gt. kmax) then
         stop 'genintrk: KMAX too small'
      endif
!
!   Make the RK integrals: IA <= IB, IA <= IC, IA <= ID, IB <= ID
!   When GEN is false, sweep through to find dimension
!
      gen = .false.

  999 n = 0
      do k = 0, j2max
         do ia = 1, nw
            do ib = ia, nw
               do ic = ia, nw
                  if (triangrk(nkl(ia),k,nkl(ic))) then
                     do id = ib, nw
                        if (triangrk(nkl(ib),k,nkl(id))) then
                           n = n + 1
                           if (gen .and. (mod(n,nprocs) .eq. myid)) then
                              indteirk(n) = ((ia*key+ib)*key+ic)*key+id
                              valteirk(n) = slater (ia,ib,ic,id,k)
                           endif
                        endif
                     enddo
                  endif
               enddo
            enddo
         enddo
         kstart(k+1) = n + 1
      enddo
!
!     Allocate memory for integral book keeping
!
      if (.not. gen) then
         call alloc (indteirk,n,'INDTEIRK', 'GENINTRK')
         call alloc (valteirk,n,'VALTEIRK', 'GENINTRK')

! Initialization is necessary in the mpi version

         do i = 1, n
            indteirk(i) = 0
            valteirk(i) = 0.d0
         enddo

         if (myid .eq. 0)                                               &
     &      print *, 'Allocating space for ',n,' Rk integrals'

         gen = .true.
         goto 999
      endif

      return
      end
