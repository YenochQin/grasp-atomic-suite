!***********************************************************************
!                                                                      *
      subroutine getrmp
!                                                                      *
!   Interactively  determines the  list of radiation multipolarities   *
!   and parities. This is loadad into COMMON/OSC6/.                    *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:16:10   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use offd_c
      use osc_c, only: nkp, kp
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ndkp, istart, i, iend, lenth, ios, mult
      real, dimension(3) :: cnum
      logical :: lelec, lmagn, yes
      character :: record*256, reci
!-----------------------------------------------
!
!
!   Initial allocation for PNTRKP
!
      ndkp = 1
      call alloc (kp, ndkp, 'KP', 'GETRMP' )
      if (mpi_rank .eq. 0) then
!
!   Entry message
!
    1    continue
         write (6, *) 'Enter the list of transition specifications'
         write (6, *) ' e.g.,  E1,M2  or  E1 M2  or  E1;M2 :'
!
!   Initialise NKP
!
    2    continue
         read (*, '(A)') record
         nkp = 0
!
!   Parse RECORD from left to right
!
         istart = 0
         i = 1
    3    continue
         reci = record(i:i)
         if (reci/=' ' .and. reci/=',' .and. reci/=';') then
            if (istart == 0) istart = i
         else
            if (istart /= 0) then
               iend = i - 1
               reci = record(istart:istart)
               if (reci == 'E') then
                  lelec = .true.
                  lmagn = .false.
               else if (reci == 'M') then
                  lelec = .false.
                  lmagn = .true.
               else
                  write (6, *) 'GETRMP: Transitions must be of type'
                  write (6, *) ' E or type M; reenter ...'
                  go to 2
               endif
               lenth = iend - istart
               if (lenth /= 1) then
                  write (6, *) 'GETRMP: Transition multipolarities'
                  write (6, *) ' must be integers between 1 and 9;'
                  write (6, *) ' reenter ...'
                  go to 2
               endif
               reci = record(iend:iend)
               read (reci, '(1I1)', iostat=ios) mult
               if (ios /= 0) then
                  write (6, *) 'GETRMP: Unable to decode multipolarity'
                  write (6, *) ' '//reci//'; reenter ...'
                  go to 2
               endif
               nkp = nkp + 1
               if (nkp > ndkp) then
                  call ralloc (kp, nkp, 'KP', 'GETRMP' )
                  ndkp = nkp
               endif
               if (lelec) then
                  kp(nkp) = mult*(-1)**mult
               else if (lmagn) then
                  kp(nkp) = mult*(-1)**(mult + 1)
               endif
               istart = 0
            endif
         endif
!
         if (i < 256) then
            i = i + 1
            go to 3
         endif
!
         if (nkp == 0) go to 1
       endif !myid=0
       call mpi_bcast(nkp,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
!
!   Trim array to the exact size
!
      if (ndkp /= nkp) call ralloc (kp, nkp, 'KP', 'GETRMP')
      call mpi_bcast(kp(1),nkp,mpi_integer,0,mpi_comm_world, mpi_error_code)
!
!   If M1 or E2 inquire if the transitions are between levels
!   with different J quantum numbers.
!
      if(mpi_rank == 0) then
         do i = 1, nkp
            if (kp(i) == 1) then
               write (*, *)                                         &
               'M1 transitions only between levels with different J?'
               yes = getyn()
               if (yes) then
                  noffd1 = 1
               else
                  noffd1 = 0
               endif
            endif
            if (kp(i) /= 2) cycle
            write (*, *)                                            &
            'E2 transitions only between levels with different J?'
            yes = getyn()
            if (yes) then
               noffd2 = 1
            else
               noffd2 = 0
            endif
         end do
      endif !myid=0
      call mpi_bcast(noffd1,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast(noffd2,1,mpi_integer,0,mpi_comm_world, mpi_error_code)
!
      return
      end subroutine getrmp
