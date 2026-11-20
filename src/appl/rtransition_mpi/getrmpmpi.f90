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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use array_memory_management, only: alloc, ralloc
      use transition_off_diagonal_flags, only: noffd1, noffd2
      use transition_matrix_workspace, only: nkp, kp
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use, intrinsic :: iso_fortran_env, only: real32
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ndkp, istart, i, iend, lenth, ios, mult
      real(real32) :: cnum(3)
      logical :: lelec, lmagn, yes
      character :: record*256, reci
!-----------------------------------------------
!
!
!   Initial allocation for PNTRKP
!
      ndkp = 1
      call alloc (kp, ndkp, "KP", "GETRMP" )
      if (mpi_rank == 0) then
!
!   Entry message
!
    1    continue
         write (6, *) "Enter the list of transition specifications"
         write (6, *) " e.g.,  E1,M2  or  E1 M2  or  E1;M2 :"
!
!   Initialise NKP
!
    2    continue
         read (*, "(A)") record
         nkp = 0
!
!   Parse RECORD from left to right
!
         istart = 0
         i = 1
    3    continue
         reci = record(i:i)
         if (reci/=" " .and. reci/="," .and. reci/=";") then
            if (istart == 0) istart = i
         else
            if (istart /= 0) then
               iend = i - 1
               reci = record(istart:istart)
               if (reci == "E") then
                  lelec = .true.
                  lmagn = .false.
               else if (reci == "M") then
                  lelec = .false.
                  lmagn = .true.
               else
                  write (6, *) "GETRMP: Transitions must be of type"
                  write (6, *) " E or type M; reenter ..."
                  goto 2
               end if
               lenth = iend - istart
               if (lenth /= 1) then
                  write (6, *) "GETRMP: Transition multipolarities"
                  write (6, *) " must be integers between 1 and 9;"
                  write (6, *) " reenter ..."
                  goto 2
               end if
               reci = record(iend:iend)
               read (reci, "(1I1)", iostat=ios) mult
               if (ios /= 0) then
                  write (6, *) "GETRMP: Unable to decode multipolarity"
                  write (6, *) " "//reci//"; reenter ..."
                  goto 2
               end if
               nkp = nkp + 1
               if (nkp > ndkp) then
                  call ralloc (kp, nkp, "KP", "GETRMP" )
                  ndkp = nkp
               end if
               if (lelec) then
                  kp(nkp) = mult*(-1)**mult
               else if (lmagn) then
                  kp(nkp) = mult*(-1)**(mult + 1)
               end if
               istart = 0
            end if
         end if
!
         if (i < 256) then
            i = i + 1
            goto 3
         end if
!
         if (nkp == 0) goto 1
       end if !myid=0
       call mpi_bcast(nkp,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
!
!   Trim array to the exact size
!
      if (ndkp /= nkp) call ralloc (kp, nkp, "KP", "GETRMP")
      call mpi_bcast(kp(1),nkp,mpi_integer,0,mpi_comm_world, mpi_error_code)
!
!   If M1 or E2 inquire if the transitions are between levels
!   with different J quantum numbers.
!
      if(mpi_rank == 0) then
         do i = 1, nkp
            if (kp(i) == 1) then
               write (*, *)                                         &
               "M1 transitions only between levels with different J?"
               yes = getyn()
               if (yes) then
                  noffd1 = 1
               else
                  noffd1 = 0
               end if
            end if
            if (kp(i) /= 2) cycle
            write (*, *)                                            &
            "E2 transitions only between levels with different J?"
            yes = getyn()
            if (yes) then
               noffd2 = 1
            else
               noffd2 = 0
            end if
         end do
      end if !myid=0
      call mpi_bcast(noffd1,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast(noffd2,1,mpi_integer,0,mpi_comm_world, mpi_error_code)
!
      return
      end subroutine getrmp
