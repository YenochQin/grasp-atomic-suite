!***********************************************************************
!
      subroutine getoldwt(ndef, ncmin, level_weights)
!
!   Interactively determines the weights for EOL calculation.
!   It's modified to always ask the question for the weight
!
!   Call(s) to: [LIB92]: GETYN.
!
!   Written by Xinghong He                Last revision: 19 Mar 1999
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ndef
      integer , intent(in) :: ncmin
      real(real64)  :: level_weights(ncmin)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n159
      real(real64) :: sumwgt
!-----------------------------------------------

! Standard weights: ncmin=1, OL calculation

      if (ncmin == 1) then
         level_weights(1) = -1.0_real64
         return
      end if

! Select a method to assign level weights for ncmin > 1 case

      write (stderr_unit, *) "level weights (1 equal;  5 standard;  9 user)"

      ! Let user try 10 times to get the correct input. 10 is BIG
      ! enough since the idea here is to allow user mistakes and
      ! at the same time to avoid an infinity loop.

      do i = 1, 10
         read (stdin_unit, *) n159
         if (n159==1 .or. n159==5 .or. n159==9) exit
         write (stderr_unit, *) "Input not correct, do it again. tried=", i
      end do

      if (ndef==0) then
        write(734,*) n159,"! level weights"
      end if

      if (i > 10) stop

!------------------------------------------------------------------

      select case (n159)                         ! Equal weight
      case (1)
         level_weights = -2.0_real64
      case (5)                                   ! Standard weight
         level_weights = -1.0_real64
      case (9)                                   ! User-input weight

  123    continue
         write (stderr_unit, *) "Enter the (relative) weights of the", ncmin, &
            " levels :"
         read (stdin_unit, *) (level_weights(i),i=1,ncmin)

         sumwgt = 0.0_real64
         do i = 1, ncmin
            if (level_weights(i) <= 0.0_real64) then
               write (stderr_unit, *) "Weights must exceed 0;"
               goto 123
            else
               sumwgt = sumwgt + level_weights(i)
            end if
         end do

         if (ndef == 0) then
            write(734,*) (level_weights(i), i = 1, ncmin)
         end if

         sumwgt = 1.0_real64/sumwgt
         level_weights = sumwgt*level_weights

      case default
         write (stderr_unit, *) "Impossible ! Because it was guarded"
         stop
      end select

      return
      end subroutine getoldwt
