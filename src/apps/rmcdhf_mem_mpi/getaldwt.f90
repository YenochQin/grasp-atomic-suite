!***********************************************************************
!
      subroutine getaldwt(csf_count, level_weights)
!
!   Interactively determines the weights.
!
!   Call(s) to: [LIB92]: ITJPO
!
!   Written by Xinghong He                Last revision: 19 Mar 1999
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use itjpo_i, only: itjpo
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: csf_count
      real(real64)  :: level_weights(csf_count)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n159
      real(real64) :: sumwgt, ftjpoi
!-----------------------------------------------

! Select a method to assign level weights for NCF > 1 case

      write (stderr_unit, *) "level weights (1 equal;  5 standard;  9 user)"

      ! Let user try 10 times to get the correct input. 10 is BIG
      ! enough since the idea here is to allow user mistakes and
      ! at the same time to avoid an infinity loop.

      do i = 1, 10
         read (stdin_unit, *) n159
         if (n159==1 .or. n159==5 .or. n159==9) exit
         write (stderr_unit, *) "Input not correct, do it again. tried=", i
      end do

      if (i > 10) stop "Must be running un-attended"

!------------------------------------------------------------------

      select case (n159)                         ! Equal weight
      case (1)
         level_weights = 1.0_real64
         sumwgt = dble(csf_count)
      case (5)                                   ! Standard weight
         sumwgt = 0.0_real64
         do i = 1, csf_count
            ftjpoi = dble(itjpo(i))
            level_weights(i) = ftjpoi
            sumwgt = sumwgt + ftjpoi
         end do
      case (9)                                   ! User-input weight

  123    continue
         write (stderr_unit, *) "Enter the (relative) weights of the", csf_count, &
            " levels :"
         read (stdin_unit, *) (level_weights(i),i=1,csf_count)

         sumwgt = 0.0_real64
         do i = 1, csf_count
            if (level_weights(i) <= 0.0_real64) then
               write (stderr_unit, *) "Weights must exceed 0;"
               goto 123
            else
               sumwgt = sumwgt + level_weights(i)
            end if
         end do

      case default
         write (stderr_unit, *) "Impossible ! Because it was guarded"
         stop
      end select

      sumwgt = 1.0_real64/sumwgt
      level_weights = sumwgt*level_weights

      return
      end subroutine getaldwt
