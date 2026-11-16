!***********************************************************************
!
      subroutine getaldwt(ncf, wt)
!
!   Interactively determines the weights.
!
!   Call(s) to: [LIB92]: ITJPO
!
!   Written by Xinghong He                Last revision: 19 Mar 1999
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use itjpo_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ncf
      real(real64)  :: wt(ncf)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n159
      real(real64) :: sumwgt, ftjpoi
!-----------------------------------------------

! Select a method to assign level weights for NCF > 1 case

      write (stderr_unit, *) 'level weights (1 equal;  5 standard;  9 user)'

      ! Let user try 10 times to get the correct input. 10 is BIG
      ! enough since the idea here is to allow user mistakes and
      ! at the same time to avoid an infinity loop.

      do i = 1, 10
         read (stdin_unit, *) n159
         if (n159==1 .or. n159==5 .or. n159==9) exit
         write (stderr_unit, *) 'Input not correct, do it again. tried=', i
      end do

      if (i > 10) stop 'Must be running un-attended'

!------------------------------------------------------------------

      select case (n159)                         ! Equal weight
      case (1)
         wt = 1._real64
         sumwgt = dble(ncf)
      case (5)                                   ! Standard weight
         sumwgt = 0._real64
         do i = 1, ncf
            ftjpoi = dble(itjpo(i))
            wt(i) = ftjpoi
            sumwgt = sumwgt + ftjpoi
         end do
      case (9)                                   ! User-input weight

  123    continue
         write (stderr_unit, *) 'Enter the (relative) weights of the', ncf, &
            ' levels :'
         read (stdin_unit, *) (wt(i),i=1,ncf)

         sumwgt = 0._real64
         do i = 1, ncf
            if (wt(i) <= 0._real64) then
               write (stderr_unit, *) 'Weights must exceed 0;'
               go to 123
            else
               sumwgt = sumwgt + wt(i)
            endif
         end do

      case default
         write (stderr_unit, *) 'Impossible ! Because it was guarded'
         stop
      end select

      sumwgt = 1._real64/sumwgt
      wt = sumwgt*wt

      return
      end subroutine getaldwt
