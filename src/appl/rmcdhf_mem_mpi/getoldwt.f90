!***********************************************************************
!
      subroutine getoldwt(ndef, ncmin, wt)
!
!   Interactively determines the weights for EOL calculation.
!   It's modified to always ask the question for the weight
!
!   Call(s) to: [LIB92]: GETYN.
!
!   Written by Xinghong He                Last revision: 19 Mar 1999
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use iounit_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ndef
      integer , intent(in) :: ncmin
      real(double)  :: wt(ncmin)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n159
      real(double) :: sumwgt
!-----------------------------------------------

! Standard weights: ncmin=1, OL calculation

      if (ncmin == 1) then
         wt(1) = -1.d0
         return
      endif

! Select a method to assign level weights for ncmin > 1 case

      write (istde, *) 'level weights (1 equal;  5 standard;  9 user)'

      ! Let user try 10 times to get the correct input. 10 is BIG
      ! enough since the idea here is to allow user mistakes and
      ! at the same time to avoid an infinity loop.

      do i = 1, 10
         read (istdi, *) n159
         if (n159==1 .or. n159==5 .or. n159==9) exit
         write (istde, *) 'Input not correct, do it again. tried=', i
      end do

      if (ndef.eq.0) then
        write(734,*) n159,'! level weights'
      end if

      if (i > 10) stop

!------------------------------------------------------------------

      select case (n159)                         ! Equal weight
      case (1)
         wt = -2.d0
      case (5)                                   ! Standard weight
         wt = -1.d0
      case (9)                                   ! User-input weight

  123    continue
         write (istde, *) 'Enter the (relative) weights of the', ncmin, &
            ' levels :'
         read (istdi, *) (wt(i),i=1,ncmin)

         sumwgt = 0.d0
         do i = 1, ncmin
            if (wt(i) <= 0.d0) then
               write (istde, *) 'Weights must exceed 0;'
               go to 123
            else
               sumwgt = sumwgt + wt(i)
            endif
         end do
         sumwgt = 1.d0/sumwgt
         wt = sumwgt*wt

      case default
         write (istde, *) 'Impossible ! Because it was guarded'
         stop
      end select

      return
      end subroutine getoldwt
