!***********************************************************************
!                                                                      *
      subroutine setqna(ja, jb)
!                                                                      *
!   This generates the  arrays  defining  the quantum numbers of the   *
!   states involved in the  matrix  element  linking  configurations   *
!   labelled by JA, JB.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: ICHOP, IQ, JCUP, JQS.                         *
!                                                                      *
!                                           Last update: 30 Oct 1987   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:39   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use m_c
      use orb_c,           only: ncf, nw
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iq_i
      use jqs_i
      use ichop_i
      use jcup_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ja
      integer  :: jb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, k, jcnt, jcntop, jw1, jw2, jw
!-----------------------------------------------
!
!
!   List parameters defining all shells in both configurations, whether
!   participating or not
!
      do j = 1, nw
         nq1(j) = iq(j,ja)
         nq2(j) = iq(j,jb)
         do k = 1, 3
            jjq1(k,j) = jqs(k,j,ja)
            jjq2(k,j) = jqs(k,j,jb)
         end do
      end do
!
!   Define coupling schemes: set JLIST array to define those shells
!   which are open in either configuration, and KLIST array to locate
!   the rest. Exclude shells which are empty in both configurations
!
      npeel = 0
      ncore = 0
      do j = 1, nw
         if (ichop(j,ja)==(-1) .and. ichop(j,jb)==(-1)) cycle
         if (ichop(j,ja)==1 .and. ichop(j,jb)==1) then
            ncore = ncore + 1
            klist(ncore) = j
         else
            npeel = npeel + 1
            jlist(npeel) = j
         endif
      end do
!
!   Return if not more than one shell is open
!
      if (npeel <= 1) return
!
!   Set arrays of coupling angular momenta interpolating closed
!   shells where necessary. Left hand side first ...
!
      jcnt = 1
      jcntop = 0
      jw1 = jlist(1)
      jw2 = jlist(2)
      if (ichop(jw1,ja) /= 0) then
         jjc1(1) = jqs(3,jw2,ja)
         if (ichop(jw2,ja) == 0) jcntop = 1
      else
         jcntop = 1
         if (ichop(jw2,ja) == 0) then
            jjc1(1) = jcup(jcnt,ja)
            jcnt = jcnt + 1
         else
            jjc1(1) = jqs(3,jw1,ja)
         endif
      endif
!
      do j = 3, npeel
         jw = jlist(j)
         if (ichop(jw,ja) /= 0) then
            jjc1(j-1) = jjc1(j-2)
         else
            if (jcntop /= 0) then
               jjc1(j-1) = jcup(jcnt,ja)
               jcnt = jcnt + 1
            else
               jjc1(j-1) = jqs(3,jw,ja)
            endif
            jcntop = jcntop + 1
         endif
      end do
!
!   ... and repeat for right hand side
!
      jcnt = 1
      jcntop = 0
      jw1 = jlist(1)
      jw2 = jlist(2)
      if (ichop(jw1,jb) /= 0) then
         jjc2(1) = jqs(3,jw2,jb)
         if (ichop(jw2,jb) == 0) jcntop = 1
      else
         jcntop = 1
         if (ichop(jw2,jb) == 0) then
            jjc2(1) = jcup(jcnt,jb)
            jcnt = jcnt + 1
         else
            jjc2(1) = jqs(3,jw1,jb)
         endif
      endif
!
      do j = 3, npeel
         jw = jlist(j)
         if (ichop(jw,jb) /= 0) then
            jjc2(j-1) = jjc2(j-2)
         else
            if (jcntop /= 0) then
               jjc2(j-1) = jcup(jcnt,jb)
               jcnt = jcnt + 1
            else
               jjc2(j-1) = jqs(3,jw,jb)
            endif
            jcntop = jcntop + 1
         endif
      end do
!
      return
      end subroutine setqna
