!***********************************************************************
!                                                                      *
      subroutine setj(is, js, ks, ns, kj23)
!                                                                      *
!   Sets the tables required by  the recoupling  coefficient package   *
!
!
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:36   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use l1_c
      use l2_c
      use m_c,      only: jjc1, jjc2, jjq1, jjq2, jlist
      use couple_c, only: mja, nja, j1, j2, j3, mangm
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ns
      integer, intent(inout) :: kj23
      integer, dimension(2,2), intent(in) :: is
      integer, dimension(2,2), intent(in) :: js
      integer, dimension(2,2), intent(in) :: ks
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, ij, i, ns1, ns2, ns3, ij1, ij2, j
!-----------------------------------------------
!
!
!   1.0  Set J1 array
!
      ii = 0
      if (ns > 0) then
         j1(:ns) = jbq1(3,jlist(:ns))
         ii = ns
      endif
      if (ns /= 1) then
         ns1 = ns - 1
         if (ns1 > 0) then
            j1(ii+1:ns1+ii) = jjc1(:ns1)
            ii = ns1 + ii
         endif
         if (ns1 > 0) then
            j1(ii+1:ns1+ii) = jjc2(:ns1)
            ii = ns1 + ii
         endif
      endif
      do i = 1, 2
         ii = ii + 1
         ij = is(i,1)
         j1(ii) = jjq1(3,ij)
         if (i==1 .and. is(1,1)==is(2,1)) j1(ii) = jtq1(3)
         j1(ii+4) = ks(i,1)
      end do
      do i = 1, 2
         ii = ii + 1
         ij = is(i,2)
         j1(ii) = jjq2(3,ij)
         if (i==1 .and. is(1,2)==is(2,2)) j1(ii) = jtq2(3)
         j1(ii+4) = ks(i,2)
      end do
!
!   2.0  Set J2, J3 arrays if not already available
!
      ns2 = max(4,ns + 2)
      if (kj23 <= 0) then
!
         do i = 4, ns2
            j2(i,1) = ns + i - 4
            j2(i,2) = i - 2
            j2(i,3) = ns + i - 3
            j3(i,1) = j2(i,1) + ns - 1
            j3(i,2) = i - 2
            j3(i,3) = j2(i,3) + ns - 1
         end do
         j2(4,1) = 1
         j3(4,1) = 1
!
!   At this stage, the entries in rows corresponding to active
!   shells are set incorrectly.
!
!   3.0  Set rows 1 through 3
!
         ns3 = 3*ns
         j2(1,1) = ns3 + 5
         j2(1,2) = ns3 + 7
         j2(1,3) = ns3 + 3
         j2(2,1) = js(1,1)
         j2(2,2) = ns3 + 3
         j2(2,3) = ns3 - 1
         j2(3,1) = js(2,1)
         j2(3,2) = ns3 + 4
         j2(3,3) = ns3
!
         j3(1,1) = ns3 + 7
         j3(1,2) = ns3 + 4
         j3(1,3) = ns3 + 6
         j3(2,1) = js(1,2)
         j3(2,2) = ns3 + 5
         j3(2,3) = ns3 + 1
         j3(3,1) = js(2,2)
         j3(3,2) = ns3 + 6
         j3(3,3) = ns3 + 2
!
!   4.0  Set remaining resultants
!
         ij1 = js(1,1)
         ij2 = js(2,1)
         if (ij2 > 1) j2(ij2+2,2) = j2(3,3)
         if (ij2 == 1) j2(4,1) = j2(3,3)
         if (ij1 == ij2) then
            j2(3,1) = j2(2,3)
         else
!
            if (ij1 > 1) j2(ij1+2,2) = j2(2,3)
            if (ij1 == 1) j2(4,1) = j2(2,3)
         endif
!
         ij1 = js(1,2)
         ij2 = js(2,2)
         if (ij2 > 1) j3(ij2+2,2) = j3(3,3)
         if (ij2 == 1) j3(4,1) = j3(3,3)
         if (ij1 == ij2) then
            j3(3,1) = j3(2,3)
         else
!
            if (ij1 > 1) j3(ij1+2,2) = j3(2,3)
            if (ij1 == 1) j3(4,1) = j3(2,3)
         endif
!
!   All arrays now set. Put up flag KJ23.
!
         kj23 = 1
         mja = ns3 + 7
         nja = ns + 3
!
!   5.0  Save J2, J3 and return
!
         j2s(:ns2,:) = j2(:ns2,:)
         j3s(:ns2,:) = j3(:ns2,:)
         return
      endif
!
!   6.0  Reset J2, J3 from buffers if KJ23 has been set
!
      j2(:ns2,:) = j2s(:ns2,:)
      j3(:ns2,:) = j3s(:ns2,:)
      return
!
      end subroutine setj
