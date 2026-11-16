!***********************************************************************
!                                                                      *
    real(kind(0.0d0))  function vpintf (ia,ib)
!                                                                      *
!   Computes nuclear vacuum polarization integrals.                    *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
     use vast_kind_param, only: double
     use parameter_def,   only: keyorb
     use debug_c
     use grid_c
     use ncdist_c
     use orb_c
     use tatb_c, only: mtp, ta, tb
     use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
     use quad_i
     implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ia, ib
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer :: k
!------------------------------------------------
!
      mtp = min (mf(ia),mf(ib))
      !TA(1) = 0.0D 00  !! Why would this be wrong?!!
      ta(1)  = 0
      do 1 k = 2,mtp
         ta(k) = (pf(k,ia)*pf(k,ib)+qf(k,ia)*qf(k,ib))                  &
                 *zdist(k)
    1 continue
      call quad (vpintf)
!
      if (ldbpr(9)) write (99,300) np(ia),nh(ia),np(ib),nh(ib),vpintf
!
      return
!
  300 format (/'VPINTF: V (',1i2,1a2,',',1i2,1a2,') = ',1pd19.12)
!
      end function vpintf
