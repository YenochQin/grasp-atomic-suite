!***********************************************************************
!                                                                      *
      logical function triangbreit2 (ia,ib,ic,id,l)
!                                                                      *
!     Evaluate the triangle relations for Breit integrals of type 2    *
!     See Grant and McKenzier JPB 13 (1980), 2675-2681 formula 5       *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orb_c, only: nkl, nkj
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ia,ib,ic,id,l
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: la,lb,lc,ld,ja,jb,jc,jd
!-----------------------------------------------
!
      la = nkl(ia)
      lb = nkl(ib)
      lc = nkl(ic)
      ld = nkl(id)
      ja = nkj(ia)
      jb = nkj(ib)
      jc = nkj(ic)
      jd = nkj(id)

      triangbreit2 = .false.

      if ((mod(l+la+lb,2).eq.0).and.(abs(ja-jb).le.2*l).and.           &
         (abs(ja+jb).ge.2*l)) then
         if ((mod(l+lc+ld,2).eq.0).and.(abs(jc-jd).le.2*l).and.        &
            (abs(jc+jd).ge.2*l)) then
            triangbreit2 = .true.
            return
         end if
      end if

      return
      end function triangbreit2
