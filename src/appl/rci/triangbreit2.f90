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
      use orb_c, only: orbital_l, orbital_two_j
      implicit none (type, external)
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
      la = orbital_l(ia)
      lb = orbital_l(ib)
      lc = orbital_l(ic)
      ld = orbital_l(id)
      ja = orbital_two_j(ia)
      jb = orbital_two_j(ib)
      jc = orbital_two_j(ic)
      jd = orbital_two_j(id)

      triangbreit2 = .false.

      if ((mod(l+la+lb,2)==0).and.(abs(ja-jb)<=2*l).and.           &
         (abs(ja+jb)>=2*l)) then
         if ((mod(l+lc+ld,2)==0).and.(abs(jc-jd)<=2*l).and.        &
            (abs(jc+jd)>=2*l)) then
            triangbreit2 = .true.
            return
         end if
      end if

      return
      end function triangbreit2
