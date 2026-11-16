!***********************************************************************
!                                                                      *
      logical function triangbreit1 (ia,ib,ic,id,k)
!                                                                      *
!     Evaluate the triangular relation for the Breit integral of       *
!     See paper by I. Grant and B McKenzie, JPB 13, (1980) 2671-2681   *
!     formula 5                                                        *
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
      integer, intent(in) :: ia,ib,ic,id,k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: la,lb,lc,ld,ja,jb,jc,jd,l
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

      triangbreit1 = .false.

      do l = k-1,k+1
         if ((mod(k+la+lb,2).eq.1).and.(abs(ja-jb).le.2*l).and.     &
            (abs(ja+jb).ge.2*l)) then
            if ((mod(k+lc+ld,2).eq.1).and.(abs(jc-jd).le.2*l).and.  &
               (abs(jc+jd).ge.2*l)) then
               triangbreit1 = .true.
               return
            end if
         end if
      end do

      return
      end function triangbreit1
