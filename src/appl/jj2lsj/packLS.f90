!
!***********************************************************************
!                                                                      *
      subroutine packls(m, el, q, couple, str)
!                                                                      *
!   Subroutine written by Bin LIU                                      *
!                                                                      *
!   Rules for encoding                                                 *
!   1. All blanks deleted                                              *
!   2. If Qi=1, omit Qi                                                *
!   3. If Qi=1 or Qi>=4l+1, omit ALFAi                                 *
!   4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi.          *
!                                                                      *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01
!...Switches:
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lval_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)            :: m
!GG      CHARACTER(LEN=64), INTENT(OUT) :: STR
      character(len=164), intent(out) :: str
      integer, intent(in)            :: q(*)
      character(len=3), intent(in)   :: el(*), couple(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer          :: full1, full, const, j, i, k, n, n1, nn, icore
      integer          :: l, previous_q
      character(len=1) :: ch1
      character(len=3) :: ch3
!-----------------------------------------------
!  FULL   :  4l+2
!  CONST  :   constant for converting lowercase to uppercase
!  CH*    :  temporary variables
!
      const = ichar('a') - ichar('A')
      str = ' '
      j = 1
      icore = 0
      full = 0; n = 0; full1 = 0; n1 = 0
!
!   -----  begin to encode  -----
!
      do i = 1, m
!        Fortran does not guarantee short-circuit evaluation of .AND.;
!        never reference Q(0), even when a condition starts with I /= 1.
         previous_q = 0
         if (i > 1) previous_q = q(i-1)
         if (q(i) == 0) cycle
         if (el(i)(1:1) == ' ') then
            str(j:j+1) = el(i)(2:3)
            k = 2
         else
            str(j:j+2) = el(i)
            k = 3
            if (el(i)(3:3) == ' ') k = 2
         endif
         ch1 = str(j+1:j+1)
         if (ch1>='A' .and. ch1<='Z') str(j+1:j+1) = char(ichar(ch1) + const)
         if(i /= 1) full1 = full
         if(i /= 1) n1 = n
         l = lval(ch1)
         full = 4*l + 2
!
!  -----  convert Qi into character  -----
!
         j = j + k
         n = q(i)
         if (n > 14) stop 'Too many electrons'
         if(i /= 1 .and. full1 /= n1) icore = icore + full1 - n1
!
!  -----  If Qi<>1, add Qi
!           If Qi<4l+1, add TERMi for the shell -----
!
         if (n /= 1) then
            str(j:j) = '('
            j = j + 1
            nn = n
            if (n > 9) then
               str(j:j) = '1'
               j = j + 1
               nn = n - 10
            endif
            str(j:j+1) = char(ichar('0') + nn)//')'
            j = j + 2
            if(i/=1 .and. icore<2 .and. previous_q<=1 .and. n==2) then
               continue
            else if(i/=1 .and. icore==0 .and. previous_q==full1 .and. n<3  &
                                                      .and. m==i) then
               continue
            else if(i/=1 .and. icore==0 .and. previous_q==full1 .and. l<2  &
                                                      .and. m==i) then
               continue
            else if(i/=1 .and. icore==0 .and. previous_q==full1 .and. l==2 &
                                          .and. n > 7 .and. m==i) then
               continue
            else if(i/=1 .and. icore==0 .and. previous_q==full1 .and. l==3 &
                                         .and. n > 11 .and. m==i) then
               continue
            else if (n<(full-1) .and. m/=1) then
               ch3 = couple(i)
               ch1 = ch3(2:2)
               if (ch1>='a' .and. ch1<='z') ch3(2:2) = char(ichar(ch1) - const)
               if(len_trim(ch3) == 2) then
                  str(j:j+1) = trim(ch3)
                  j = j + 2
               else if(l == 2 .and. n < 3) then
                  str(j:j+1) = ch3(1:2)
                  j = j + 2
               else if(l == 2 .and. n > 7) then
                  str(j:j+1) = ch3(1:2)
                  j = j + 2
               else if(l == 2 .and. n == 3) then
                  if(ch3(1:1) == '2' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 2 .and. n == 7) then
                  if(ch3(1:1) == '2' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 2 .and. n == 4) then
                  if(ch3(1:1) == '3' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'P') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'S') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 2 .and. n == 6) then
                  if(ch3(1:1) == '3' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'P') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'S') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 2 .and. n == 5) then
                  if(ch3(1:1) == '2' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 3 .and. n < 3) then
                  str(j:j+1) = ch3(1:2)
                  j = j + 2
               else if(l == 3 .and. n > 11) then
                  str(j:j+1) = ch3(1:2)
                  j = j + 2
               else if(l == 3 .and. n == 3) then
                  if(ch3(1:1) == '2' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 3 .and. n == 11) then
                  if(ch3(1:1) == '2' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 3 .and. n == 4) then
                  if(ch3(1:1) == '3' .and. ch3(2:2) == 'K') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'I') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'P') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'L') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'I') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'S') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 3 .and. n == 10) then
                  if(ch3(1:1) == '3' .and. ch3(2:2) == 'K') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'I') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '3' .and. ch3(2:2) == 'P') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'L') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'I') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '1' .and. ch3(2:2) == 'S') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 3 .and. n == 5) then
                  if(ch3(1:1) == '4' .and. ch3(2:2) == 'K') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'I') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'P') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'M') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'L') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'K') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'I') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'P') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else if(l == 3 .and. n == 9) then
                  if(ch3(1:1) == '4' .and. ch3(2:2) == 'K') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'I') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '4' .and. ch3(2:2) == 'P') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'M') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'L') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'K') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'I') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'H') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'G') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'F') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'D') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else if(ch3(1:1) == '2' .and. ch3(2:2) == 'P') then
                    str(j:j+2) = ch3
                    j = j + 3
                  else
                    str(j:j+1) = ch3(1:2)
                    j = j + 2
                  end if
               else
                  str(j:j+2) = ch3
                  j = j + 3
               end if
            endif
         endif
!
!  -----  If i=1 or Qi=4l+2 and i<>m,
!           insert '.'; else _RESULTANTi.  -----
!
         if(i/=1 .and. previous_q==full1 .and. n==1 .and. m/=i .and.      &
                                                       icore<2) then
            continue
         else if(i/=1 .and. previous_q==full1 .and. n==full-1 .and.       &
                                            m/=i .and. icore<2) then
            continue
         else if(i/=1 .and. previous_q==full1 .and. icore==0 .and.        &
                                                          m/=i) then
            continue
         else if (i/=1 .and. n/=full .or. m==i) then
            ch3 = couple(m+i-1)
            ch1 = ch3(2:2)
            if (ch1>='a' .and. ch1<='z') ch3(2:2)=char(ichar(ch1)-const)
            k = 2
            if (m == 1) k = 3
            str(j:j+k) = '_'//ch3(1:k)
            j = j + k + 1
         endif
         if (i == m) cycle
         str(j:j) = '.'
         j = j + 1
      end do
!
!>>>>>    Because of a compiler error on the SUN, the following is
!         needed to have the string printed correctly.
!         STR = STR(1:J)
      return
      end subroutine packls
