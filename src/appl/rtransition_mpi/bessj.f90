!***********************************************************************
!                                                                      *
      subroutine bessj(w)
!                                                                      *
!   This routine evaluates Bessel fuctions J K ( W*R/C ) at the grid   *
!   points for K=L-1,L,L+1 and stores  them in the  arrays BJ(..,1),   *
!   BJ(..,2),BJ(..,3) respectively. It uses a power series expansion   *
!   for  small  r and switches to sin/cos expansion when more than 4   *
!   terms in power series are required.                                *
!                                         Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use bess_c, only: bj, tc, td
      use debug_c, only: ldbpr
      use grid_c
      use osc_c, only: lk, kk, l=>lk
      use def_c, only: speed_of_light_au, vacuum_speed_of_light_au, pi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64) , intent(in) :: w
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i1, iw, nn, iend, iprod, i, jchan, j, iswap, modnn4, &
                 jj, jjj, ii
      real(real64) :: epsi, dfn, wr, wa, xbess1, s1, ssn, scn, sn, cn, b, skeep
!-----------------------------------------------
!
!
      epsi = 1.0e-5_real64
!
      iw = 1
      nn = l - 1
      if (kk /= 0) then
         iw = 2
         nn = l
      end if
    1 continue
      iend = 2*nn + 1
      iprod = 1
      do i = 1, iend, 2
         iprod = iprod*i
      end do
      dfn = iprod
      jchan = radial_point_count
      bj(1,iw) = 1.0_real64
      do j = 2, radial_point_count
         wr = w*radial_grid(j)
         wa = -wr*wr*0.5_real64
         xbess1 = 1.0_real64
         s1 = 1.0_real64
         do i = 1, 4
            xbess1 = xbess1*wa/dble(i*(2*(nn + i) + 1))
            s1 = s1 + xbess1
            if (abs(xbess1) < abs(s1)*epsi) goto 4
         end do
         jchan = j
         exit
    4    continue
         bj(j,iw) = s1*wr**nn/dfn
      end do
!
!   Use sin/cos expansion when power series takes longer
!   than 4 terms to converge
!
      if (jchan < radial_point_count) then
         iswap = 0
         modnn4 = mod(nn - 1,4) + 1
         select case (modnn4)
!
!   NN = 1, 5, 9, ...
!
         case default
            ssn = -1.0_real64
            scn = 1.0_real64
            iswap = 1
!
!   N = 2, 6, 10,....
!
         case (2)
            ssn = -1.0_real64
            scn = -1.0_real64
!
!   NN = 3, 7, 11,...
!
         case (3)
            ssn = 1.0_real64
            scn = -1.0_real64
            iswap = 1
!
!   NN = 0, 4, 8,...
!
         case (4)
            ssn = 1.0_real64
            scn = 1.0_real64
         end select
!
   13    continue
         do j = jchan, radial_point_count
            wa = w*radial_grid(j)
            if (iswap <= 0) then
               sn = ssn*sin(wa)
               cn = scn*cos(wa)
            else
               sn = ssn*cos(wa)
               cn = scn*sin(wa)
            end if
            i = -1
            s1 = 0.0_real64
            i = i + 1
            i1 = i
            do i = i1, nn
               if (i == 0) then
                  b = 1.0_real64/wa
               else
                  b = b*dble((nn + i)*(nn - i + 1))/dble(2*i)/wa
               end if
               s1 = s1 + b*sn
               skeep = sn
               sn = cn
               cn = -skeep
            end do
            bj(j,iw) = s1
         end do
      end if
      if (nn>=l + 1 .or. kk==1) then
!
!   Print out Bessel functions if (debug) option set
!
         if (ldbpr(16)) then
            do jj = 1, 3
               jjj = l - 2 + jj
               write (99, 300) jjj, (bj(ii,jj),ii=1,radial_point_count)
            end do
         end if
!
!   All done
!
! zou
         do jj = 1, 3
            jjj = l - 2 + jj
            bj(:radial_point_count,jj) = bj(:radial_point_count,jj)*(speed_of_light_au/vacuum_speed_of_light_au)**jjj
         end do
! zou
         return
      else
         nn = nn + 1
         iw = iw + 1
         goto 1
      end if
!
  300 format(/," Bessel function of order ",i3,/(1p,7d18.10))
      return
!
      end subroutine bessj
