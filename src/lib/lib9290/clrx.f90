!***********************************************************************
!                                                                      *
      real(real64) function clrx (kappaa, k, kappab)
!                                                                      *
!   The value of CLRX is the 3-j symbol:                               *
!                                                                      *
!                    ( JA        K        JB  )                        *
!                    ( 1/2       0       -1/2 )                        *
!                                                                      *
!   The  K'S are kappa angular quantum numbers. The formula is taken   *
!   from D M Brink and G R Satchler, <Angular Momentum>, second edi-   *
!   tion (Oxford: Clarendon press, 1968), p 138.   The logarithms of   *
!   the first  MFACT  factorials must be available in  COMMON/FACTS/   *
!   for this program to function correctly. Note that  N!  is stored   *
!   in FACT(N+1)                                                       *
!                                                                      *
!   No subroutines called.                                             *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use facts_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: kappaa, k, kappab
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ka, kb, kapkb, kabkp, kamkb, kbmka
      real(real64) :: exptrm
!-----------------------------------------------
!
!
!   Determine the absolute values of the kappas
!
      ka = abs(kappaa)
      kb = abs(kappab)
!
!   Perform the triangularity check
!
      if (abs(ka - kb)<=k .and. ka+kb-1>=k) then
!
!   Triangularity satisfied; compute the 3j coefficient
!
!   Begin with the logarithm of the square of the leading term
!
         exptrm = -log(dble(ka*kb))
!
!   Compute the logarithm of the square root of the leading term
!   and the factorial part that doesn't depend on the parity of
!   KA+KB+K (the delta factor)
!
         kapkb = ka + kb
         kabkp = kapkb + k
         kamkb = ka - kb
         kbmka = kb - ka
         exptrm = 0.5_real64*(exptrm + gam(kapkb-k)+gam(kamkb+k+1)+gam(kbmka+k+1)-&
            gam(kabkp+1))
!
!   The remainder depends on the parity of KA+KB+K
!
         if (mod(kabkp,2) == 0) then
!
!   Computation for even parity case
!
!   Include the phase factor: a minus sign if necessary
!
            if (mod(3*kabkp/2,2) == 0) then
               clrx = 1.0_real64
            else
               clrx = -1.0_real64
            endif
!
!   Include the contribution from the factorials
!
            exptrm = exptrm + gam((kabkp+2)/2) - gam((kapkb-k)/2) - gam((kamkb+&
               k+2)/2) - gam((kbmka+k+2)/2)
!
         else
!
!   Computation for odd parity case
!
!   Include the phase factor: a minus sign if necessary
!
            if (mod((3*kabkp - 1)/2,2) == 0) then
               clrx = 1.0_real64
            else
               clrx = -1.0_real64
            endif
!
!   Include the contribution from the factorials
!
            exptrm = exptrm + gam((kabkp+1)/2) - gam((kapkb-k+1)/2) - gam((&
               kamkb+k+1)/2) - gam((kbmka+k+1)/2)
!
         endif
!
!   Final assembly
!
         clrx = clrx*exp(exptrm)
!
      else
!
!   Triangularity violated; set the coefficient to zero
!
         clrx = 0.0_real64
!
      endif
!
      return
!
      end function clrx
