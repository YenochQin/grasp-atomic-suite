!***********************************************************************
!                                                                      *
      subroutine breid(ja, jb, ja1, ipca, jb1)
!                                                                      *
!   Computes closed shell contributions - aaaa and exchange only.      *
!                                                                      *
!   Call(s) to: [LIB92]: CLRX, CXK, ITRIG, TALK, SNRC.                 *
!                                                                      *
!                                           LAST UPDATE: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use breit_core_occupations, only: icore, max_orbitals
      use numeric_constants, only: eleven, eps, four, half, one, seven, ten, tenth, three, two, zero
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use csf_pair_angular_state, only: jjc1, jjc2, jjq1, jjq2, jlist, klist, max_orbitals, ncore, npeel, nq1, nq2
      use orbital_csf_basis, only: orbital_principal_n, orbital_kappa
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use breit_radial_rank_ranges, only: snrc
      use relativistic_three_j, only: clrx
      use breit_angular_coefficient_buffer, only: talk
      use angular_momentum_triangle, only: itrig
      use breit_radial_integral_coefficients, only: cxk
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ja
      integer, intent(in) :: jb
      integer, intent(in) :: ja1
      integer, intent(in) :: ipca
      integer, intent(in) :: jb1
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: numax = 20
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: js(4), kaps(4), ks(4)
      integer :: ia1, ib1, isg, nqs1, nqs2, i, nd1, nd2, ne1, ne2, ibrd, ibre, &
         n, nu, k, kap1, itype, mu, ip, ipp, kk, nup1
      real(real64) :: cone(7, 20)
      real(real64) :: s(12)
      real(real64) :: const, gam, dksks, dnunu1, coef, proc, prod
!-----------------------------------------------
!
!   1.0  Initialization
!
      if (ipca == 2) then
         ia1 = klist(ja1)
      else
         ia1 = jlist(ja1)
      end if
      ib1 = klist(jb1)
!
      isg = 1
      if (ja == jb) then
         if (icore(ia1)/=0 .and. icore(ib1)/=0) then
            if (ja > 1) return
            isg = -1
         end if
      end if
!
      js(1) = ia1
      js(2) = ib1
      js(3) = ia1
      js(4) = ib1
      nqs1 = nq1(ia1)
      nqs2 = nq2(ib1)
      do i = 1, 4
         kaps(i) = 2*orbital_kappa(js(i))
         ks(i) = iabs(kaps(i))
      end do
      const = nqs1*nqs2
      if (ibug2 /= 0) write (99, '("BREID: orbitals ",2i3)') ia1, ib1
!
!   2.0  Set range of tensor indices
!
      call snrc (js, kaps, ks, nd1, nd2, ne1, ne2, ibrd, ibre)
      if (ibug2 /= 0) write (99, '(2x,"ND1 ND2 NE1 NE2 IBRD IBRE ",6i5)') nd1, nd2, ne1, ne2, ibrd, ibre
      if (ia1 == ib1) then
!
!   3.0 Calculate aaaa interaction
!
         do n = 1, nd2
            nu = nd1 + 2*(n - 1)
            k = nu
            if (mod(k,2) /= 1) return
            kap1 = kaps(1)/2
            gam = clrx(kap1,nu,kap1)
            dksks = ks(1)*ks(1)
            dnunu1 = nu*(nu + 1)
            coef = const*two*dksks*gam*gam/dnunu1
            if (ibug2 /= 0) write (99, '(2x,"aaaa contribution: NU,GAM,COEF",i5,2(3x,1p,d15.8))') nu, gam, coef
            itype = isg*4
            call talk (ja, jb, nu, ia1, ia1, ia1, ia1, itype, coef)
         end do
         return
      end if
!
!   Calculate exchange interactions
!
      if (ibre < 0) return
      if (ne2 > numax) then
         write (*, '("BREID: Dimension error for NUMAX.")')
         stop
      end if
!
      cone(:,:ne2) = zero
!
      proc = -const/dble(ks(1)*ks(2))
!
!   Negative sign arises from Pauli phase factor
!
      do n = 1, ne2
         nu = ne1 + 2*(n - 1)
         k = nu
         ip = (ks(1)-ks(2))/2 + k
         ipp = ip + 1
         if (nu /= 0) then
            kk = k + k + 1
            if (itrig(ks(1),ks(2),kk) /= 0) then
               prod = proc
               if (mod(ip,2) /= 0) prod = -prod
               call cxk (s, js, kaps, nu, k, ibre, 2)
               if (ibug2 /= 0) write (99, '(2x,"PROD = ",1p,d15.8,/," S",7d15.8)') prod, (s(mu),mu=1,3)
               cone(:3,n) = cone(:3,n) + prod*s(:3)
            end if
!
            k = nu - 1
            kk = k + k + 1
            if (itrig(ks(1),ks(2),kk) /= 0) then
               prod = proc
               if (mod(ipp,2) /= 0) prod = -prod
               call cxk (s, js, kaps, nu, k, ibre, 2)
               if (ibug2 /= 0) write (99, '(2x,"PROD = ",1p,d15.8,/," S",7d15.8)') prod, (s(mu),mu=1,3)
               cone(:3,n) = cone(:3,n) + prod*s(:3)
!
            end if
         end if
         if (n == ne2) exit
         k = nu + 1
         kk = k + k + 1
         prod = proc
         if (mod(ipp,2) /= 0) prod = -prod
         call cxk (s, js, kaps, nu, k, ibre, 2)
         if (ibug2 /= 0) write (99, '(2x,"PROD = ",1p,d15.8,/," S",7d15.8)') prod, (s(mu),mu=1,7)
         cone(:,n) = cone(:,n) + prod*s(:7)
      end do
!
!   4.0  Output results
!
      do n = 1, ne2
         nu = ne1 + 2*(n - 1)
         itype = isg*5
         call talk (ja, jb, nu, ib1, ia1, ib1, ia1, itype, cone(1,n))
         call talk (ja, jb, nu, ia1, ib1, ib1, ia1, itype, cone(2,n))
         call talk (ja, jb, nu, ia1, ib1, ia1, ib1, itype, cone(3,n))
         if (n == ne2) cycle
         nup1 = nu + 1
         itype = isg*6
         call talk (ja, jb, nup1, ia1, ib1, ia1, ib1, itype, cone(4,n))
         call talk (ja, jb, nup1, ib1, ia1, ib1, ia1, itype, cone(5,n))
         call talk (ja, jb, nup1, ia1, ib1, ib1, ia1, itype, cone(6,n))
         call talk (ja, jb, nup1, ib1, ia1, ia1, ib1, itype, cone(7,n))
      end do
      return
!
      return
!
      end subroutine breid
