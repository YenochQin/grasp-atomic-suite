!***********************************************************************
!                                                                      *
      subroutine biotr1(pi, qi, nli, ninshli, pf, qf, nlf, ninshlf, &
             ngrid, mxl, scr, lscr, ntestg, cishl, cici, cfshl, cfci)
!                                                                      *
! Generate Matrices for rotating radial functions and                  *
! for counter rotating CI coefficients                                 *
!                                                                      *
! This subroutine is a modification of the corresponding MCHF          *
! routine                                                              *
!                                                                      *
! Per Jonsson,     Lund June 1996                                      *
!                                                                      *
! This subroutine has been modified to support transformations         *
! in cases like 3s2p1d 3s2p                                            *
!                                                                      *
! Per Jonsson,     Lund Feb  1997                                      *
!                                                                      *
! =====                                                                *
! Input                                                                *
! =====                                                                *
!                                                                      *
! PI,QI   : Input radial shells for initial state                      *
!         : (large and small component)                                *
! NLI     : Number of shells per kappa for initial state               *
! NINSHLI : Number of inactive shells  per kappa for initial state     *
!                                                                      *
! PF,QF   : Input radial shells for final   state                      *
!         : (large and small component)                                *
! NLF     : Number of shells per kappa for final   state               *
! NINSHLF : Number of inactive shells  per kappa for final state       *
!                                                                      *
! NGRID   : Number of gridpoints                                       *
! MXL     : Number of kappa quantum numbers                            *
! SCR     :                                                            *
! LSCR    : Total length of scratch space.                             *
! NTESTG  : Global PRINT flag : = 0 => complete silence                *
!                               = 1 => test and PRINT overlap matrices,*
!                                      PRINT header and  t matrix      *
!                            .gt. 1 => Hope you now what you are doing *
!                                                                      *
! ======                                                               *
! Output                                                               *
! ======                                                               *
!                                                                      *
! CISHL : Rotation matrix for initial state shells                     *
! CICI  : Rotation matrix for initial state CI coefficients            *
! CFSHL : Rotation matrix for final state shells                       *
! CFCI  : Rotation matrix for final state CI coefficients              *
! PI,QI : Initial shells in biorthogonal basis                         *
!         biorthogonal basis                                           *
! PF,QF : Final   shells in biorthogonal basis                         *
!                                                                      *
!                                                                      *
! Inactive Shells : The functions of the inactive shells must be       *
!                   be supplied to the program, and NLI,NLF            *
!                   must refer to the total number of                  *
!                   occupied shells ( inactive+active)                 *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ifnmnx_i
      use ielsum_i
!      USE gets_I
!      USE wrtmat_I
!      USE copvec_I
!      USE invmat_I
!      USE ulla_I
!      USE trpmat_I
!      USE matml4_I
!      USE scalve_I
!      USE setvec_I
!      USE pamtmt_I
      use copvec_i, only: copvec
      use gets_i, only: gets
      use invmat_i, only: invmat
      use matml4_i, only: matml4
      use pamtmt_i, only: pamtmt
      use scalve_i, only: scalve
      use setvec_i, only: setvec
      use trpmat_i, only: trpmat
      use ulla_i, only: ulla
      use wrtmat_i, only: wrtmat
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ngrid,  mxl
      integer, intent(in) :: lscr,  ntestg
      integer, dimension(mxl)  :: nli, ninshli, nlf, ninshlf
!     The following have dimensions (NGRID:*)
!GG      REAL(DOUBLE), DIMENSION(:,:), pointer :: PI, QI, PF, QF
      real(double), dimension(ngrid,*) :: pi, qi, pf, qf
!     The following have fixed dimensions (NGRID, NLMAX)
      real(double), dimension(*)  :: cishl, cfshl
!     The work array
!     REAL(DOUBLE), DIMENSION(LWORK1) :: work
!     Th following have fixed dimensions (20*NLMAX*NLMAX)
      real(double), dimension(*) :: cici, cfci
      real(double), dimension(*) :: scr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntestl, ntest, ili, ilf, nlimx, nlfmx, nlifmx, nlti, &
                 nltf, kfree, klpi, klpf, klpif, klstot, klsif, klsifi, &
                 klci, klcf, klscr, l, iioff, ifoff, ni, nf, nifmn, iii,&
                 jjj, klpmx, klpmn, nmx, nmn, klcmx, klcmn, ndiff, j, &
                 klti, i, kltf
      real(double) :: tii, tiii
!-----------------------------------------------
!
      ntestl = 0
      ntest = max(ntestl,ntestg)
      if (ntest >= 1) then
         write (6, *)
         write (6, *) '                   *************************'
         write (6, *) '                   *   Entering BIOTR1     *'
         write (6, *) '                   *************************'
         write (6, *)
      endif
!
!. Scratch should at least be of length
!
      ili = 1
      ilf = 1
!
!. Largest number of shells of a given symmetry
!
      nlimx = ifnmnx(nli,mxl,1)
      nlfmx = ifnmnx(nlf,mxl,1)
      nlifmx = max(nlimx,nlfmx)
      if (ntest >= 10) write (6, *) ' NLIMX,NLFMX NLIFMX ', nlimx, nlfmx, &
         nlifmx
!
!. Total numner of shells
!
      nlti = ielsum(nli,mxl)
      nltf = ielsum(nlf,mxl)
      if (ntest >= 10) write (6, *) ' NLTI NLTF', nlti, nltf
!
! Scratch space for orbital rotations
      kfree = 1
!
      klpi = kfree
      kfree = kfree + nlimx*ngrid
      write (6, *) ' In biotrn: KLPI    = ', klpi
      write (6, *) '            KFREE   = ', kfree
!
      klpf = kfree
      kfree = kfree + nlfmx*ngrid
      write (6, *) ' In biotrn: KLPF    = ', klpf
      write (6, *) '            KFREE   = ', kfree
!
      klpif = kfree
      kfree = kfree + nlifmx*ngrid
      write (6, *) ' In biotrn: KLPIF   = ', klpif
      write (6, *) '            KFREE   = ', kfree
!
!. Total overlap matrix
!
      klstot = kfree
      kfree = kfree + nlti*nltf
      write (6, *) ' In biotrn: KLSTOT  = ', klstot
      write (6, *) '            KFREE   = ', kfree
!
      klsif = kfree
      kfree = kfree + nlifmx**2
      write (6, *) ' In biotrn: KLSIF   = ', klsif
      write (6, *) '            KFREE   = ', kfree
!
      klsifi = kfree
      kfree = kfree + nlifmx**2
      write (6, *) ' In biotrn: KLSIFI  = ', klsifi
      write (6, *) '            KFREE   = ', kfree
!
      klci = kfree
      kfree = kfree + nlifmx**2
      write (6, *) ' In biotrn: KLCI    = ', klci
      write (6, *) '            KFREE   = ', kfree
!
      klcf = kfree
      kfree = kfree + nlifmx**2
      write (6, *) ' In biotrn: KLCF    = ', klcf
      write (6, *) '            KFREE   = ', kfree
!
      klscr = kfree
      kfree = kfree + nlifmx**2 + nlifmx*(nlifmx + 1)
      write (6, *) ' In biotrn: KLSCR   = ', klscr
      write (6, *) '            KFREE   = ', kfree
      write (6, *) '         =>  FREE =     ', kfree
!
!. Check length of scratch
!
      if (lscr <= kfree - 1) then
       write (6,*)' BIOTR1 in trouble ! '
       write (6,*)' Increase dimension of scratch before call to BIOTR1'
       write (6,*)' Current and required length (LSCR,KFREE-1)', lscr, &
            kfree - 1
       stop 'Increase LWORK before call to BIOTR1'
       endif
!
!. Obtain overlap matrix
!
      call gets (scr(klstot), nlti,nltf)

      do l = 1, mxl
         if (ntest >= 5) then
            write (6, *) '   L = ', l
            write (6, *) '       Orbital rotation...'
         endif
!
!. Offsets for given L in shell matrices
!
         if (l == 1) then
            iioff = 1
            ifoff = 1
         else
            iioff = iioff + nli(l-1)**2
            ifoff = ifoff + nlf(l-1)**2
         endif
         if (ntest >= 1) then
           write (6, *)
           write (6, *) &
           ' BIOTRN : Information on transformations of shells with L =',l
           write (6, *)
         endif
!
! =========================================================
! 1 : Obtain Biorthogonal forms of initial and final shells
! =========================================================
!
!. The overlap matrix can be written
!
!     * * * * * * *
!     *       *   *
!     *   S   * X *
!     *       *   *
!     * * * * * * *
!
! where S is the quadratic subblock.
! The basis functions corresponding to the quadratic
! subblock are made biorthogonal with an UL decomposition.
! The remaining basis
! functions are made biorthogonal by choosing  the
! transformation
!            * * *
!            *   *
!            * Y *
!            *   *
!            * * *
!            * 1 *
!            * * *
!
!
! With Y = -S-1*X
!
         ni = nli(l)
         nf = nlf(l)
!
!.1.1 : Obtain shells of given L in proper order for
!       biorthogonal treatment
!
         if (l == 1) then
            ili = 1
            ilf = 1
         else
            ili = ili + nli(l-1)
            ilf = ilf + nlf(l-1)
         endif

!
! 1.2 obtain biorthogonal of the first min(ni,nf) shells
!
         nifmn = min(ni,nf)
!
!. Overlap matrix SIF = Integral (PI(I)*PF(J))
!
!ww Per change to support cases like 3s2p1d 3s2p. The belonging endif
!ww Is just at the end

         if (nifmn <= 0) cycle

         do iii = 1, nifmn
            do jjj = 1, nifmn
               scr(klsif+(jjj-1)*nifmn+iii-1) = scr(klstot-1+(jjj+ilf-1-1)*nlti&
                  +iii+ili-1)
            end do
         end do
         if (ntest >= 15) then
            write (6, *) ' Overlap matrix '
            call wrtmat (scr(klsif), nifmn, nifmn, nifmn, nifmn)
         endif
!
! Obtain upper triangular CI and CF so CI(T) S CF = 1
! or CF CI(T) = S-1, which corresponds to an UL decomposition
!
!. Invert S
!
         call copvec (scr(klsif), scr(klsifi), nifmn**2)
         call invmat (scr(klsifi), scr(klci), nifmn, nifmn)

!. UL decompose
         call copvec (scr(klsifi), scr(klsif), nifmn**2)
         call ulla (scr(klsif), scr(klcf), scr(klci), nifmn, scr(klscr))
         call trpmat (scr(klci), nifmn, nifmn, scr(klscr))
         call copvec (scr(klscr), scr(klci), nifmn**2)
!
!. The transformation matrix between the first NIFMX
!. shells is now known, biorthogonalize remaining orbitals
!
         if (ni/=nf .and. ni/=0 .and. nf/=0) then
            if (ni > nf) then
               klpmx = klpi
               klpmn = klpf
               nmx = ni
               nmn = nf
               klcmx = klci
               klcmn = klcf
            else
               klpmx = klpf
               klpmn = klpi
               nmx = nf
               nmn = ni
               klcmx = klcf
               klcmn = klci
            endif
            ndiff = nmx - nmn
!
! Y = -S-1 * X
!. overlap X between remaining orbitals and the other set
!
            if (ni > nf) then
!
! I columns F rows
!
               do iii = nmn + 1, nmx
                  do jjj = 1, nf
                     scr(klsif+(iii-nmn-1)*nf+jjj-1) = scr(klstot-1+(jjj+ilf-1-&
                        1)*nlti+iii+ili-1)
                  end do
               end do
            else if (nf > ni) then
! F columns I rows
               do jjj = nmn + 1, nmx
                  do iii = 1, ni
                     scr(klsif+(jjj-nmn-1)*ni+iii-1) = scr(klstot-1+(jjj+ilf-1-&
                        1)*nlti+iii+ili-1)
                  end do
               end do
            endif
!
            if (ni > nf) then
               call trpmat (scr(klsifi), nmn, nmn, scr(klscr))
               call copvec (scr(klscr), scr(klsifi), nmn**2)
            endif
            call matml4 (scr(klscr), scr(klsifi), scr(klsif), nmn, ndiff, nmn, &
               nmn, nmn, ndiff, 0)
            call scalve (scr(klscr), -1.0d0, nmn*ndiff)
            call copvec (scr(klscr), scr(klsif), nmn*ndiff)
!
! Construct complete CMX
!
            call setvec (scr(klscr), 0.0d0, nmx**2)
            do j = 1, nmx
               if (j <= nifmn) then
                  call copvec (scr(klcmx+(j-1)*nifmn), scr(klscr+(j-1)*nmx), &
                     nmn)
               else
                  call copvec (scr(klsif+(j-nmn-1)*nmn), scr(klscr+(j-1)*nmx), &
                     nmn)
                  scr(klscr-1+(j-1)*nmx+j) = 1.0d0
               endif
            end do
!
            call copvec (scr(klscr), scr(klcmx), nmx**2)
         endif
!ww Pertest
!         ENDIF
!
!. The two upper triangular matrices CI and CF are now known
!. Transfer to permanent arrays
!
         call copvec (scr(klci), cishl(iioff), ni**2)
         call copvec (scr(klcf), cfshl(ifoff), nf**2)
!
!. Rotate the large component of the shells
!
         call copvec (pi(1,ili), scr(klpi), ni*ngrid)
         call copvec (pf(1,ilf), scr(klpf), nf*ngrid)
!
         write (*, *) 'Transformation matrices initial'
         call wrtmat (scr(klci), ni, ni, ni, ni)
         call matml4 (scr(klpif), scr(klpi), scr(klci), ngrid, ni, ngrid, ni, &
            ni, ni, 0)
         call copvec (scr(klpif), pi(1,ili), ni*ngrid)
         write (*, *) 'Transformation matrices final'
         call wrtmat (scr(klcf), nf, nf, nf, nf)
         call matml4 (scr(klpif), scr(klpf), scr(klcf), ngrid, nf, ngrid, nf, &
            nf, nf, 0)
         call copvec (scr(klpif), pf(1,ilf), nf*ngrid)
!
!. Rotate the small component of the shells
!
         call copvec (qi(1,ili), scr(klpi), ni*ngrid)
         call copvec (qf(1,ilf), scr(klpf), nf*ngrid)
!
         call matml4 (scr(klpif), scr(klpi), scr(klci), ngrid, ni, ngrid, ni, &
            ni, ni, 0)
         call copvec (scr(klpif), qi(1,ili), ni*ngrid)
         call matml4 (scr(klpif), scr(klpf), scr(klcf), ngrid, nf, ngrid, nf, &
            nf, nf, 0)
         call copvec (scr(klpif), qf(1,ilf), nf*ngrid)
!
         if (ntest >= 1) then
            write (6, *) ' Test of overlap of biorthonormal functions'
! F columns I rows
            do jjj = 1, nf
               do iii = 1, ni
                  scr(klsif+(jjj-1)*ni+iii-1) = scr(klstot-1+(jjj+ilf-1-1)*nlti&
                     +iii+ili-1)
               end do
            end do
            call matml4 (scr(klscr), scr(klci), scr(klsif), ni, nf, ni, ni, ni&
               , nf, 1)
            call matml4 (scr(klsif), scr(klscr), scr(klcf), ni, nf, ni, nf, nf&
               , nf, 0)
            write (6, *) &
              ' new overlap matrix ( should be 1 on diag, 0 elsewhere )'
            call wrtmat (scr(klsif), ni, nf, ni, nf)
         endif

         if (ntest >= 1) then
            write (6, *)
            write (6, *) ' Orbital Rotation matrix for I state'
            call wrtmat (cishl(iioff), ni, ni, ni, ni)
            write (6, *) ' Orbital Rotation matrix for F state'
            call wrtmat (cfshl(ifoff), nf, nf, nf, nf)
            write (6, *)
         endif
!
!. Matrix for counterrotation of CI coefficients, initial state
!
         klti = klsif
         call pamtmt (scr(klci), scr(klti), scr(klscr), ni)
         do i = 1, ni
            tii = scr(klti-1+(i-1)*ni+i)
            tiii = 1.0d0/tii
            call scalve (scr(klti+(i-1)*ni), tiii, i - 1)
         end do
         call copvec (scr(klti), cici(iioff), ni*ni)
!
!. Matrix for counterrotation of CI coefficients, Final state
!
         kltf = klsif
         call pamtmt (scr(klcf), scr(kltf), scr(klscr), nf)
         do i = 1, nf
            tii = scr(kltf-1+(i-1)*nf+i)
            tiii = 1.0d0/tii
            call scalve (scr(kltf+(i-1)*nf), tiii, i - 1)
         end do
         call copvec (scr(kltf), cfci(ifoff), nf*nf)
         if (ntest < 1) cycle
         write (6, *)
         write (6, *) ' CI-Rotation matrix for I state'
         call wrtmat (cici(iioff), ni, ni, ni, ni)
         write (6, *) ' CI-Rotation matrix for F state'
         call wrtmat (cfci(ifoff), nf, nf, nf, nf)
         write (6, *)

      end do
!
      return
      end subroutine biotr1
