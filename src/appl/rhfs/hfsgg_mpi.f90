!***********************************************************************
!                                                                      *
      subroutine hfsgg_mpi
!                                                                      *
!   MPI version of HFSGG. Each rank processes a strided subset of the  *
!   outer IC loop and the accumulated hyperfine constants are reduced   *
!   over all MPI ranks before the root rank writes the output files.    *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnw
      use memory_man
      use decide_c
      use def_c,           only: aumamu, cvac, empam, rbcm, aucm,     &
                                 ccms, b1=>aumamu
      use eigv_c
      use foparm_c
      use jlabl_c,               labj=>jlbr, labp=>jlbp
      use mpi_c,           only: myid, nprocs, mpi_comm_world
      use nsmdat_c,        only: sqn, dmomnm, qmomb,                  &
                                 hfsi=>sqn, hfsd=>dmomnm, hfsq=>qmomb
      use orb_c
      use prnt_c
      use opt6_c
      use syma_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rinthf_i
      use rint_i
      use matelt_i
      use ispar_i
      use itjpo_i
      use itrig_i
      use oneparticlejj_i
      use gracah1_i
      use mpiu_i, only: gdrsummpi_root
      implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: cutoff = 1.0d-10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ffmin, ffmax, ff, i, j, kt, ipt, ic, ir, itjpoc, itjpor&
         , idiff, ia, ib, k, kk, loc1, loc2, ii, jj, jjii, jb, ja, jjb, jja, iflag, &
         occ_ic, occ_ir, ndiff, diff_pos, diff_neg
      real(double), dimension(nnnw) :: tshell
      real(double), dimension(2,nnnw,nnnw) :: rintme, amelt, elemfac
      real(double), dimension(nnnw,nnnw) :: rintgj, rintdgj, gjmelt, dgjmelt, &
         gjfac, dgjfac
!     .. Local pointer arrays
      real(double), dimension(:,:), pointer :: hfc
      real(double),  dimension(:), pointer :: gjc, dgjc
      real(double) ::  apart, gjpart, dgjpart, elemnt,  elemntgj, elemntdgj,&
         contr, contrgj, contrdgj, aumhz, barnau, dnmau, gfac, hfac, fj, &
         gja1, afa1, afa2, bfa1, bfa2, bfa3, gj, dgj, tilde1, &
         tilde2, factor1, factor2, rac1, rac2, hfselt1, hfselt2, evecic1, evecir1,&
         evecic2, evecir2
!-----------------------------------------------
!
!
!   Matrix elements smaller than CUTOFF are not accumulated
!
!
!   Allocate storage for local arrays
!
      call alloc (hfc, 5, nvec*nvec, 'HFC', 'HFS_MPI')
      call alloc (gjc, nvec*nvec, 'GJC', 'HFS_MPI')
      call alloc (dgjc, nvec*nvec, 'DGLC', 'HFS_MPI')
!
!   Initialise
!
      hfc(:,:nvec*nvec) = 0.0d00
!
      gjc(:nvec*nvec) = 0.0d00
      dgjc(:nvec*nvec) = 0.0d00
      gjfac(:,:) = 0.0d00
      dgjfac(:,:) = 0.0d00

!
!   Calculate and save the radial integrals and angular
!   matrix elements for the two multipolarities
!
      do kt = 1, 2
         do i = 1, nw
            do j = 1, nw
               if (kt == 1) then
                  rintme(kt,i,j) = rinthf(i,j,-2)
                  rintgj(i,j) = rinthf(i,j,1)
               else
                  rintme(kt,i,j) = rint(i,j,-3)
                  rintdgj(i,j) = rint(i,j,0)
               endif
               call matelt (i, kt, j, apart, gjpart, dgjpart)
               amelt(kt,i,j) = apart
               elemfac(kt,i,j) = apart*rintme(kt,i,j)
               if (kt /= 1) cycle
               gjmelt(i,j) = gjpart
               dgjmelt(i,j) = dgjpart
            end do
         end do
      end do
      do i = 1, nw
         do j = 1, nw
            gjfac(i,j) = gjmelt(i,j)*rintgj(i,j)
            dgjfac(i,j) = dgjmelt(i,j)*rintdgj(i,j)
         end do
      end do
!
!   Set the parity of the one-body operators
!
      ipt = 1
!
!   Sweep through the Hamiltonian matrix to determine the
!   diagonal and off-diagonal hyperfine constants
!
      do ic = myid + 1, ncf, nprocs
         itjpoc = itjpo(ic)
!
         do ir = 1, ncf
!
!   If LFORDR is .TRUE., a `first order' calculation is indicated;
!   only the CSFs with serial numbers exceeding IC are treated specially
!   only diagonal elements are evaluated for the `first order' CSFs
!
            if (lfordr .and. ic>iccut .and. ic/=ir) cycle
!
            itjpor = itjpo(ir)
            idiff = itjpoc - itjpor
!
!   Loop over the multipolarities
!
            do kt = 1, 2
!
!   Initialise the accumulator
!
               elemnt = 0.0d00
               elemntgj = 0.0d00
               elemntdgj = 0.0d00
!
!   Consider 3 cases
!               (k)
!   (1) < J || T   || J >    , k = 1,2  and IR >= IC
!
!               (k)
!   (2) < J || T   || J-1 >  , k = 1,2
!
!               (k)
!   (3) < J || T   || J-2 >  , k = 2
!
               if (.not.(idiff==0 .and. ir>=ic .or. idiff==2 .or. idiff==4&
                   .and. kt==2)) cycle
!
               if (itrig(itjpoc, itjpor, 2*kt + 1) == 0) then
                  cycle
               endif
               if (ipt /= 0 .and. ispar(ic)*ispar(ir)*ipt /= 1) then
                  cycle
               endif
               ndiff = 0
               diff_pos = 0
               diff_neg = 0
               do i = 1, nw
                  occ_ic = iqa(i,ic)
                  occ_ir = iqa(i,ir)
                  if (occ_ic == occ_ir) cycle
                  ndiff = ndiff + 1
                  if (ndiff > 2) exit
                  if (iabs(occ_ic - occ_ir) > 1) then
                     ndiff = 3
                     exit
                  endif
                  if (occ_ic > occ_ir) then
                     diff_pos = i
                  else
                     diff_neg = i
                  endif
               end do
               if (ndiff /= 0 .and. ndiff /= 2) then
                  cycle
               endif
               if (ndiff == 2) then
                  if (abs(elemfac(kt,diff_pos,diff_neg)) <= cutoff) then
                     if (.not.(kt == 1 .and. idiff == 0 .and. &
                         (abs(gjfac(diff_pos,diff_neg)) > cutoff .or. &
                          abs(dgjfac(diff_pos,diff_neg)) > cutoff))) then
                        cycle
                     endif
                  endif
               endif
               call oneparticlejj(kt,ipt,ic,ir,ia,ib,tshell)
!GG               CALL TNSRJJ (KT, IPT, IC, IR, IA, IB, TSHELL)
!
!   Accumulate the contribution from the one-body operators;
!
               if (ia /= 0) then
                  if (ia == ib) then
                     if (kt/=1 .or. idiff/=0) then
                        do ia = 1, nw
                           if (abs(tshell(ia)) <= cutoff) cycle
                           if (abs(elemfac(kt,ia,ia)) <= cutoff) cycle
                           elemnt = elemnt + elemfac(kt,ia,ia)*tshell(ia)
                           cycle
                        end do
                     else
                        do ia = 1, nw
                           if (abs(tshell(ia)) <= cutoff) cycle
                           if (abs(elemfac(kt,ia,ia)) <= cutoff .and. &
                               abs(gjfac(ia,ia)) <= cutoff .and. &
                               abs(dgjfac(ia,ia)) <= cutoff) cycle
                           elemnt = elemnt + elemfac(kt,ia,ia)*tshell(ia)
                           elemntgj = elemntgj + gjfac(ia,ia)*tshell(ia)
                           elemntdgj = elemntdgj + dgjfac(ia,ia)*tshell(ia)
                        end do
                     endif
                  else
                     if (abs(tshell(1)) > cutoff) then
                        elemnt = elemnt + elemfac(kt,ia,ib)*tshell(1)
                        if (kt==1 .and. idiff==0) then
                           elemntgj = elemntgj + gjfac(ia,ib)*tshell(1)
                           elemntdgj = elemntdgj + dgjfac(ia,ib)*tshell(1)
                        endif
                     endif
                  endif
               endif
!
!   Multiply with the configuration expansion coefficients and add the
!   contributions from the matrix elements to obtain total contributions
!
               if (abs(elemnt) <= cutoff .and. abs(elemntgj) <= cutoff .and. &
                   abs(elemntdgj) <= cutoff) then
                  cycle
               endif
               do k = 1, nvec
                  loc1 = (k - 1)*ncf
                  evecic1 = evec(ic + loc1)
                  evecir1 = evec(ir + loc1)
                  do kk = 1, nvec
                     loc2 = (kk - 1)*ncf
                     evecic2 = evec(ic + loc2)
                     evecir2 = evec(ir + loc2)
                     if (idiff==0 .and. ir/=ic) then
                        contr = elemnt*(evecic1*evecir2 + evecir1*evecic2)
                        contrgj = elemntgj*(evecic1*evecir2 + evecir1*evecic2)
                        contrdgj = elemntdgj*(evecic1*evecir2 + evecir1*evecic2)
                     else
                        contr = elemnt*evecic1*evecir2
                        contrgj = elemntgj*evecic1*evecir2
                        contrdgj = elemntdgj*evecic1*evecir2
                     endif
!
!   Magnetic dipole and the two operators of the g_j factor
!
                     if (kt == 1) then
                        if (idiff == 0) then
                           hfc(1,nvec*(k-1)+kk) = hfc(1,nvec*(k-1)+kk) + contr
                           gjc(nvec*(k-1)+kk) = gjc(nvec*(k-1)+kk) + contrgj
                           dgjc(nvec*(k-1)+kk) = dgjc(nvec*(k-1)+kk) + contrdgj
                        else if (itjpoc - itjpor == 2) then
                           hfc(2,nvec*(k-1)+kk) = hfc(2,nvec*(k-1)+kk) + contr
                        endif
!
!   Electric quadrupole
!
                     else if (kt == 2) then
                        select case (idiff)
                        case (0)
                           hfc(3,nvec*(k-1)+kk) = hfc(3,nvec*(k-1)+kk) + contr
                        case (2)
                           hfc(4,nvec*(k-1)+kk) = hfc(4,nvec*(k-1)+kk) + contr
                        case (4)
                           hfc(5,nvec*(k-1)+kk) = hfc(5,nvec*(k-1)+kk) + contr
                        end select
                     endif
                  end do
               end do
!
            end do
!
         end do
      end do
!
      call gdrsummpi_root(hfc(:,:), 5*nvec*nvec)
      call gdrsummpi_root(gjc(:), nvec*nvec)
      call gdrsummpi_root(dgjc(:), nvec*nvec)
!
      if (myid /= 0) go to 900
!
!   These are the conversion factors to obtain the hyperfine
!   constants in MHz
!
      aumhz = aucm*ccms*1.0d-06
      barnau = 1.0d-24/rbcm**2
      dnmau = b1/(2.0d00*cvac*empam)
!
      gfac = aumhz*dnmau*hfsd/hfsi
      hfac = aumhz*2.0d00*hfsq*barnau
!
!   Output the hyperfine interaction constants
!
      write (24, 302)
      write (29, 402)
!
      do i = 1, nvec
         do ii = 1, nvec
!
            jj = iatjpo(i)
            jjii = iatjpo(ii)
!
            if (.not.(jj==jjii .and. jjii>1 .or. jj>jjii)) cycle
!
            fj = 0.5d00*dble(jj - 1)
!
            gja1 = sqrt(1.0d00/(fj*(fj + 1.0d00)))

            afa1 = gfac*gja1

            if (jj == 2) then
               afa2 = 0.d000
            else
               afa2 = gfac*sqrt(1.0d00/(fj*(2.0d00*fj - 1.0d00)))
            endif
            bfa1 = hfac*sqrt((fj*(2.0d00*fj - 1.0d00))/((fj + 1.0d00)*(2.0d00*&
               fj + 3.0d00)))
            if (jj == 2) then
               bfa2 = 0.0d00
            else
               bfa2 = 0.25d00*hfac*sqrt((fj*(fj - 1.0d00))/((fj + 1.0d00)*(&
                  2.0d00*fj - 1.0d00)))
            endif
            if (jj == 4) then
               bfa3 = 0.0d00
            else
               bfa3 = 0.125d00*hfac*sqrt((fj*(fj - 1.0d000)*(2.0d00*fj - 1.0d00&
                  ))/(2.0d00*fj - 3.0d00))
            endif
!
!   Diagonal (J,J) A and B factors
!
            select case (jj - jjii)
            case (0)
               if (i <= ii) then
!                GJ = CVAC*GJA1*GJC(NVEC*(I-1)+II)
!                DGJ = 0.001160D0*GJA1*DGJC(NVEC*(I-1)+II)

                  write (24, 303) ivec(i), labj(jj), labp((iaspar(i)+3)/2), &
                     ivec(ii), labj(jjii), labp((iaspar(ii)+3)/2), afa1*hfc(1,&
                     nvec*(i-1)+ii), bfa1*hfc(3,nvec*(i-1)+ii)
!
!   Output diagonal hfs and g_j factors to file <name>.h or <name>.ch
!
                  if (i == ii) then
                     gj = cvac*gja1*gjc(nvec*(i-1)+ii)
                     dgj = 0.001160d0*gja1*dgjc(nvec*(i-1)+ii)
                     write (29, 403) ivec(i), labj(jj), labp((iaspar(i)+3)/2), &
                        afa1*hfc(1,nvec*(i-1)+ii), bfa1*hfc(3,nvec*(i-1)+ii), &
                        gj, dgj, gj + dgj
                  endif
               endif
!
!   Off diagonal (J,J-1) A and B factors
!
            case (2)
               write (24, 303) ivec(i), labj(jj), labp((iaspar(i)+3)/2), ivec(&
                  ii), labj(jjii), labp((iaspar(ii)+3)/2), afa2*hfc(2,nvec*(i-1&
                  )+ii), bfa2*hfc(4,nvec*(i-1)+ii)
!
!   Off diagonal (J,J-2) B factor
!
            case (4)
               write (24, 303) ivec(i), labj(jj), labp((iaspar(i)+3)/2), ivec(&
                  ii), labj(jjii), labp((iaspar(ii)+3)/2), 0.0d00, bfa3*hfc(5,&
                  nvec*(i-1)+ii)
            end select
!
         end do
      end do
!
!   These are the factors needed to obtain the F-dependent hyperfine
!   matrix elements in Hartrees
!
      tilde1 = sqrt((hfsi + 1.0d00)/hfsi)*hfsd*0.5d00*b1/(empam*cvac)
      if (hfsi > 0.6d00) then
         tilde2 = sqrt((hfsi + 1.0d00)*(2.0d00*hfsi + 3.d0)/((2.0d00*hfsi - &
            1.0d00)*hfsi))
         tilde2 = tilde2*hfsq*0.5d00*barnau
      else
         tilde2 = 0.0d00
      endif
!
!     II = 2*I
!
      ii = nint(2.0d00*hfsi)
!
!   Calculate the F-dependent matrix elemnts
!
!   Loop over the states
!
      do jb = 1, nvec
         do ja = 1, nvec
            jjb = iatjpo(jb) - 1
            jja = iatjpo(ja) - 1
            if (.not.(jja==jjb .and. jja>0 .or. jjb>jja)) cycle
!
!   Determine the possible F quantum numbers for the matrix element
!
            ffmin = max(abs(jja - ii),abs(jjb - ii))
            ffmax = min(jja + ii,jjb + ii)
!
!   Loop over the possible F quantum numbers
!
            iflag = 0
            do ff = ffmin, ffmax, 2
!
!   Phase factor
!
               if (mod((ii + jja - ff)/2,2) == 1) then
                  factor1 = -1.0d00
               else
                  factor1 = 1.0d00
               endif
               factor2 = factor1*sqrt((dble(jjb) + 1.0d00)*(dble(ii) + 1.0d00))
!
!   Determine the Racah W coefficients.
!
!GG               CALL DRACAH (II, JJA, II, JJB, FF, 2, RAC1)
              call gracah1 (ii,jja,ii,jjb,ff,2,rac1)
!GG               CALL DRACAH (II, JJA, II, JJB, FF, 4, RAC2)
              call gracah1 (ii,jja,ii,jjb,ff,4,rac2)
!
!   Obtain and output matrix elements for J,J
!
               if (jja - jjb == 0) then
                  hfselt1 = factor2*rac1*hfc(1,nvec*(jb-1)+ja)*tilde1
                  hfselt2 = factor2*rac2*hfc(3,nvec*(jb-1)+ja)*tilde2
                  if (abs(hfselt1 + hfselt2) > cutoff*10.0d-05) then
                     if (iflag == 0) write (24, 304)
                     iflag = 1
                     write (24, 305) ivec(jb), labj(jjb+1), labp((iaspar(jb)+3)&
                        /2), ivec(ja), labj(jja+1), labp((iaspar(ja)+3)/2), &
                        labj(ff+1), hfselt1 + hfselt2
                  endif
!
!   Obtain and output matrix elements for J,J-1
!
               else if (abs(jja - jjb) == 2) then
                  hfselt1 = factor2*rac1*hfc(2,nvec*(jb-1)+ja)*tilde1
                  hfselt2 = factor2*rac2*hfc(4,nvec*(jb-1)+ja)*tilde2
                  if (abs(hfselt1 + hfselt2) > cutoff*10.0d-05) then
                     if (iflag == 0) write (24, 304)
                     iflag = 1
                     write (24, 305) ivec(jb), labj(jjb+1), labp((iaspar(jb)+3)&
                        /2), ivec(ja), labj(jja+1), labp((iaspar(ja)+3)/2), &
                        labj(ff+1), hfselt1 + hfselt2
                  endif
!
!   Obtain and output matrix elements for J,J-2
!
               else if (abs(jja - jjb) == 4) then
                  hfselt1 = 0.0d00
                  hfselt2 = factor2*rac2*hfc(5,nvec*(jb-1)+ja)*tilde2
                  if (abs(hfselt1 + hfselt2) > cutoff*10.0d-05) then
                     if (iflag == 0) write (24, 304)
                     iflag = 1
                     write (24, 305) ivec(jb), labj(jjb+1), labp((iaspar(jb)+3)&
                        /2), ivec(ja), labj(jja+1), labp((iaspar(ja)+3)/2), &
                        labj(ff+1), hfselt1 + hfselt2
                  endif
               else
                  hfselt1 = 0.0d00
                  hfselt2 = 0.0d00
               endif
            end do
         end do
      end do
!
  900 continue
      call dalloc (hfc, 'HFC', 'HFS_MPI')
      call dalloc (gjc, 'GJC', 'HFS_MPI')
      call dalloc (dgjc, 'DGJC', 'HFS_MPI')
      return
!
  302 format(/,/,' Interaction constants:'/,/,&
         ' Level1  J Parity  Level2  J Parity',8x,'A (MHz)',13x,'B (MHz)'/)
  402 format(/,/,' Interaction constants:'/,/,' Level1  J Parity ',8x,'A (MHz)'&
         ,13x,'B (MHz)',13x,'g_J',14x,'delta g_J',11x,'total g_J'/)
  303 format(1x,1i3,5x,2a4,2x,1i3,5x,2a4,1p,2d20.10)
  403 format(1x,1i3,5x,2a4,1p,5d20.10)
  304 format(/,/,' Matrix elements:'/,/,&
         ' Level1  J Parity  Level2  J Parity    F ',4x,'Matrix element (a.u.)'&
         /)
  305 format(1x,1i3,5x,2a4,2x,1i3,5x,2a4,1x,a4,4x,1p,1d20.10)
      return
!
      end subroutine hfsgg_mpi
