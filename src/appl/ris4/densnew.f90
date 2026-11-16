!***********************************************************************
!                                                                      *
subroutine densnew(doit,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
!       SUBROUTINE DENSNEW(DOIT,DINT1,DINT2,DINT3,                      &
!                          DINT4,DINT5,DINT6,DINT7,                     &
!                          DINT1VEC,DENS1VEC)
!                                                                      *
!   IF angular coefficients must be calculated                         *
!   This routine controls combines the radial and angular parts for the*
!   calculation of the NMS parameter, the electron density at the      *
!   origin and radial expectation values.
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, GETYN                         *
!                        ITJPO, ONESCALAR                              *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!   Modified by C. Naz\'e  Feb. 2012                                   *
!   Modified by J. Ekman   Nov. 2013                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: orbital_index_radix, max_orbitals, max_radial_points
      use debug_c
      use decide_c
      use def_c
      use eigv_c
      use foparm_c
      use grid_c
      use jlabl_c
      use npar_c
      use orb_c
      use blk_c
      use prnt_c
      use teilst_c
      use buffer_c
      use ris_c
      use syma_c
      use prnt_c,           only : nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i
      use convrt_i
      use getyn_i
      use itjpo_i
      use onescalar_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint1, dint2, &
                                                        dint3, dint4, &
                                                        dint5, dint6, &
                                                        dint7
      integer, intent(in) :: doit
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(max_orbitals) :: tshell, tshell_s
      real(double), dimension(:), pointer :: emt1, emt2, emt3, emt4,   &
                                             emt5, emt6
      real(double) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5,     &
                      elemnt6, elemnt7, contri1, contri2, contri3,     &
                      contri4, contri5, contri6, contri7
      logical :: vsh, nucde, smssh, yes
      character :: cnum*11, ck*2
      integer, dimension(max_orbitals) :: ia_s
      integer :: ka, iopar, incor, ic, lcnum, itjpoc, ir, ia, ib, i, j
      integer :: l, loc, ncontr, lab, ii
!-----------------------------------------------
!
! DOIT: IF DOIT=1 angular coefficients will be stored after creation
! DINT1 contain the density
! DINT2 contain the uncorrected NMS parameter: K^1_NMS
! DINT3 contain the expect. value <r>
! DINT4 contain the expect. value <r2>
! DINT5 contain the expect. value <r-1>
! DINT6 contain the expect. value <r-2>
! DINT7 contain the sum of NMS parameters: K^1_NMS+K^2_NMS+K^3_NMS
!
!
!   Set the rank (zero) and parity (even) for the one-particle
!   coefficients
!
      ka = 0
      iopar = 1
      incor = 1
!
!   Allocate storage for the arrays in BUFFER
!
      call alcbuf (1)
!
!   Sweep through the Hamiltonian matrix to determine the
!   sms parameter
!
     do 13 ic = 1,csf_count
!
!   Output IC on the screen to show how far the calculation has preceede
!
        call convrt (ic,cnum,lcnum)
        if (mod(ic,100).eq.0) then
          print *, 'Column '//cnum(1:lcnum)//' complete;'
        end if
!
        itjpoc = itjpo (ic)

        do 12 ir = ic,csf_count
!
!   Matrix elements are diagonal in J
!
          if (itjpo(ir) .eq. itjpoc) then
!
!   Initialise the accumulator
!
            elemnt1 = 0.0d00
            elemnt2 = 0.0d00
            elemnt3 = 0.0d00
            elemnt4 = 0.0d00
            elemnt5 = 0.0d00
            elemnt6 = 0.0d00
            elemnt7 = 0.0d00
!
!   Call the MCT package to compute T coefficients
!
           call onescalar(ic,ir,ia,ib,tshell)
!GG            CALL TNSRJJ (KA,IOPAR,IC,IR,IA,IB,TSHELL)
            if (ia .ne. 0) then
              if (ia .eq. ib) then
                    ncontr = 0
                do 8 ia = 1,orbital_count
                  if (abs (tshell(ia)) .gt. cutoff) then
                    ncontr = ncontr + 1
                    tshell_s(ncontr) = tshell(ia)
                    ia_s(ncontr) = ia
                    elemnt1 = elemnt1 + dint1(ia,ia)*tshell(ia)
                    elemnt2 = elemnt2 + dint2(ia,ia)*tshell(ia)
                    elemnt3 = elemnt3 + dint3(ia,ia)*tshell(ia)
                    elemnt4 = elemnt4 + dint4(ia,ia)*tshell(ia)
                    elemnt5 = elemnt5 + dint5(ia,ia)*tshell(ia)
                    elemnt6 = elemnt6 + dint6(ia,ia)*tshell(ia)
                    elemnt7 = elemnt7 + dint7(ia,ia)*tshell(ia)
                  endif
    8           continue
                if (doit.eq.1) write(50) ic,ir,ncontr
                do i = 1,ncontr
                  lab = ia_s(i)*(key + 1)
                  if (doit.eq.1) write(50) tshell_s(i),lab
                end do
              else
                if (abs (tshell(1)) .gt. cutoff) then
                  if (orbital_kappa(ia).eq.orbital_kappa(ib)) then
                    if (doit.eq.1) write(50) ic,ir,1
                    lab = ia*key + ib
                    if (doit.eq.1) write(50) tshell(1),lab
                    elemnt1 = elemnt1 + dint1(ia,ib)*tshell(1)
                    elemnt2 = elemnt2 + dint2(ia,ib)*tshell(1)
                    elemnt3 = elemnt3 + dint3(ia,ib)*tshell(1)
                    elemnt4 = elemnt4 + dint4(ia,ib)*tshell(1)
                    elemnt5 = elemnt5 + dint5(ia,ib)*tshell(1)
                    elemnt6 = elemnt6 + dint6(ia,ib)*tshell(1)
                    elemnt7 = elemnt7 + dint7(ia,ib)*tshell(1)
                  endif
                endif
              endif
            endif
            do 9 j = 1,nvec
              loc = (j-1)*csf_count
              contri1 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt1
              contri2 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt2
              contri3 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt3
              contri4 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt4
              contri5 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt5
              contri6 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt6
              contri7 = csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)*elemnt7
              if (ir.ne.ic) then
                contri1 = 2.0d00 * contri1
                contri2 = 2.0d00 * contri2
                contri3 = 2.0d00 * contri3
                contri4 = 2.0d00 * contri4
                contri5 = 2.0d00 * contri5
                contri6 = 2.0d00 * contri6
                contri7 = 2.0d00 * contri7
              endif
              dens1(j) = dens1(j) + contri1
              dens2(j) = dens2(j) + contri2
              dens3(j) = dens3(j) + contri3
              dens4(j) = dens4(j) + contri4
              dens5(j) = dens5(j) + contri5
              dens6(j) = dens6(j) + contri6
              dens7(j) = dens7(j) + contri7
    9       continue
          endif
   12   continue
   13 continue

       if (doit.eq.1) write(50) -1
!
! Empty the buffer and close file
      if (doit.eq.1) close(50)
!
!   Deallocate storage for the arrays in BUFFER
      call alcbuf (3)
      return
      end subroutine densnew
