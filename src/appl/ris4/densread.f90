!***********************************************************************
!                                                                      *
subroutine densread(dint1,dint2,dint3,dint4,dint5,dint6,dint7)
!      SUBROUTINE DENSREAD(DINT1,DINT2,DINT3,                           &
!                          DINT4,DINT5,DINT6,                           &
!                          DINT7,DINT1VEC,DENS1VEC)
!                                                                      *
!   IF angular coefficients already exist                              *
!   This routine controls combines the radial and angular parts for the*
!   calculation of the NMS parameter, the electron density at the      *
!   origin and radial expectation values.                              *
!                                                                      *
!   Written by C\'edric Naz\'e                                         *
!                                                                      *
!                                         Last revision: Feb. 2011     *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: keyorb, nnnw, nnnp
      use prnt_c
      use ris_c
      use orb_c
      use eigv_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), dimension(nnnw,nnnw), intent(in) :: dint1, dint2,  &
                                          dint3, dint4, dint5, dint6,  &
                                          dint7
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(nnnw) :: tshell_r
      real(double) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5
      real(double) :: elemnt6, elemnt7
      real(double) :: contri1, contri2, contri3, contri4, contri5
      real(double) :: contri6, contri7
      integer :: icold, irold, ios, ia, ib, ic, ir, i, j, l, loc, lab
      integer :: ncount
!-----------------------------------------------
!
! DINT1 contain the density
! DINT2 contain the uncorrected NMS parameter: K^1_NMS
! DINT3 contain the expect. value <r>
! DINT4 contain the expect. value <r2>
! DINT5 contain the expect. value <r-1>
! DINT6 contain the expect. value <r-2>
! DINT7 contain the sum of NMS parameters: K^1_NMS+K^2_NMS+K^3_NMS
!

      icold = 0
      irold = 0

      rewind (50)
   16 read (50,iostat = ios) ic,ir,ncount

      if (ios .eq. 0) then
!*      print*, 'ic', IC,IR,NCOUNT
!
!   Read successful; decode the labels of I(ab)
!
!      Initialise
        if ((ic.ne.icold).or.(ir.ne.irold)) then
!*            ISPARC = ISPAR (IC)
!*            ITJPOC = ITJPO (IC)
!*            ITJPOR = ITJPO (IR)
!*            IDIFF = ITJPOC - ITJPOR
            icold = ic
            irold = ir
!*        DO I = 1,NELMNT
          elemnt1 = 0.0d00
          elemnt2 = 0.0d00
          elemnt3 = 0.0d00
          elemnt4 = 0.0d00
          elemnt5 = 0.0d00
          elemnt6 = 0.0d00
          elemnt7 = 0.0d00
!*        ENDDO
        endif
        do i= 1,ncount
           read(50) tshell_r(i),lab
           ia = lab/key
           ib = mod(lab,key)
!*           LOC = INDEX(I)
           elemnt1 = elemnt1 + dint1(ia,ib)*tshell_r(i)
           elemnt2 = elemnt2 + dint2(ia,ib)*tshell_r(i)
           elemnt3 = elemnt3 + dint3(ia,ib)*tshell_r(i)
           elemnt4 = elemnt4 + dint4(ia,ib)*tshell_r(i)
           elemnt5 = elemnt5 + dint5(ia,ib)*tshell_r(i)
           elemnt6 = elemnt6 + dint6(ia,ib)*tshell_r(i)
           elemnt7 = elemnt7 + dint7(ia,ib)*tshell_r(i)
        enddo
!
!   Return to the start of the loop
!

!*      ICI = 0
!*      DO 21 I = 1,NELMNT
!*        IRI = IROW(I)
!*         IF (I .GT. IENDC(ICI)) ICI = ICI+1
         do 22 j = 1,nvec
            loc = (j-1)*ncf
            contri1 = evec(ic+loc)*evec(ir+loc)*elemnt1
            contri2 = evec(ic+loc)*evec(ir+loc)*elemnt2
            contri3 = evec(ic+loc)*evec(ir+loc)*elemnt3
            contri4 = evec(ic+loc)*evec(ir+loc)*elemnt4
            contri5 = evec(ic+loc)*evec(ir+loc)*elemnt5
            contri6 = evec(ic+loc)*evec(ir+loc)*elemnt6
            contri7 = evec(ic+loc)*evec(ir+loc)*elemnt7
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
   22    continue
!*   21 CONTINUE

      goto 16
      endif
      return
      end subroutine densread
