!***********************************************************************
!                                                                      *
      subroutine ti1tv(ciin,ncsf,nciv,i,l,t,nshl,ciout,ntestg)
!                                                                      *
!   Calculate the action of the  operator                              *
!   Sum( j) T(jl) E(jl,il) on a set of vectors                         *
!                                                                      *
!   only j<i terms ! included
!                                                                      *
!   Modification of the TI1TV program for MCHF                         *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!   =====                                                              *
!   Input                                                              *
!   =====                                                              *
!   CIIN  : Input CI vectors                                           *
!   NCSCF : Length of CI expansion                                     *
!   NCIV  : Number of CI vectors                                       *
!   I     : Shell number                                               *
!   L     : L value of shells of this excitation                       *
!   NSHL  : Number of shells with this L                               *
!   T     : The integrals T(*,il)                                      *
!                                                                      *
!   ======                                                             *
!   Output                                                             *
!   ======                                                             *
!                                                                      *
!   CIOUT : List of output CI vectors                                  *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use parameter_def,   only: keyorb
      use mcpdata_c
      use sbdat1_c
      use orbord_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setvec_i
      use wrtmat_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ncsf, nciv, nshl, i, l, ntestg
      real(double), dimension(ncsf, nciv), intent(in) :: ciin
      real(double), dimension(nshl), intent(in) :: t
      real(double), dimension(ncsf, nciv), intent(out) :: ciout
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: key=keyorb
      integer :: ntestl, ntest, nfound, k, ia, ib, ir, ira, irb, ifirst, &
                 j, ileft, iright, ivec, ielmnt
      real(double) :: racah, x
!
      ntestl = 00
      ntest = max(ntestl,ntestg)
      ntest = 000
!
      if(ntest.ge.10) write(6,*) ' Entering TI1TV'

      call setvec(ciout,0.0d0,nciv*ncsf)
!
!.  Obtain address of first coupling coefficient for h(jl,il) :  IFIRST
!.  Obtain number of        coupling coefficient for h(jl,il) :  NFOUND
!
!        NFOUND : Number of coefficients obtained
!
      nfound = 0
      do k = 1,nintg
        ia = intgrl(k)/key
        ib = mod(intgrl(k),key)
!
!   Determine the ordernumber of the orbitals I, IA and IB within the
!   shell with symmetry L
!
        ir = nshlpp(l,nshlp(l,i))
        ira = nshlpp(l,ia)
        irb = nshlpp(l,ib)

        if (ir.eq.irb.and.ira.lt.irb) then
          if (k.eq.1) then
            ifirst = 1
          else
            ifirst = intptr(k-1) + 1
          endif
          nfound = intptr(k) - ifirst + 1

   15     do ielmnt = 1, nfound
            racah = cnn(ifirst-1+ielmnt)
            j = ira
            x = t(j)*racah
!ww            IF (NORDII.EQ.0.AND.NORDFF.EQ.0) THEN
!
!  Expression for normal orbital ordering
!
              ileft = jann(ifirst-1+ielmnt)
              iright = jbnn(ifirst-1+ielmnt)
!ww            ELSEIF (NORDII.EQ.1.AND.NORDFF.EQ.1) THEN
!
!  Expression for reverse orbital ordering
!
!ww              IRIGHT = JANN(IFIRST-1+IELMNT)
!ww              ILEFT = JBNN(IFIRST-1+IELMNT)
!ww            ELSE
!ww              WRITE(*,*) 'SOMETHING WRONG'
!ww              STOP
!ww            ENDIF
            do ivec =1, nciv
              ciout(ileft,ivec) =                                       &
     &        ciout(ileft,ivec) + x*ciin(iright,ivec)
            enddo
          enddo
        endif
      enddo
!
      if(ntest.ge.100) then
        write(6,*)
        write(6,*) ' Input and output vectors from TI1TV I,L',i,l
        call wrtmat(ciin,ncsf,nciv,ncsf,nciv)
        write(6,*)
        call wrtmat(ciout,ncsf,nciv,ncsf,nciv)
        write(6,*)
      end if
!
      if(ntest.ge.10) write(6,*) ' LEAVING TI1TV'

      return
      end
