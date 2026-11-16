!***********************************************************************
!                                                                      *
      subroutine smsnew(doit,vint,vint2)
!                                                                      *
!   This routine controls the main sequence of routine calls for the   *
!   calculation  of the  sms parameter, the electron density at the    *
!   origin.
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, GETYN                         *
!                        ITJPO,                                        *
!               [SMS92]: RINTISO, RINTDENS, VINTI                      *
!               [LIBRANG]: RKCO_GG                                     *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!                                         Last revision: 10 Nov 1995   *
!                                                                      *
!   Modified by C. Naz\'e  Mai. 2012                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: keyorb, nnnw
      use debug_c
      use prnt_c
      use orb_c
      use buffer_c
      use eigv_c
      use ris_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i
      use itjpo_i
      use convrt_i, only: convrt
      use rkco_gg_i, only: rkco_gg
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), dimension(nnnw,nnnw), intent(in) :: vint, vint2
      integer, intent(in) :: doit
!-----------------------------------------------
!  E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      external     :: cord
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character*11 :: cnum
      character*4  :: jlbl,labj,labp
      character*2  :: ck
      logical      :: getyn,firstt,vsh,nucde,smssh,yes
      logical      :: lfordr,ltrans,lvp,lse,lnms,lsms
      real(double) :: contri, contrik1, vcoeff
      integer      :: i,j,k,lab,loc,ic,ir,itjpoc,iia,iib,iic,iid,incor
      integer      :: lcnum
!-----------------------------------------------
!
      incor = 1
!
!   Allocate storage for the arrays in BUFFER
!
      call alcbuf (1)
!
!   Sweep through the Hamiltonian matrix to determine the
!   sms parameter
!
      do 13 ic = 1,ncf
!
!   Output IC on the screen to show how far the calculation has preceede
!
        call convrt (ic,cnum,lcnum)
        if (mod(ic,100).eq.0) then
          print *, 'Column '//cnum(1:lcnum)//' complete;'
        end if
!
        itjpoc = itjpo (ic)
        do 12 ir = ic,ncf
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs
!
!   Initialize
!
!
!   Matrix elements are diagonal in J
!
          if (itjpo(ir) .eq. itjpoc) then
            nvcoef = 0
!GG            CALL RKCO (IC,IR,COR,CORD,INCOR)
            call rkco_gg (ic, ir, cord, incor, 1)
!
            do 11 i = 1,nvcoef
              vcoeff = coeff(i)
              if (abs (vcoeff) .gt. cutoff) then
                iia = label(1,i)
                iib = label(2,i)
                iic = label(3,i)
                iid = label(4,i)
                k   = label(5,i)
!
!   Only K = 1   LABEL(5,I) .EQ. 1
!
                if (label(5,i) .eq. 1) then
                  if (ldbpa(2)) then
                    write (99,309) k,ic,ir,                            &
                    np(iia),nh(iia),np(iib),nh(iib),                   &
                    np(iic),nh(iic),np(iid),nh(iid),vcoeff
                  endif
                  if(doit.eq.1) write(51) ic,ir
!** Storage sequence
                  lab  = ((iia*key + iic)*key+iib)*key+iid
                  if(doit.eq.1) then
                    write(51) vcoeff,lab
                  endif
!***
                  do 10 j = 1,nvec
                    loc = (j-1)*ncf
                    contrik1 = - evec(ic+loc)*evec(ir+loc)             &
                           * vcoeff                                    &
                           * vint (label(1,i),label(3,i))              &
                           * vint (label(2,i),label(4,i))
                    contri = - evec(ic+loc)*evec(ir+loc)               &
                           * vcoeff                                    &
                           * ( vint2(label(1,i),label(3,i))            &
                           * vint (label(2,i),label(4,i))              &
                           + vint2(label(2,i),label(4,i))              &
                           * vint (label(1,i),label(3,i)) )/2.0d00
                    if (ir.ne.ic) then
                      contri = 2.0d00 * contri
                      contrik1 = 2.0d00 * contrik1
                    endif
                    smsc1(j) = smsc1(j) + contrik1
                    smsc2(j) = smsc2(j) + contri
   10             continue
                endif
              endif
   11       continue
          endif
   12   continue
   13 continue
      if(doit.eq.1) write(51) -1
!
! Empty the buffer and close file
      if(doit.eq.1) close(51)
!
!   Deallocate storage for the arrays in BUFFER
!
      call alcbuf (3)
      return
  309 format (' V^[(',1i2,')]_[',1i3,',',1i3,']',                      &
         ' (',1i2,1a2,',',1i2,1a2,';',                                 &
              1i2,1a2,',',1i2,1a2,') = ',1pd19.12)
!
      end subroutine smsnew
