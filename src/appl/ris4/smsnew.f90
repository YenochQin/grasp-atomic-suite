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
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: orbital_index_radix, max_orbitals
      use debug_c
      use prnt_c
      use orb_c
      use buffer_c
      use eigv_c
      use ris_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i, only: alcbuf
      use itjpo_i, only: itjpo
      use convrt_i, only: convrt
      use rkco_gg_i, only: rkco_gg
      use cord_i, only: cord
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
      integer, intent(in) :: doit
!-----------------------------------------------
!  E x t e r n a l   F u n c t i o n s
!-----------------------------------------------

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=11) :: cnum
      character(len=4)  :: jlbl,labj,labp
      character(len=2)  :: ck
      logical      :: getyn,firstt,vsh,nucde,smssh,yes
      logical      :: lfordr,ltrans,lvp,lse,lnms,lsms
      real(real64) :: contri, contrik1, vcoeff
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
      do 13 ic = 1,csf_count
!
!   Output IC on the screen to show how far the calculation has preceede
!
        call convrt (ic,cnum,lcnum)
        if (mod(ic,100)==0) then
          print *, "Column "//cnum(1:lcnum)//" complete;"
        end if
!
        itjpoc = itjpo (ic)
        do 12 ir = ic,csf_count
!
!   Call the MCP package to generate V coefficients; ac and bd
!   are the density pairs
!
!   Initialize
!
!
!   Matrix elements are diagonal in J
!
          if (itjpo(ir) == itjpoc) then
            angular_coefficient_count = 0
!GG            CALL RKCO (IC,IR,COR,CORD,INCOR)
            call rkco_gg (ic, ir, cord, incor, 1)
!
            do 11 i = 1,angular_coefficient_count
              vcoeff = angular_coefficients(i)
              if (abs (vcoeff) > cutoff) then
                iia = angular_coefficient_labels(1,i)
                iib = angular_coefficient_labels(2,i)
                iic = angular_coefficient_labels(3,i)
                iid = angular_coefficient_labels(4,i)
                k   = angular_coefficient_labels(5,i)
!
!   Only K = 1   LABEL(5,I) .EQ. 1
!
                if (angular_coefficient_labels(5,i) == 1) then
                  if (ldbpa(2)) then
                    write (99,309) k,ic,ir,                            &
                    orbital_principal_n(iia),orbital_symmetry_labels(iia),orbital_principal_n(iib), &
                       & orbital_symmetry_labels(iib),                   &
                    orbital_principal_n(iic),orbital_symmetry_labels(iic),orbital_principal_n(iid), &
                       & orbital_symmetry_labels(iid),vcoeff
                  end if
                  if(doit==1) write(51) ic,ir
!** Storage sequence
                  lab  = ((iia*key + iic)*key+iib)*key+iid
                  if(doit==1) then
                    write(51) vcoeff,lab
                  end if
!***
                  do 10 j = 1,nvec
                    loc = (j-1)*csf_count
                    contrik1 = - csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)             &
                           * vcoeff                                    &
                           * vint (angular_coefficient_labels(1,i),angular_coefficient_labels(3,i))              &
                           * vint (angular_coefficient_labels(2,i),angular_coefficient_labels(4,i))
                    contri = - csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)               &
                           * vcoeff                                    &
                           * ( vint2(angular_coefficient_labels(1,i),angular_coefficient_labels(3,i))            &
                           * vint (angular_coefficient_labels(2,i),angular_coefficient_labels(4,i))              &
                           + vint2(angular_coefficient_labels(2,i),angular_coefficient_labels(4,i))              &
                           * vint (angular_coefficient_labels(1,i),angular_coefficient_labels(3,i)) )/2.0_real64
                    if (ir/=ic) then
                      contri = 2.0_real64 * contri
                      contrik1 = 2.0_real64 * contrik1
                    end if
                    smsc1(j) = smsc1(j) + contrik1
                    smsc2(j) = smsc2(j) + contri
   10             continue
                end if
              end if
   11       continue
          end if
   12   continue
   13 continue
      if(doit==1) write(51) -1
!
! Empty the buffer and close file
      if(doit==1) close(51)
!
!   Deallocate storage for the arrays in BUFFER
!
      call alcbuf (3)
      return
  309 format (" V^[(",1i2,")]_[",1i3,",",1i3,"]",                      &
         " (",1i2,1a2,",",1i2,1a2,";",                                 &
              1i2,1a2,",",1i2,1a2,") = ",1pd19.12)
!
      end subroutine smsnew
