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
      use grasp_parameters,    only: orbital_index_radix, max_orbitals
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_principal_n, orbital_symmetry_labels
      use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
          angular_coefficient_labels, angular_coefficients
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use ris_c, only: dens1, dens2, dens3, dens4, dens5, dens6, dens7, smsc1, smsc2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use angular_coefficient_buffer, only: alcbuf
      use csf_quantum_queries, only: itjpo
      use integer_format, only: convrt
      use two_electron_angular_coefficients, only: rkco_gg
      use closed_shell_mcp, only: cord
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
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=11) :: number_text
      character(len=4)  :: jlbl,labj,labp
      character(len=2)  :: ck
      logical      :: getyn,firstt,vsh,nucde,smssh,yes
      logical      :: lfordr,ltrans,lvp,lse,lnms,lsms
      real(real64) :: contri, contrik1, vcoeff
      integer      :: i,j,k,lab,loc,csf_index,ir,itjpoc,iia,iib,iic,iid,include_closed_shell_terms
      integer      :: text_length
!-----------------------------------------------
!
      include_closed_shell_terms = 1
!
!   Allocate storage for the arrays in BUFFER
!
      call alcbuf (1)
!
!   Sweep through the Hamiltonian matrix to determine the
!   sms parameter
!
      do 13 csf_index = 1,csf_count
!
!   Output IC on the screen to show how far the calculation has preceede
!
        call convrt (csf_index,number_text,text_length)
        if (mod(csf_index,100)==0) then
          print *, "Column "//number_text(1:text_length)//" complete;"
        end if
!
        itjpoc = itjpo (csf_index)
        do 12 ir = csf_index,csf_count
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
            call rkco_gg (csf_index, ir, cord, include_closed_shell_terms, 1)
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
                    write (99,'(" V^[(",1i2,")]_[",1i3,",",1i3,"]"," (",1i2,1a2,",",1i2,1a2,";",1i2,1a2,",",1i2,&
                       &1a2,") = ",1pd19.12)') k,csf_index,ir,                            &
                    orbital_principal_n(iia),orbital_symmetry_labels(iia),orbital_principal_n(iib), &
                       & orbital_symmetry_labels(iib),                   &
                    orbital_principal_n(iic),orbital_symmetry_labels(iic),orbital_principal_n(iid), &
                       & orbital_symmetry_labels(iid),vcoeff
                  end if
                  if(doit==1) write(51) csf_index,ir
!** Storage sequence
                  lab  = ((iia*orbital_label_radix + iic)*orbital_label_radix+iib)*orbital_label_radix+iid
                  if(doit==1) then
                    write(51) vcoeff,lab
                  end if
!***
                  do 10 j = 1,nvec
                    loc = (j-1)*csf_count
                    contrik1 = - csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)             &
                           * vcoeff                                    &
                           * vint (angular_coefficient_labels(1,i),angular_coefficient_labels(3,i))              &
                           * vint (angular_coefficient_labels(2,i),angular_coefficient_labels(4,i))
                    contri = - csf_mixing_coefficients(csf_index+loc)*csf_mixing_coefficients(ir+loc)               &
                           * vcoeff                                    &
                           * ( vint2(angular_coefficient_labels(1,i),angular_coefficient_labels(3,i))            &
                           * vint (angular_coefficient_labels(2,i),angular_coefficient_labels(4,i))              &
                           + vint2(angular_coefficient_labels(2,i),angular_coefficient_labels(4,i))              &
                           * vint (angular_coefficient_labels(1,i),angular_coefficient_labels(3,i)) )/2.0_real64
                    if (ir/=csf_index) then
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
!
      end subroutine smsnew
