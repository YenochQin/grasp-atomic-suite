!***********************************************************************
!                                                                      *
      subroutine mcpin (name,ik,ntestg,inpci)
      use biorthogonal_subshell_layout, only: biorthogonal_angular_unit
      use calculation_file_units, only: mixing_input_unit
!                                                                      *
!   This routine controls the computation  and storage of the values   *
!   and all indices of the angular coefficients                        *
!                                                                      *
!                                                                      *
!                   T  (ab)                                            *
!                    rs                                                *
!                                                                      *
!   k is the multipolarity of a two-particle Coulomb integral. a, b,   *
!   are orbital sequence numbers.  r and s are configuration           *
!   state function indices.                                            *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, TNSRJJ                         *
!               [GENMCP]: QSORT.                                       *
!                                                                      *
!   Written by Per Jonsson                Last revision: JUne 1996     *
!   Bug corrected 2005-10-18                                           *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: orbital_index_radix, max_orbitals
      use array_memory_management, only: alloc, dalloc
      use atomic_calculation_state, only: electron_count
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use first_order_csf_cutoffs, only: iccut
      use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
      use orbital_csf_basis,           only: csf_count, orbital_count
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
      use biorthogonal_transition_selection, only: itjqj1, itjqj2, jqj1, jqj2, ntrans
      use biorthogonal_orbital_order, only: nordff, nordii
      use mcp_coefficient_storage, only: cnn, intgrl, intptr, jann, jbnn, ncoeff, nintg
      use biorthogonal_subshell_layout,         only: kamax,nshlii,nshlff,nshlppii,nshlppff,&
                                 ninii,ninff,nshlpii,nshlpff,nakinvii, &
                                 nakinvff,nlmax
      use biorthogonal_subshell_mapping,        only: nshlp, nshlpp
      use biorthogonal_orbital_rotations,         only: cici, cfci
      use transition_basis_blocks, only: ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
          nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE tnsrjj_I
!      USE cor_I
!      USE cord_I
!      USE itjpo_I
      use citrag_i, only: citrag
      use getmix_i, only: getmix
      implicit none (type, external)
      integer :: biorthogonal_mixing_output_unit
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=24), intent(in) :: name
      integer, intent(in) :: ik, ntestg, inpci
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!   Set the encoding key and its square
      integer, parameter :: orbital_label_radix = orbital_index_radix, keysq = orbital_label_radix*orbital_label_radix, nvmax=100
      real(real64), parameter :: cutoff = 1.0e-10_real64

      logical ::  f0int,lincr,is_restart,comp,avail
      character(len=20) :: number_text
      character(len=2)  :: ck
!
      real(real64) :: shell_coefficients(max_orbitals)
      real(real64) ::sc(nvmax)
      integer :: nakinv(max_orbitals)
      integer :: llistt(nlmax), nshl(nlmax), ninl(nlmax)
      real(real64) :: cirot(20*nlmax*nlmax)
      real(real64) :: evsc(10000)
!
      real(real64), pointer :: scr(:), ciout(:)
      integer :: ntest, ntestl, i, j, iblk, ncfd, nwd, kamaxd,  &
                 l, iioff, nvecsize
!
      ntestl = 00
      ntest = max(ntestg,ntestl)
      ntest = 00
!
!   Set up data for the initial and final state case respectively
!
      if (ik==1) then
        do i = 1,nlmax
          nshl(i) = nshlii(i)
          ninl(i) = ninii(i)
        end do

        do i = 1,max_orbitals
          nakinv(i) = nakinvii(i)
        end do

        do j = 1,nlmax
          do i = 1,nlmax
            nshlp(i,j) = nshlpii(i,j)
          end do
        end do

        do j = 1,max_orbitals
          do i = 1,nlmax
            nshlpp(i,j) = nshlppii(i,j)
          end do
        end do

        do i = 1,20*nlmax*nlmax
          cirot(i) = cici(i)
        end do

      else if (ik==2) then

        do i = 1,nlmax
          nshl(i) = nshlff(i)
          ninl(i) = ninff(i)
        end do

        do i = 1,max_orbitals
          nakinv(i) = nakinvff(i)
        end do

        do j = 1,nlmax
          do i = 1,nlmax
            nshlp(i,j) = nshlpff(i,j)
          end do
        end do

        do j = 1,max_orbitals
          do i = 1,nlmax
            nshlpp(i,j) = nshlppff(i,j)
          end do
        end do

        do i = 1,20*nlmax*nlmax
          cirot(i) = cfci(i)
        end do
      end if

      rewind (biorthogonal_angular_unit)
      do 1000 iblk = 1, nblock
!
!   Read the CI vectors
!
      call getmix(name,inpci,iblk)
!
!   Allocate memory for use in CITRAG
!
      call alloc (scr, csf_count*nvec, "scr", "MCPIN")
      call alloc (ciout, csf_count*nvec, "CIOUT", "MCPIN")
!
!   if data available on file read the datafile. Else sort the
!   MCP data into inegral based lists.
!
      read(biorthogonal_angular_unit) ncfd,nwd,kamaxd
      do l = 1,kamax

!************
!
!. Offset for given L in shell matrices
!
          if (l==1) then
            iioff = 1
          else
            iioff = iioff + nshl(l-1)** 2
          end if
!  Corrected PER J

!**************

          read(biorthogonal_angular_unit) nintg,ncoeff
          if(ncoeff*nintg==0) cycle
!
!   Allocate memory. Note that in this case we can,
!   since the data is sorted, allocate less
!   memory than is done in qqsort.
!
          call alloc (jann,ncoeff, "JANN", "MCPIN")
          call alloc (jbnn,ncoeff, "JBNN", "MCPIN")
          call alloc (intgrl,nintg, "INTGRL", "MCPIN")
          call alloc (cnn,ncoeff, "CNN", "MCPIN")
          call alloc (intptr,nintg, "INTPTR", "MCPIN")

          do i = 1,nintg
            read(biorthogonal_angular_unit) intgrl(i),intptr(i)
          end do
          do i = 1,ncoeff
            read(biorthogonal_angular_unit) cnn(i),jann(i),jbnn(i)
          end do
!*
!*. Offset for given L in shell matrices
!*
!        IF (L.EQ.1) THEN
!          IIOFF = 1
!        ELSE
!          IIOFF = IIOFF + NSHL(L-1)** 2
!        END IF
!   BUG THE OFFSET SHOULD BE BEFORE THE CYCLE STATEMENT!!
!
! Transform for given L
!
        if (nshl(l)>0)                                              &
     &  call citrag(csf_mixing_coefficients(:),csf_count,nvec,l,nshl(l),                    &
     &              cirot(iioff),ninl(l),ntestl,ciout, scr )
!
!   Deallocate storage for this kappa.
!
        call dalloc (jann, "JANN", "MCPIN")
        call dalloc (jbnn, "JBNN", "MCPIN")
        call dalloc (intgrl, "INTGRL", "MCPIN")
        call dalloc (cnn, "CNN", "MCPIN")
        call dalloc (intptr, "INTPRT", "MCPIN")

      end do

      call dalloc (scr, "SCR", "MCPIN")
      call dalloc (ciout, "CIOUT", "MCPIN")
!
!   Write the rotated CI vectors on file
!
      if(iblk==1) then
        j = index(name," ")
        if (inpci==0) then
          open (newunit=biorthogonal_mixing_output_unit,file=name(1:j-1)//".cbm",form="UNFORMATTED",  &
     &    status="UNKNOWN")
        else
          open (newunit=biorthogonal_mixing_output_unit,file=name(1:j-1)//".bm",form="UNFORMATTED",   &
     &    status="UNKNOWN")
        end if

        write(biorthogonal_mixing_output_unit) "G92MIX"
        write(biorthogonal_mixing_output_unit) electron_count,ncftot,orbital_count,nvectot,nvecsize,nblock
      end if

      write(biorthogonal_mixing_output_unit) iblk,csf_count,nvec,iatjpo(1),iaspar(1)
      write(biorthogonal_mixing_output_unit) (ivec(i),i = 1,nvec)
!     WRITE(31) (IATJPO(I),IASPAR(I),I = 1,NVEC)
      write(biorthogonal_mixing_output_unit) hamiltonian_energy_shift,(shifted_eigenvalues(i),i = 1,nvec)
      write(biorthogonal_mixing_output_unit) ((csf_mixing_coefficients(i+(j-1)*csf_count),i = 1,csf_count),j = 1,nvec)

      call dalloc (shifted_eigenvalues, "EVAL", "MCPIN")
      call dalloc (csf_mixing_coefficients, "EVEC", "MCPIN")
      call dalloc (ivec, "IVEC", "MCPIN")
      call dalloc (iatjpo, "IATJPO", "MCPIN")
      call dalloc (iaspar, "IASPAR", "MCPIN")
 1000 continue
      close (biorthogonal_angular_unit)
!
!   Close the  mixing  file
!
      close (mixing_input_unit)
      close (biorthogonal_mixing_output_unit)
!
      return
      end subroutine mcpin
