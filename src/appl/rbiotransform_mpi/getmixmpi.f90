!***********************************************************************
!                                                                      *
      subroutine getmix(name, inpci, iblk)
!                                                                      *
!   Open, check, load data from and close the rscf.mix file.           *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, LENGTH, OPENFL.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 25 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man, only: alloc, byte, double
      use def_c, only: double, electron_count
      use orb_c, only: csf_count, orbital_count, csf_subshell_occupations
      use eigv_c, only: double, hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use prnt_c, only: byte, double, ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use syma_c, only: double, iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use blk_c, only: double, ncfblk, nevinblk, ncfinblk, two_j, nblock, idxblk, nelectot, ncftot, nwtot, nvectot, &
          nvecsiztot, nblock1, nblocki, nblockf, ncfi, ncff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: inpci
      integer, intent(in) :: iblk
      character, intent(in) :: name*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ik, j, ios, nb, iatjp, iaspa, i
      character :: g92mix*6
!-----------------------------------------------

      ik = 30
      if (iblk == 1) then
         j = index(name," ")
         if (inpci == 0) then
            open(unit=ik, file=name(1:j-1)//".cm", form="UNFORMATTED", status=&
               "OLD")
         else
            open(unit=ik, file=name(1:j-1)//".m", form="UNFORMATTED", status=&
               "OLD")
         end if

         read (ik, iostat=ios) g92mix
         if (ios/=0 .or. g92mix/="G92MIX") then
            write (6, *)                                              &
                    "File", ik, "Not a GRASP MIXing Coefficients File;"
            close(ik)
            stop
         end if
!
         read (ik) nelectot, ncftot, nwtot, nvectot, nvecsiztot, nblock1
         if (electron_count/=nelectot .or. orbital_count/=nwtot) then
!     :    (NCF .NE. NCFT) .OR.
            write (6, *) "File", ik, "is not an"
            write (6, *) " appropriate to Coefficients file"
            close(ik)
            stop
         end if
      end if
!
!   Load data from the  rscf.mix  file
!
      write (6, *) "Loading MIXing Coefficients File ..."
!
      read (ik) nb, csf_count, nvec, iatjp, iaspa
      call alloc (shifted_eigenvalues, nvec, "EVAL", "GETMIX")
      call alloc (csf_mixing_coefficients, csf_count*nvec, "EVEC", "GETMIX")
      call alloc (ivec, nvec, "IVEC", "GETMIX")
      call alloc (iatjpo, nvec, "IATJPO", "GETMIX")
      call alloc (iaspar, nvec, "IASPAR", "GETMIX")
!
!   These arrays are deallocated in mcp
!
      read (ik) (ivec(i),i=1,nvec)
!     READ (IK) (IATJPO(I),IASPAR(I),I = 1,NVEC)
      iatjpo(:nvec) = iatjp
      iaspar(:nvec) = iaspa
      read (ik) hamiltonian_energy_shift, (shifted_eigenvalues(i),i=1,nvec)
      read (ik) ((csf_mixing_coefficients(i + (j - 1)*csf_count),i=1,csf_count),j=1,nvec)
!
      write (6, *) " ... load complete;"
!
!   Close the  rscf.mix  file
!
!     CLOSE (IK)
!
      return
      end subroutine getmix
