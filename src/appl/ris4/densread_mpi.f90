!***********************************************************************
!                                                                      *
      subroutine densread_mpi(dint1,dint2,dint3,dint4,dint5,dint6,dint7)
!                                                                      *
!   MPI version of DENSREAD. Each rank scans the full angular file but *
!   only accumulates records whose IC belongs to its strided subset.    *
!                                                                      *
!***********************************************************************
      use vast_kind_param,  only: double
      use parameter_def,    only: orbital_index_radix, max_orbitals
      use debug_c,          only: cutoff
      use prnt_c
      use ris_c
      use orb_c
      use eigv_c
      use mpi_c,            only: mpi_rank, mpi_process_count
      implicit none
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint1, dint2,  &
                                          dint3, dint4, dint5, dint6,  &
                                          dint7
      integer, parameter :: key = orbital_index_radix
      real(double), dimension(max_orbitals) :: tshell_r
      real(double) :: elemnt1, elemnt2, elemnt3, elemnt4, elemnt5
      real(double) :: elemnt6, elemnt7
      real(double) :: contri1, contri2, contri3, contri4, contri5
      real(double) :: contri6, contri7
      real(double) :: evpair, pair_scale
      integer :: ios, ia, ib, ic, ir, i, j, loc, lab, ncount
      logical :: ownrec

      rewind (50)
   16 read (50,iostat = ios) ic,ir,ncount

      if (ios .eq. 0) then
         ownrec = mod(ic-1, mpi_process_count) .eq. mpi_rank
         elemnt1 = 0.0d00
         elemnt2 = 0.0d00
         elemnt3 = 0.0d00
         elemnt4 = 0.0d00
         elemnt5 = 0.0d00
         elemnt6 = 0.0d00
         elemnt7 = 0.0d00

         do i= 1,ncount
            read(50) tshell_r(i),lab
            if (.not. ownrec) cycle
            ia = lab/key
            ib = mod(lab,key)
            elemnt1 = elemnt1 + dint1(ia,ib)*tshell_r(i)
            elemnt2 = elemnt2 + dint2(ia,ib)*tshell_r(i)
            elemnt3 = elemnt3 + dint3(ia,ib)*tshell_r(i)
            elemnt4 = elemnt4 + dint4(ia,ib)*tshell_r(i)
            elemnt5 = elemnt5 + dint5(ia,ib)*tshell_r(i)
            elemnt6 = elemnt6 + dint6(ia,ib)*tshell_r(i)
            elemnt7 = elemnt7 + dint7(ia,ib)*tshell_r(i)
         enddo

         if (ownrec) then
            if (abs(elemnt1) .le. cutoff .and. &
                abs(elemnt2) .le. cutoff .and. &
                abs(elemnt3) .le. cutoff .and. &
                abs(elemnt4) .le. cutoff .and. &
                abs(elemnt5) .le. cutoff .and. &
                abs(elemnt6) .le. cutoff .and. &
                abs(elemnt7) .le. cutoff) goto 16
            if (ir .eq. ic) then
               pair_scale = 1.0d00
            else
               pair_scale = 2.0d00
            endif
            do j = 1,nvec
               loc = (j-1)*csf_count
               evpair = pair_scale * csf_mixing_coefficients(ic+loc)*csf_mixing_coefficients(ir+loc)
               contri1 = evpair*elemnt1
               contri2 = evpair*elemnt2
               contri3 = evpair*elemnt3
               contri4 = evpair*elemnt4
               contri5 = evpair*elemnt5
               contri6 = evpair*elemnt6
               contri7 = evpair*elemnt7
               dens1(j) = dens1(j) + contri1
               dens2(j) = dens2(j) + contri2
               dens3(j) = dens3(j) + contri3
               dens4(j) = dens4(j) + contri4
               dens5(j) = dens5(j) + contri5
               dens6(j) = dens6(j) + contri6
               dens7(j) = dens7(j) + contri7
            end do
         endif

         goto 16
      endif

      return
      end subroutine densread_mpi
