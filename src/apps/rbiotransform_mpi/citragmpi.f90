!***********************************************************************
!                                                                      *
!                                                                      *
      subroutine citrag(input_mixing_coefficients,csf_count,level_count,l,subshell_count,t,nin,ntestg, &
         & transformed_mixing_coefficients,scr)
!                                                                      *
!   Calculate the action of the  operator                              *
!                                                                      *
!    T(L) = T(nl) T(n-1 L ) .... T(1 l)                                *
!                                                                      *
!   where                                                              *
!                                                                      *
!    T(i l ) = ( sum(k=0,2l+1) t(i l ) ** k ) (t_ilil ** E(il,il))     *
!    t(i l ) = sum(j.lt.i) t(j,i) E(lj,li)                             *
!    E(lj,li) total symmetric shell excitation operator                *
!                                                                      *
!   Modification of the CITRA routine for MCHF                         *
!                                                                      *
!   =====                                                              *
!   Input                                                              *
!   =====                                                              *
!   CIIN      : Input CI vectors ( Destroyed in the process )          *
!   NCSCF     : Length of CI expansion                                 *
!   NCIV      : Number of CI vectors                                   *
!   L         : L value of shells in excitations                       *
!   NSHL      : Number of shells with this L                           *
!   T         : List of coefficients                                   *
!   NIN       : Number of inactive shells                              *
!                                                                      *
!   ======                                                             *
!   Output                                                             *
!   ======                                                             *
!                                                                      *
!   CIOUT : List of output CI vectors                                  *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use biorthogonal_subshell_mapping, only: max_orbitals, nlmax, nshlp, nshlpp
      use orbital_csf_basis, only: max_orbitals, orbital_kappa
      use mpi_state, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_sum
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use wrtmat_i, only: wrtmat
      use scalve_i, only: scalve
!      USE copvec_I
      use tiinig_i, only: tiinig
      use ti1tv_i, only: ti1tv
      use vecsum_i, only: vecsum
      use copvec_i, only: copvec
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: csf_count
      integer, intent(in) :: level_count
      integer, intent(in) :: l
      integer, intent(in) :: subshell_count
      integer, intent(in) :: nin
      integer, intent(in) :: ntestg
      real(real64)  :: input_mixing_coefficients(csf_count,level_count)
      real(real64), intent(in) :: t(subshell_count,subshell_count)
      real(real64)  :: transformed_mixing_coefficients(csf_count,level_count)
      real(real64)  :: scr(csf_count,level_count)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntestl, ntest, iin, ipot, i, n
      real(real64) :: factor, tii, xnfaci
      real(real64)  :: scrtmp(csf_count,level_count)
!-----------------------------------------------
!

      ntestl = 0
      ntest = max(ntestg,ntestl)
      ntest = 0

!
      if (ntest >= 10) then
         write (6, *)
         write (6, *) " ***************"
         write (6, *) " Entering CITRAG"
         write (6, *) " ***************"
         write (6, *)
      end if
      if (ntest >= 100) then
         write (6, *) " Input CI vectors "
         call wrtmat (input_mixing_coefficients, csf_count, level_count, csf_count, level_count)
         write (6, *) " Transformation matrix T"
         call wrtmat (t, subshell_count, subshell_count, subshell_count, subshell_count)
      end if
!
!. Factor from inactive shells
!
      if (nin /= 0) then
         factor = 1.0_real64
         do iin = 1, nin
            factor = factor*t(iin,iin)
         end do
!
!       IPOT = 2*(2*L+1)  (number of m_lm_s. This should be replaced
!                          by (2j+1) corresponding to L)
!
         ipot = 2*iabs(orbital_kappa(nshlp(l,iin)))
         factor = factor**ipot
         call scalve (input_mixing_coefficients, factor, level_count*csf_count)
      end if
      if (nin == subshell_count) call copvec (input_mixing_coefficients, transformed_mixing_coefficients, level_count*csf_count)
!
      do i = nin + 1, subshell_count
         if (ntest >= 100) write (6, *) " Loop I,L = ", i, l
!
!. The diagonal contribution
!
         tii = t(i,i)
         call tiinig (input_mixing_coefficients, csf_count, level_count, i, l, tii, transformed_mixing_coefficients, ntestg)
!        IF (LWORK2.LT.NCIV*NCSF) THEN
!          WRITE(*,*) 'In CITRAG: Dimension of LWORK2 must be',
!     &               'increased to at least',NCIV*NCSF
!        ENDIF
         call copvec (transformed_mixing_coefficients, scr, level_count*csf_count)
!
!.  Off diagonal contributions
!
         xnfaci = 1.0_real64
         do n = 1, 2*iabs(orbital_kappa(nshlp(l,i)))
            if (ntest >= 100) write (6, *) " Loop N = ", n
!
!   T ** (N-1) is supposed to be in SCR, copy to CIIN
!   and apply S
!
            call copvec (scr, input_mixing_coefficients, level_count*csf_count)
            call ti1tv (input_mixing_coefficients,csf_count,level_count,i,l,t(1,i),subshell_count,scrtmp,ntestg)
!            CALL MPI_ALLREDUCE(SCRtmp(1,1),SCR(1,1),NCIV*NCSF,         &
            call mpi_allreduce(scrtmp,scr,level_count*csf_count,         &
            mpi_double_precision, mpi_sum, mpi_comm_world, mpi_error_code)
            xnfaci = xnfaci/float(n)
            call vecsum (transformed_mixing_coefficients, transformed_mixing_coefficients, scr, 1.0_real64, xnfaci, &
               & level_count*csf_count)
         end do
         call copvec (transformed_mixing_coefficients, input_mixing_coefficients, level_count*csf_count)
!
      end do
!
      if (ntest >= 100) then
         write (6, *) " Output CI vectors L = ", l
         call wrtmat (transformed_mixing_coefficients, csf_count, level_count, csf_count, level_count)
      end if
!
      return
      end subroutine citrag
