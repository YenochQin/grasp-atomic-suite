      module densread_seltz_mpi_i
      interface
      subroutine densread_seltz_mpi (dint1,dint2,dint3,dint4,dint5,dint6, &
         dint7,dint1vec,dens1vec,nrnuc)
      use vast_kind_param, only: double
      use parameter_def, only: nnnw
      use prnt_c, only: nvec
      integer, intent(in) :: nrnuc
      real(double), dimension(nnnw,nnnw), intent(in) :: dint1, dint2,  &
                                          dint3, dint4, dint5, dint6,  &
                                          dint7
      real(double), dimension(nvec,nrnuc), intent(out) :: dens1vec
      real(double), dimension(nnnw,nnnw,nrnuc), intent(in) :: dint1vec
      end subroutine
      end interface
      end module
