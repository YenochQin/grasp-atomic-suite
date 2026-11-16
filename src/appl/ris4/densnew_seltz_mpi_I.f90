      module densnew_seltz_mpi_i
      interface
      subroutine densnew_seltz_mpi (doit,dint1,dint2,dint3,dint4,dint5,dint6, &
         dint7,dint1vec,dens1vec,nrnuc)
      use vast_kind_param, only: double
      use parameter_def, only: max_orbitals
      use prnt_c, only: nvec
      integer, intent(in) :: doit, nrnuc
      real(double), dimension(max_orbitals,max_orbitals), intent(in) :: dint1, dint2, &
                                                        dint3, dint4, &
                                                        dint5, dint6, &
                                                        dint7
      real(double), dimension(nvec,nrnuc), intent(out) :: dens1vec
      real(double), dimension(max_orbitals,max_orbitals,nrnuc), intent(in) :: dint1vec
      end subroutine
      end interface
      end module
