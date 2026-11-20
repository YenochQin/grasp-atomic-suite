!
!***********************************************************************
!                                                                      *
      program jj2lsj_omp
!                                                                      *
!     This MAIN program controls the transformation of atomic states,  *
!     which are given in a jj-coupled CSF basis, into an LS-coupled    *
!     basis. The program requires a jj-coupled basis in standard order *
!     where, if both subshells | n j = l-1/2> and | n j = l+1/2>       *
!     of a given shell (nl) occurs, they always follow successively    *
!     in this order. The LS-coupled basis, moreover, is given in       *
!     the same sequence of shells.                                     *
!                                                                      *
!     All LS-jj transformation coefficients are precalculated and      *
!     'stored' in the modules rabs_lsj_data_1, rabs_lsj_data_2 and     *
!     rabs_lsj_data_3.                                                 *
!                                                                      *
!     Calls:  FACTT, SETISO, JJ2LSJ, starttime, stoptime.              *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                                  May 2011   *
!     VILNIUS                                               May 2017   *
!                                                                      *
!     Modified by G. Gaigalas and C. Cychen                     2021   *
!     Modified by G. Gaigalas                                   2022   *
!     Modified by G. Gaigalas                               Jan 2024   *
!     Modified by G. Gaigalas                               Jan 2025   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use jj2lsj_profile, only: profile_finish, profile_start
      use jj2lsj_code, only: jj2lsj, profile_finish, profile_start
      use log_factorial_table, only: factt
      use isotope_data_input, only: setiso
      use run_timer_start, only: starttime
      use run_timer_finish, only: stoptime
      use physical_constants, only: setcon
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer           :: start_clock_count
!-----------------------------------------------
      call profile_start()
      call starttime (start_clock_count, "jj2lsj_omp")
      print *, " "
      print *, "jj2lsj_omp: Transformation of ASFs from a jj-coupled CSF basis"
      print *, "        into an LS-coupled CSF basis  (Fortran 95 version)"
      print *, "        (C) Copyright by   G. Gaigalas and Ch. F. Fischer,"
      print *, "        (2025)."
      print *, "        Input files: name_descriptors.parquet, name_header.toml, name.mix/name.cimix"
      print *, "         (optional)  name.lsj.T"
      print *, "        Ouput files: name.lsj.lbl,"
      print *, "         (optional)  name.lsj.c, name.lsj.j,"
      print *, "                     name.uni.lsj.lbl, name.uni.lsj.sum,"
      print *, "                     name.lsj.T"
      print *, " "
!
!  Set up the table of logarithms of factorials
      call setcon
      call factt
!
      call setiso("isodata")
      call jj2lsj
!
      call profile_finish()
      call stoptime (start_clock_count, "jj2lsj_omp")
      stop
      end program jj2lsj_omp
