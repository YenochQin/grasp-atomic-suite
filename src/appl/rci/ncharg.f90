!***********************************************************************
!                                                                      *
      subroutine ncharg
!                                                                      *
!   This routine evaluates the nuclear charge density, and stores it   *
!   in the  common  array  ZDIST .                                     *
!                                                                      *
!   Call(s) to: [LIB92]: ES.                                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use def_c,           only: pi,  nuclear_charge, precis
      use grid_c
      use npar_c
      use ncdist_c
      use tatb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use es_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: c, a, cba, pi2, abc, abc2, abc3, s2mcba, s3mcba, en, &
         znorm, extrm, zdisti
      logical :: form1, form2
!-----------------------------------------------
!
!
!   Initialize array to zero
!
      zdist(:radial_point_count) = 0.0d00
!
!   Fermi charge distribution
!
      if (nuclear_model_parameter_count == 2) then
         c = nuclear_model_parameters(1)
         a = nuclear_model_parameters(2)
         cba = c/a
         pi2 = pi*pi
         abc = a/c
         abc2 = abc*abc
         abc3 = abc2*abc
         call es ((-cba), s2mcba, s3mcba)
         en = 1.0d00 + pi2*abc2 - 6.0d00*abc3*s3mcba
         znorm = 3.0d00*nuclear_charge/(4.0d00*pi*en*c**3)
         form1 = .true.
         form2 = .false.
         do i = 1, radial_point_count
            if (form1) then
               extrm = exp((radial_grid(i)-c)/a)
               zdist(i) = znorm/(1.0d00 + extrm)
               if (1.0d00/extrm <= precis) then
                  form1 = .false.
                  form2 = .true.
               endif
            else if (form2) then
               zdisti = znorm*exp((-(radial_grid(i)-c)/a))
               if (abs(zdisti) > 0.0d00) then
                  zdist(i) = zdisti
               else
                  mtp = i
                  exit
               endif
            endif
         end do
      endif
!
      return
!
      end subroutine ncharg
