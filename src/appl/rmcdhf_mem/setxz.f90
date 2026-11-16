!***********************************************************************
!                                                                      *
      subroutine setxz(j)
!                                                                      *
!   This subprogram sets the inhomogeneous terms to zero.              *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grid_c
      use int_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!
      xu(:radial_point_count) = 0.0_real64
      xv(:radial_point_count) = 0.0_real64
!
      return
      end subroutine setxz
