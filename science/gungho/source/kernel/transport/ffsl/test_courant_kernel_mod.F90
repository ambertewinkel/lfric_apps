!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Test kernel in the z direction. (Cell centred Courant number calculation)
!> @details Adapted from the Piecewise Constant Method applied to the z direction.
!!          Original file: ffsl_flux_z_constant_kernel_mod.F90

module test_courant_kernel_mod

use argument_mod,                   only : arg_type,              &
                                           GH_FIELD, GH_REAL,     &
                                           GH_READ, GH_WRITE,     &
                                           GH_SCALAR, CELL_COLUMN
use fs_continuity_mod,              only : W3, W2v
use constants_mod,                  only : r_tran, i_def, l_def, EPS_R_TRAN
use kernel_mod,                     only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: test_courant_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                  &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W3)   & ! Courant cell centre
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: test_courant_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: test_courant_code

contains

!> @brief Computes the mass flux for FFSL using PCM in the z direction.
!> @param[in]     nlayers   Number of layers
!> @param[in]     dep_dist  The vertical departure points
!> @param[in,out] C_w3      The 1D vertical Courant number at cell centres
!> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine test_courant_code( nlayers,   &
                                      dep_dist,  &
                                      C_w3,      &
                                      ndf_w2v,   &
                                      undf_w2v,  &
                                      map_w2v,   &
                                      ndf_w3,    &
                                      undf_w3,   &
                                      map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: ndf_w3
  real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
  real(kind=r_tran),   intent(inout) :: C_w3(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

  ! Internal variables
  integer(kind=i_def) :: k, w2v_idx, w3_idx

  real(kind=r_tran)   :: displacement, displacement_m1

  w2v_idx = map_w2v(1)
  w3_idx = map_w3(1)

  ! Calculate Courant number at cell centers ! not quite correct, I need to
  ! base it off the local cell centred volume, I am not sure how the departure
  ! point 'Courant number' is calculated.
  displacement_m1 = 0.0_r_tran
  do k = 0, nlayers - 1
    displacement = dep_dist(w2v_idx + k + 1) ! signed Courant number at W2v
    if (k == nlayers - 1) then 
      displacement = 0.0_r_tran ! update for top face
    end if
    C_w3(w3_idx + k) = 0.5_r_tran*(ABS(displacement_m1) + ABS(displacement))
    displacement_m1 = displacement ! update for lower face
  end do

end subroutine test_courant_code

end module test_courant_kernel_mod
