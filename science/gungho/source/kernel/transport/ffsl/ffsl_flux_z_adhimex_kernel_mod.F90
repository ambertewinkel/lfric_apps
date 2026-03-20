!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   AdHImEx kernel in the z direction.
!> @details Adapted from the Piecewise Constant Method applied to the z direction.
!!          Original file: ffsl_flux_z_constant_kernel_mod.F90

module ffsl_flux_z_adhimex_kernel_mod

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
type, public, extends(kernel_type) :: ffsl_flux_z_adhimex_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                  &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! flux
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! detj
       arg_type(GH_SCALAR, GH_REAL,    GH_READ)        & ! dt
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: ffsl_flux_z_adhimex_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: ffsl_flux_z_adhimex_code
public :: fifth_order_adhimex


contains

!> @brief Computes the mass flux for FFSL using AdHImEx in the z direction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] flux      The flux to be computed
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
!> @param[in]     field     The field to construct the flux
!> @param[in]     detj      Volume of cells
!> @param[in]     dt        Time step
!> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine ffsl_flux_z_adhimex_code( nlayers,    &
                                      flux,      &
                                      dep_dist,  &
                                      field,     &
                                      detj,      &
                                      dt,        &
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
  real(kind=r_tran),   intent(inout) :: flux(undf_w2v)
  real(kind=r_tran),   intent(in)    :: field(undf_w3)
  real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
  real(kind=r_tran),   intent(in)    :: detj(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  real(kind=r_tran),   intent(in)    :: dt

  ! Internal variables
  integer(kind=i_def) :: k, s, i_s, w2v_idx, w3_idx

  integer(kind=i_def), parameter :: nstages = 5

  real(kind=r_tran)   :: courant(nlayers)             ! Courant number at cell centres
  real(kind=r_tran)   :: implness_w2v(nlayers + 1)    ! implicitness at faces (div 1D)
  real(kind=r_tran)   :: implness_w3(nlayers)         ! implicitness at centres
  real(kind=r_tran)   :: a_ex(nstages, nstages)       ! explicit Butcher tableau
  real(kind=r_tran)   :: a_im(nstages, nstages)       ! implicit Butcher tableau
  real(kind=r_tran)   :: field_s(nlayers)             ! stage field
  real(kind=r_tran)   :: field_s_w2v(nlayers + 1)     ! interpolated field_s at faces
  real(kind=r_tran)   :: c_field_s_w2v(nlayers + 1)   ! interpolated field_s*C at faces
  real(kind=r_tran)   :: f_ex_adv(nstages, nlayers)   ! ex advective divergence
  real(kind=r_tran)   :: f_im_adv(nstages, nlayers)   ! im advective divergence
  real(kind=r_tran)   :: f_ex_con(nstages, nlayers)   ! ex conservative divergence
  real(kind=r_tran)   :: f_im_con(nstages, nlayers)   ! im conservative divergence
  real(kind=r_tran)   :: zero

  w2v_idx = map_w2v(1)
  w3_idx = map_w3(1)
  zero = 0.0_r_tran

  ! Calculate Courant number and implicitness - assumes uniform vertical grid
  courant(1) = 0.5_r_tran*(ABS(dep_dist(w2v_idx)) &
                               + ABS(dep_dist(w2v_idx + 1)))
  implness_w2v(1) = zero
  do k = 1, nlayers - 1
    courant(k + 1) = 0.5_r_tran*(ABS(dep_dist(w2v_idx + k)) + ABS(dep_dist(w2v_idx + k + 1)))
    implness_w2v(k + 1) = 1.0_r_tran - 1.0_r_tran/(1.0_r_tran + 0.7_r_tran*(MAX( &
                          1.4_r_tran, courant(k), courant(k + 1)) - 1.4_r_tran))
  end do
  implness_w2v(nlayers + 1) = zero
  implness_w3 = MAX( implness_w2v(1 : nlayers), implness_w2v(2 : nlayers + 1) )

  ! Set up Butcher tableau (remember column-major order of reshape)
  a_ex = reshape((/ zero, zero, zero, zero, zero,                               &
                    zero, zero, 1.0_r_tran, 0.25_r_tran, 1.0_r_tran/6.0_r_tran, &
                    zero, zero, zero, 0.25_r_tran, 1.0_r_tran/6.0_r_tran,       &
                    zero, zero, zero, zero, 2.0_r_tran/3.0_r_tran,              &
                    zero, zero, zero, zero, zero /), shape(a_ex))
  a_im = reshape((/ zero, 0.5_r_tran, 0.5_r_tran, 0.5_r_tran, 0.5_r_tran,       &
                    zero, zero, zero, zero, zero,                               &
                    zero, zero, zero, zero, zero,                               &
                    zero, zero, zero, zero, zero,                               &
                    zero, zero, zero, zero, 0.5_r_tran /), shape(a_im))

  ! Setting all f_... elements to zero
  f_ex_adv = zero
  f_im_adv = zero
  f_ex_con = zero
  f_im_con = zero
  flux(w2v_idx : w2v_idx + nlayers) = zero

  do s = 1, nstages
    ! Calculate rhs (depends on stage for constancy)
    field_s = field(w3_idx : w3_idx + nlayers - 1)
    if ( (s == 2) .or. (s == 3) ) then ! advective
      do i_s = 1, s
        field_s = field_s + f_ex_adv(i_s,:)*a_ex(s, i_s) + f_im_adv(i_s,:)*a_im(s, i_s)
      end do
    else ! conservative
      do i_s = 1, s
        field_s = field_s + f_ex_con(i_s,:)*a_ex(s, i_s) + f_im_con(i_s,:)*a_im(s, i_s)
      end do
    end if

    ! Call matrix solver for last stage if required
    if ( (s == 5) .and. (any(implness_w2v /= zero)) ) then
      field_s = field_s ! todo
    end if

    ! Interpolate field_s to faces
    call fifth_order_adhimex( nlayers, &
                    field_s_w2v, field_s, dep_dist(w2v_idx : w2v_idx + nlayers) )
    c_field_s_w2v = dep_dist(w2v_idx : w2v_idx + nlayers)*field_s_w2v

    ! Calculate various divergences (f_...) based on the new stage field
    do k = 1, nlayers
      f_ex_adv(s, k) = - ( (1.0_r_tran - implness_w2v(k + 1))*c_field_s_w2v(k + 1) &
                               - (1.0_r_tran - implness_w2v(k))*c_field_s_w2v(k) )
      f_im_adv(s, k) = - ( implness_w2v(k + 1)*c_field_s_w2v(k + 1)                &
                               - implness_w2v(k)*c_field_s_w2v(k) )
      f_ex_con(s, k) = - (1.0_r_tran - implness_w3(k))*( c_field_s_w2v(k + 1)      &
                               - c_field_s_w2v(k) )
      f_im_con(s, k) = - implness_w3(k)*( c_field_s_w2v(k + 1)                     &
                               - c_field_s_w2v(k) )
    end do

    ! Update total flux
    flux(w2v_idx) = 0.0_r_tran
    do k = 2, nlayers
      flux(w2v_idx + k - 1) = flux(w2v_idx + k - 1) + (a_ex(nstages,s)&
                                    *(1.0_r_tran - implness_w2v(k)) + &
          a_im(nstages,s)*implness_w2v(k))*c_field_s_w2v(k)*detj(w3_idx + k - 1)/dt
    end do
    flux(w2v_idx + nlayers) = 0.0_r_tran
  end do

end subroutine ffsl_flux_z_adhimex_code

!> @brief Calculates fifth-order interpolated field at faces for AdHImEx
!> @param[in]     nl   Number of layers
!> @param[in,out] fieldh    The flux to be computed
!> @param[in]     field     The field to construct the flux
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
subroutine fifth_order_adhimex( nl,              &
                                      fieldh,    &
                                      field,     &
                                      dep_dist )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                ! nlayers
  real(kind=r_tran),   intent(in)    :: field(nl)         ! field to interpolate
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)  ! Courant (used for sign)
  real(kind=r_tran),   intent(inout) :: fieldh(nl + 1)    ! interpolated field

  ! Internal variables
  integer(kind=i_def) :: k

  ! Loop over most of the spatial domain (excluding boundaries)
  do k = 3, nl - 3
    fieldh(k + 1) = (MAX(0.0_r_tran, dep_dist(k + 1))                  & ! u>=0
      *( 2.0_r_tran*field(k - 2) - 13.0_r_tran*field(k - 1)            &
      + 47.0_r_tran*field(k)                                           &
      + 27.0_r_tran*field(k + 1) - 3.0_r_tran*field(k + 2) )           &
      - MIN(0.0_r_tran, dep_dist(k + 1))                               & ! u<0
      *( -3.0_r_tran*field(k - 1) + 27.0_r_tran*field(k)               &
      + 47.0_r_tran*field(k + 1) - 13.0_r_tran*field(k + 2)            &
      + 2.0_r_tran*field(k + 3) ))/(60.0_r_tran*ABS(dep_dist(k + 1)))
  end do

  ! Lower boundary
  fieldh(1) = (21.0_r_tran*field(1) - field(2))/20.0_r_tran              ! u>=0
  fieldh(2) = (MAX(0.0_r_tran, dep_dist(2))                            & ! u>=0
        *( 36.0_r_tran*field(1) + 27.0_r_tran*field(2)                 &
        - 3.0_r_tran*field(3) )                                        &
        - MIN(0.0_r_tran, dep_dist(2))                                 & ! u<0
        *( 24.0_r_tran*field(1) + 47.0_r_tran*field(2)                 &
        - 13.0_r_tran*field(3) + 2.0_r_tran*field(4) ))                &
        /(60.0_r_tran*ABS(dep_dist(2)))
  fieldh(3) = (MAX(0.0_r_tran, dep_dist(3))                            & ! u>=0
        *( -11.0_r_tran*field(1) + 47.0_r_tran*field(2)                &
        + 27.0_r_tran*field(3) - 3.0_r_tran*field(4) )                 &
        - MIN(0.0_r_tran, dep_dist(3))                                 & ! u<0
        *( -3.0_r_tran*field(1) + 27.0_r_tran*field(2)                 &
        + 47.0_r_tran*field(3) - 13.0_r_tran*field(4)                  &
        + 2.0_r_tran*field(5) ))/(60.0_r_tran*ABS(dep_dist(3)))

  ! Upper boundary
  fieldh(nl + 1) = (-field(nl - 1) + 21.0_r_tran*field(nl))/20.0_r_tran  ! u<=0
  fieldh(nl) = (MAX(0.0_r_tran, dep_dist(nl))                          & ! u>=0
                                 *( 2.0_r_tran*field(nl - 3)           &
                                 - 13.0_r_tran*field(nl - 2)           &
                                 + 47.0_r_tran*field(nl - 1)           &
                                 + 24.0_r_tran*field(nl) )             &
                                 - MIN(0.0_r_tran, dep_dist(nl))       & ! u<0
                                 *( -3.0_r_tran*field(nl - 2)          &
                                 + 27.0_r_tran*field(nl - 1)           &
                                 + 36.0_r_tran*field(nl) ))            &
                                 /(60.0_r_tran*ABS(dep_dist(nl)))
  fieldh(nl - 1) = (MAX(0.0_r_tran, dep_dist(nl - 1))                  & ! u>=0
                            *( 2.0_r_tran*field(nl - 4)                &
                            - 13.0_r_tran*field(nl - 3)                &
                            + 47.0_r_tran*field(nl - 2)                &
                            + 27.0_r_tran*field(nl - 1)                &
                            - 3.0_r_tran*field(nl) )                   &
                            - MIN(0.0_r_tran, dep_dist(nl - 1))        & ! u<0
                            *( -3.0_r_tran*field(nl - 3)               &
                            + 27.0_r_tran*field(nl - 2)                &
                            + 47.0_r_tran*field(nl - 1)                &
                            - 11.0_r_tran*field(nl)))                  &
                            /(60.0_r_tran*ABS(dep_dist(nl - 1)))

end subroutine fifth_order_adhimex

end module ffsl_flux_z_adhimex_kernel_mod
