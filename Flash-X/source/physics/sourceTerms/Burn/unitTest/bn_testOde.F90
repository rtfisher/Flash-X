!> @copyright Copyright 2026 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief Synthetic ODE with controlled non-finite injection, for the Burn unit test

!> @ingroup physics_Burn
!!
!! @brief A decoupled linear ODE that can be made to return one NaN or Inf
!!
!! @details
!! @anchor bn_testOde
!!
!! The steppers in BurnIntegrate take @p derivs, @p jakob and @p bjakob as
!! procedure arguments, so they can be exercised on any right-hand side.  This
!! module supplies the simplest one that is still representative -- a decoupled
!! linear decay
!!
!!     dy_i/dt = -lambda_i y_i,   lambda_i = i
!!
!! whose Jacobian is diagonal and constant.  Using a synthetic RHS rather than a
!! real network is deliberate: the defect under test lives in the generic
!! step-acceptance logic shared by every network routed through BurnIntegrate,
!! so the test should not depend on reaction rates, a Helmholtz table, or any
!! particular nuclide set.
!!
!! On top of that the module can inject exactly one non-finite value.  Call
!! bn_testArm() with BN_TEST_NAN or BN_TEST_INF and the n-th subsequent
!! evaluation of the right-hand side writes that value into dydt(1) and then
!! disarms itself.  One shot is all the test needs: the poisoned trial step must
!! be rejected, and every retry after it is clean, so a correct integrator
!! converges normally on a reduced step.
!!
!! @note The injection is deterministic -- it depends on a call counter, not on
!!       overflow, roundoff or optimisation level.

module bn_testOde

  use, intrinsic :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan, &
                                           ieee_positive_inf

  implicit none

  !! what to inject
  integer, parameter :: BN_TEST_NONE = 0
  integer, parameter :: BN_TEST_NAN  = 1
  integer, parameter :: BN_TEST_INF  = 2

  !! decay constant of mode i is bn_test_lambda*i
  real, parameter :: bn_test_lambda = 1.0e0

  integer, save :: bn_test_kind    = BN_TEST_NONE  !! what to inject, if anything
  integer, save :: bn_test_at      = 0             !! inject on this call number
  integer, save :: bn_test_ncalls  = 0             !! evaluations since the last arm
  logical, save :: bn_test_fired   = .false.       !! has the one shot been spent?

contains

  !> Arm a single injection on evaluation number @p atCall, and reset the counter
  subroutine bn_testArm(kind, atCall)
    implicit none
    integer, intent(IN) :: kind, atCall

    bn_test_kind   = kind
    bn_test_at     = atCall
    bn_test_ncalls = 0
    bn_test_fired  = .false.

  end subroutine bn_testArm

  !> Disarm, so the right-hand side is finite everywhere
  subroutine bn_testDisarm()
    implicit none

    bn_test_kind   = BN_TEST_NONE
    bn_test_at     = 0
    bn_test_ncalls = 0
    bn_test_fired  = .false.

  end subroutine bn_testDisarm

  !> @brief Right-hand side, matching the derivs_t abstract interface
  !!
  !! @param tt      independent variable (unused -- the system is autonomous)
  !! @param y       current state
  !! @param btemp   temperature (unused; present for interface conformance)
  !! @param ratdum  reaction rates (unused; present for interface conformance)
  !! @param dydt    returned derivative
  subroutine bn_testDerivs(tt, y, btemp, ratdum, dydt)
    implicit none
    real, intent(IN)  :: tt, y(:), btemp, ratdum(:)
    real, intent(OUT) :: dydt(:)

    integer :: i, n

    n = size(y)
    do i = 1, n
       dydt(i) = -bn_test_lambda * real(i) * y(i)
    end do

    bn_test_ncalls = bn_test_ncalls + 1

    !! one shot, on the requested call
    if (.not. bn_test_fired .and. bn_test_kind /= BN_TEST_NONE &
         .and. bn_test_ncalls >= bn_test_at) then
       select case (bn_test_kind)
       case (BN_TEST_NAN)
          dydt(1) = ieee_value(1.0, ieee_quiet_nan)
       case (BN_TEST_INF)
          dydt(1) = ieee_value(1.0, ieee_positive_inf)
       end select
       bn_test_fired = .true.
    end if

  end subroutine bn_testDerivs

  !> @brief Analytic Jacobian, matching the jakob_t abstract interface
  !!
  !! Diagonal and constant.  Never poisoned: the point of the test is a bad
  !! trial step, not a bad Jacobian.
  subroutine bn_testJakob(tt, y, btemp, ratdum, dfdy, nzo, nDummy)
    implicit none
    integer, intent(IN) :: nzo, nDummy
    real, intent(IN)    :: tt, btemp, ratdum(:)
    real, intent(INOUT) :: y(:)
    real, intent(OUT)   :: dfdy(:,:)

    integer :: i

    dfdy = 0.0e0
    do i = 1, nzo
       dfdy(i,i) = -bn_test_lambda * real(i)
    end do

  end subroutine bn_testJakob

  !> @brief Sparse-pointer builder, matching the bjakob_t abstract interface
  !!
  !! The Gift steppers use a dense Jacobian, so this is only present to satisfy
  !! the interface.  It reports a diagonal.
  subroutine bn_testBjakob(iloc, jloc, nzo, np)
    implicit none
    integer, intent(IN)  :: iloc(:), jloc(:), np
    integer, intent(OUT) :: nzo

    nzo = np

  end subroutine bn_testBjakob

end module bn_testOde
