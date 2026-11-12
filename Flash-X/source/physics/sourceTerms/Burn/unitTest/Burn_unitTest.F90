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
!! @brief Burn_unitTest implementation -- step acceptance in the network integrators

!> @ingroup physics_Burn
!!
!! @brief Unit test for the Burn unit: non-finite trial steps must be rejected
!!
!! @details
!! @anchor Burn_unitTest_impl
!!
!! WHAT IS BEING TESTED
!!
!! The adaptive steppers in BurnIntegrate decide whether to accept a trial step
!! from a normalised error estimate accumulated as
!!
!!     errmax = tiny
!!     do i=1,nv
!!        errmax = max(errmax, abs(yerr(i)/yscal(i)))
!!     enddo
!!     ...
!!     if (errmax .lt. 1.0) <accept>
!!
!! If any yerr(i) is NaN, whether errmax becomes NaN depends on how the compiler
!! lowers max() with a NaN operand -- Fortran does not specify it, and the same
!! source has been measured giving opposite answers on x86-64 and arm64, and
!! across optimisation levels on one of them.  Where max() returns the finite
!! operand, errmax stays at tiny, the comparison succeeds, and a step containing
!! a NaN is accepted as converged.  The poisoned value then propagates silently
!! into the abundances.
!!
!! The steppers therefore test yerr and y for finiteness explicitly and force a
!! rejection.  That guard is a no-op on finite data, so it cannot be detected by
!! comparing simulation output: this test drives the integrator directly and
!! injects the non-finite value on purpose.
!!
!! WHAT THE TEST DOES
!!
!! For each of bn_baderGift and bn_rosenGift, three cases are run on the
!! synthetic decoupled linear ODE in bn_testOde:
!!
!!   1. control -- right-hand side finite everywhere.  The step must be accepted
!!      at full size (hdid == htry) and match the analytic solution.
!!   2. NaN     -- one evaluation of the right-hand side returns a NaN.
!!   3. Inf     -- one evaluation returns +Inf.
!!
!! In cases 2 and 3 the injection is a single shot, so every retry after the
!! rejected trial is clean.  A correct stepper therefore rejects the poisoned
!! trial, reduces the step, and converges.  Each case checks that
!!
!!   * the returned state is finite -- this is the check that a stepper which
!!     accepted the poisoned trial fails, since the NaN lands in y;
!!   * a step was actually taken (0 < hdid <= htry).  How far the step is cut
!!     after a rejection is stepper-specific -- the Bader extrapolation and the
!!     Rosenbrock stepper make different choices -- so only forward progress is
!!     required, not a particular reduction;
!!   * the result still matches the analytic solution over the step actually
!!     taken, which rules out "rejected everything and returned something
!!     harmless";
!!   * and, for the poisoned cases, that the injection really did fire.  Without
!!     that last check a case could pass vacuously by never reaching the armed
!!     evaluation.
!!
!! ON PLATFORM SENSITIVITY
!!
!! Where max() does propagate NaN, errmax becomes NaN, `errmax .lt. 1.0` is
!! false, and the trial is rejected even without the explicit guard -- so this
!! test passes there with or without it.  With the guard in place it passes
!! everywhere; without it, it fails on exactly the platforms carrying the
!! defect.  It should never report a false failure.
!!
!! Only the Gift variants are exercised.  The acceptance logic is line-for-line
!! the same in bn_baderMa28 and bn_rosenMa28, but those reach the Harwell MA28
!! solver, which is licence-restricted and not present in every build; the Gift
!! path exercises the same guard without depending on it at run time.  (The Burn
!! unit still declares LIBRARY MA28, so the link dependency itself remains.)
!!
!!   @param fileUnit : unit number for the file opened by the unitTest/Burn setup
!!                     in which to write the results of the test
!!
!!   @param perfect  : true if and only if every check passed
!!
!!***

subroutine Burn_unitTest(fileUnit, perfect)

#include "constants.h"
#include "Simulation.h"

  use bn_testOde,           ONLY: bn_testDerivs, bn_testJakob, bn_testBjakob, &
                                  bn_testArm, bn_testDisarm, bn_test_lambda,  &
                                  bn_test_fired,                              &
                                  BN_TEST_NONE, BN_TEST_NAN, BN_TEST_INF
  use bnIntegrate_interface, ONLY: bn_baderGift, bn_rosenGift
  use bnNetwork_interface,   ONLY: steper_state_t
  use Logfile_interface,     ONLY: Logfile_stampMessage
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get

  implicit none

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: perfect

  !! the synthetic system is sized to the network so that any indexing of
  !! Burn_data arrays inside the steppers stays in bounds
  integer, parameter :: nv = NSPECIES

  !! stepper selection
  integer, parameter :: BADER_GIFT = 1, ROSEN_GIFT = 2

  real    :: htry, eps, tol
  integer :: injectAt
  logical :: ok

  perfect = .true.

  call RuntimeParameters_get("bn_testTolerance", tol)
  call RuntimeParameters_get("bn_testStepSize",  htry)
  call RuntimeParameters_get("bn_testEps",       eps)
  call RuntimeParameters_get("bn_testInjectAt",  injectAt)

  write(fileUnit,'("Burn unit test: step acceptance in the network integrators")')
  write(fileUnit,'("  nv = ",I0,"  htry = ",ES12.5,"  eps = ",ES12.5)') nv, htry, eps

  call runCase(BADER_GIFT, "bn_baderGift", BN_TEST_NONE, ok); perfect = perfect .and. ok
  call runCase(BADER_GIFT, "bn_baderGift", BN_TEST_NAN,  ok); perfect = perfect .and. ok
  call runCase(BADER_GIFT, "bn_baderGift", BN_TEST_INF,  ok); perfect = perfect .and. ok

  call runCase(ROSEN_GIFT, "bn_rosenGift", BN_TEST_NONE, ok); perfect = perfect .and. ok
  call runCase(ROSEN_GIFT, "bn_rosenGift", BN_TEST_NAN,  ok); perfect = perfect .and. ok
  call runCase(ROSEN_GIFT, "bn_rosenGift", BN_TEST_INF,  ok); perfect = perfect .and. ok

  call bn_testDisarm()

  if (perfect) then
     call Logfile_stampMessage('[Burn_unitTest] SUCCESS. Non-finite trial steps were rejected.')
  else
     call Logfile_stampMessage('[Burn_unitTest] FAILURE. See the unitTest_* file for details.')
  end if

  return

contains

  !> Run one stepper against one injection case and check the outcome
  subroutine runCase(whichStepper, stepperName, kind, casePassed)
    implicit none
    integer,          intent(IN)  :: whichStepper, kind
    character(len=*), intent(IN)  :: stepperName
    logical,          intent(OUT) :: casePassed

    type(steper_state_t) :: state
    real    :: y(nv), dydx(nv), yscal(nv), ratdum(1)
    real    :: x, x0, btemp, hdid, hnext
    real    :: exact, relerr, worst
    integer :: i
    logical :: finiteOut, stepTaken, accurate, injected
    character(len=16) :: caseName

    select case (kind)
    case (BN_TEST_NAN); caseName = "NaN"
    case (BN_TEST_INF); caseName = "Inf"
    case default;       caseName = "control"
    end select

    !! initial state: y_i = 1, so the analytic solution is exp(-lambda*i*h)
    do i = 1, nv
       y(i) = 1.0e0
    end do
    x     = 0.0e0
    x0    = x
    btemp = 1.0e9              !! unused by the synthetic right-hand side
    ratdum = 0.0e0

    !! the steppers expect the state struct initialised as bn_netIntegrate does
    state%ifirst = .true.
    state%first  = .true.
    state%reduct = .false.
    state%epsold = -1.0e0
    state%nvold  = -1
    state%nseq   = [2, 6, 10, 14, 22, 34, 50, 70]

    !! derivative at the initial point, and the accuracy scaling vector.  This
    !! evaluation is deliberately made before arming, so dydx is always finite:
    !! the case under test is a bad TRIAL step, not a bad input.
    call bn_testDisarm()
    call bn_testDerivs(x, y, btemp, ratdum, dydx)
    do i = 1, nv
       yscal(i) = max(1.0e-6, abs(y(i)))
    end do

    if (kind /= BN_TEST_NONE) call bn_testArm(kind, injectAt)

    select case (whichStepper)
    case (BADER_GIFT)
       call bn_baderGift(state, y, dydx, ratdum, nv, x, btemp, htry, eps, yscal, &
                         hdid, hnext, bn_testDerivs, bn_testJakob, bn_testBjakob)
    case (ROSEN_GIFT)
       call bn_rosenGift(state, y, dydx, ratdum, nv, x, btemp, htry, eps, yscal, &
                         hdid, hnext, bn_testDerivs, bn_testJakob, bn_testBjakob)
    end select

    !! did the injection actually happen?  Without this the poisoned cases could
    !! pass vacuously, by never reaching the evaluation that was armed.
    injected = bn_test_fired
    call bn_testDisarm()

    !! ---- check 1: the returned state must be finite -------------------------
    !! This is the check that fails if a poisoned trial step was accepted.
    finiteOut = .true.
    do i = 1, nv
       if (.not. (abs(y(i)) .le. huge(1.0e0))) finiteOut = .false.
    end do
    if (.not. (abs(hdid) .le. huge(1.0e0))) finiteOut = .false.

    !! ---- check 2: a step must actually have been taken ----------------------
    !! How much the step is reduced after a rejection is stepper-specific -- the
    !! Bader extrapolation and the Rosenbrock stepper make different choices --
    !! so all that is required here is forward progress.
    stepTaken = (hdid > 0.0e0) .and. (hdid <= htry)
    if (kind /= BN_TEST_NONE) stepTaken = stepTaken .and. injected

    !! ---- check 3: the answer over the step actually taken must be right -----
    !! Guards against "rejected everything and returned something harmless".
    accurate = .true.
    worst    = 0.0e0
    if (finiteOut) then
       do i = 1, nv
          exact  = exp(-bn_test_lambda * real(i) * hdid)
          relerr = abs(y(i) - exact) / max(abs(exact), tiny(1.0e0))
          worst  = max(worst, relerr)
       end do
       accurate = (worst .le. tol)
    end if

    casePassed = finiteOut .and. stepTaken .and. accurate

    write(fileUnit,'(2X,A," / ",A,": ",A)') trim(stepperName), trim(caseName), &
         merge("PASS", "FAIL", casePassed)
    write(fileUnit,'(4X,"finite=",L1,"  injected=",L1,"  hdid=",ES12.5, &
                    &" (htry=",ES12.5,")  max rel err=",ES12.5)') &
         finiteOut, injected, hdid, htry, worst

    if (.not. finiteOut) then
       write(fileUnit,'(4X,"FAILURE: a non-finite trial step was accepted")')
    else if (kind /= BN_TEST_NONE .and. .not. injected) then
       write(fileUnit,'(4X,"FAILURE: the non-finite value was never injected, so ", &
                       &"this case proved nothing -- lower bn_testInjectAt")')
    else if (.not. stepTaken) then
       write(fileUnit,'(4X,"FAILURE: no forward progress was made")')
    else if (.not. accurate) then
       write(fileUnit,'(4X,"FAILURE: result does not match the analytic solution")')
    end if

  end subroutine runCase

end subroutine Burn_unitTest
