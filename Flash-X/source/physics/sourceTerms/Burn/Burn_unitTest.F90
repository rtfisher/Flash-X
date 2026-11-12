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
!! @brief Burn_unitTest stub

!> @ingroup physics_Burn
!!
!! @brief unit test driver for the Burn unit
!!
!! @details
!! @anchor Burn_unitTest_stub
!!
!! This function is the unit test for the Burn unit.  It is invoked in the
!! setup unitTest/Burn, which REQUESTS physics/sourceTerms/Burn/unitTest to
!! replace this stub with the real implementation.  See that implementation for
!! a description of what is checked.
!!
!!   @param fileUnit : unit number for the file opened by the unitTest/Burn setup
!!                     in which to write the results of the test
!!
!!   @param perfect  : indicates the test ran without error if true
!!
!!***

subroutine Burn_unitTest(fileUnit, perfect)

  use Logfile_interface, ONLY: Logfile_stampMessage

  implicit none

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: perfect

  ! Stub version returns false
  perfect = .false.

  write(*,*) 'Burn_unitTest FAILURE:  This version is a stub'
  call Logfile_stampMessage('[Burn_unitTest] FAILURE.  This version is a stub.')

  return

end subroutine Burn_unitTest
