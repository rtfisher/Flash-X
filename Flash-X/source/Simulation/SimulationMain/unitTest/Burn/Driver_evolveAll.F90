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
!! @brief Driver_evolveAll for the Burn unit test
!!
!! @details
!! Opens the per-rank unitTest_NNNN file expected by the test harness, calls
!! Burn_unitTest, and writes the success line the harness looks for.  There is
!! no time evolution: the test drives the network integrators directly.

subroutine Driver_evolveAll()

#include "constants.h"

  use Driver_data,    ONLY: dr_globalMe
  use Burn_interface, ONLY: Burn_unitTest

  implicit none

  interface
     integer function ut_getFreeFileUnit()
     end function ut_getFreeFileUnit
  end interface

  logical :: perfect
  character(len=20) :: fileName
  integer           :: fileUnit
  integer, dimension(4) :: prNum
  integer :: temp, i

  ! stays true if no errors are found
  perfect = .true.

  temp = dr_globalMe
  do i = 1, 4
     prNum(i) = mod(temp, 10)
     temp = temp / 10
  end do

  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                          char(48+prNum(2))//char(48+prNum(1))

  fileUnit = ut_getFreeFileUnit()
  open(fileUnit, file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe

  print *, "Preparing to call Burn_unitTest"
  call Burn_unitTest(fileUnit, perfect)
  print *, "    Returned from Burn_unitTest"

  if (perfect) then
     write(fileUnit,'("all results conformed with expected values.")')
  end if

  close(fileUnit)

  return

end subroutine Driver_evolveAll
