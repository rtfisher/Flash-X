!! @file
!! @brief Simulation_init for the Burn unit test

subroutine Simulation_init()

  use Simulation_data,             ONLY: sim_meshMe, sim_initDens, sim_initTemp
  use Driver_interface,            ONLY: Driver_getMype
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get

  implicit none

# include "constants.h"

  call Driver_getMype(MESH_COMM, sim_meshMe)
  call RuntimeParameters_get('sim_initDens', sim_initDens)
  call RuntimeParameters_get('sim_initTemp', sim_initTemp)

end subroutine Simulation_init
