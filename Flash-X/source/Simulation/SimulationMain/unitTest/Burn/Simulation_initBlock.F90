!! @file
!! @brief Simulation_initBlock for the Burn unit test
!!
!! The test drives the network integrators directly and never reads the mesh,
!! so this fills the block uniformly and does no more.

subroutine Simulation_initBlock(solnData, tileDesc)

  use Simulation_data, ONLY: sim_initDens, sim_initTemp
  use Grid_tile,       ONLY: Grid_tile_t

  implicit none

# include "constants.h"
# include "Simulation.h"

  real, dimension(:,:,:,:), pointer :: solnData
  type(Grid_tile_t), intent(in)     :: tileDesc

  integer :: i, j, k, n
  integer, dimension(LOW:HIGH, MDIM) :: lim

  lim = tileDesc%limits

  do k = lim(LOW,KAXIS), lim(HIGH,KAXIS)
     do j = lim(LOW,JAXIS), lim(HIGH,JAXIS)
        do i = lim(LOW,IAXIS), lim(HIGH,IAXIS)
           solnData(DENS_VAR,i,j,k) = sim_initDens
           solnData(TEMP_VAR,i,j,k) = sim_initTemp
           do n = SPECIES_BEGIN, SPECIES_END
              solnData(n,i,j,k) = 1.0e0 / real(NSPECIES)
           end do
        end do
     end do
  end do

end subroutine Simulation_initBlock
