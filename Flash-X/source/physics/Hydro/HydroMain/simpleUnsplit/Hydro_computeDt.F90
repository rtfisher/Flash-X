!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_Unsplit/Hydro_computeDt
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  Hydro_computeDt
!!
!!
!! 
!! DESCRIPTION
!!
!!  This routine computes the timestep limiter for the simplified  Hydro solver.
!!  For more details see the documentation of the NULL implementation
!!
!!***

!!REORDER(4): U
#include "Simulation.h"
#undef FIXEDBLOCKSIZE
Subroutine Hydro_computeDt(tileDesc, &
     x, dx, uxgrid, &
     y, dy, uygrid, &
     z, dz, uzgrid, &
     blkLimits, blkLimitsGC, &
     U,  dtCheck, dtMinLoc,  &
     extraInfo)


#include "constants.h"

  use Hydro_data,       ONLY : hy_geometry, hy_cfl, hy_dref, hy_eref, &
                               hy_pref, hy_vref, hy_meshMe,           &
                               hy_useHydro, hy_updateHydroFluxes,     &
                               hy_useVaryingCFL
  use Driver_interface, ONLY : Driver_abort
  use Grid_tile,        ONLY : Grid_tile_t
  implicit none

  !! Arguments type declaration ------------------------------------------
  type(Grid_tile_t), intent(IN) :: tileDesc
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC

#ifdef FIXEDBLOCKSIZE
  ! DEV: FIXME This should not be here.  blkLimitsGC should have the correct
  !            bounds.  Also, it's not in the Unsplit version of this routine.
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif

  real, pointer         :: U(:,:,:,:)
  real,   intent(INOUT) :: dtCheck
  integer,intent(INOUT) :: dtMinLoc(5)
  real, OPTIONAL,intent(INOUT) :: extraInfo
  !! ----------------------------------------------------------------------

  integer :: i, j, k, temploc(5)
  real    :: sndspd2, delxinv, delyinv, delzinv, dt_temp, dt_ltemp


  if ((.not. hy_useHydro) .or. (.not. hy_updateHydroFluxes)) return

  dt_temp    = 0.
  temploc(:) = 0

  !! Conversion for unitSystem if needed ---------------------------------------
  associate(lo => blkLimits(LOW,  :), &
            hi => blkLimits(HIGH, :))
     do       k = lo(KAXIS), hi(KAXIS)
        do    j = lo(JAXIS), hi(JAXIS)
           do i = lo(IAXIS), hi(IAXIS)
              U(DENS_VAR,i,j,k) = U(DENS_VAR,i,j,k)/hy_dref
              U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k)/hy_eref
              U(PRES_VAR,i,j,k) = U(PRES_VAR,i,j,k)/hy_pref
              U(VELX_VAR:VELZ_VAR,i,j,k) = U(VELX_VAR:VELZ_VAR,i,j,k)/hy_vref
           end do
        end do
     end do
  end associate
  !! ---------------------------------------------------------------------------


  if (NDIM == 1) then
     !--------------------------------------------------------------------------!
     ! 1-dimensional                                                            !
     !--------------------------------------------------------------------------!
     delxinv = 1.0/dx(blkLimits(LOW,IAXIS))

     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

        sndspd2   = U(GAMC_VAR,i,1,1)*U(PRES_VAR,i,1,1)/U(DENS_VAR,i,1,1) 
        dt_ltemp  = (abs(U(VELX_VAR,i,1,1)-uxgrid(i))+sqrt(sndspd2))*delxinv

        if (dt_ltemp > dt_temp) then
           dt_temp    = dt_ltemp
           temploc(1) = i
           temploc(2) = 1
           temploc(3) = 1
           temploc(4) = tileDesc%level
           temploc(5) = hy_meshMe
        endif

     enddo

  elseif (NDIM == 2) then
     !--------------------------------------------------------------------------!
     ! 2-dimensional                                                            !
     !--------------------------------------------------------------------------!
     if (hy_geometry == CARTESIAN .OR. hy_geometry == CYLINDRICAL) then 
        ! the 'y' coordinate is not angular in Cartesian and cylindrical coords

        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
        delyinv = 1.0/dy(blkLimits(LOW,JAXIS))

        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              sndspd2  = U(GAMC_VAR,i,j,1)*U(PRES_VAR,i,j,1)/U(DENS_VAR,i,j,1)

              dt_ltemp = max((abs(U(VELX_VAR,i,j,1)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                   (abs(U(VELY_VAR,i,j,1)-uygrid(j))+sqrt(sndspd2))*delyinv)

              if (dt_ltemp > dt_temp) then
                 dt_temp    = dt_ltemp
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = 1
                 temploc(4) = tileDesc%level
                 temploc(5) = hy_meshMe
              endif

           enddo
        enddo
               
     else ! Angular coordinates in 2D: Spherical or Polar

        ! y is angular
        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))

        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              delyinv = 1.0/(x(i)*dy(j))

              sndspd2  = U(GAMC_VAR,i,j,1)*U(PRES_VAR,i,j,1)/U(DENS_VAR,i,j,1)
              dt_ltemp = max((abs(U(VELX_VAR,i,j,1)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                   (abs(U(VELY_VAR,i,j,1)-uygrid(j))+sqrt(sndspd2))*delyinv)

              if (dt_ltemp > dt_temp) then
                 dt_temp    = dt_ltemp
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = 1
                 temploc(4) = tileDesc%level
                 temploc(5) = hy_meshMe
              endif
           enddo
        enddo

     endif

  elseif (NDIM == 3) then
     !--------------------------------------------------------------------------!
     ! 3-dimensional                                                            !
     !--------------------------------------------------------------------------!

     if (hy_geometry == CARTESIAN) then

        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
        delyinv = 1.0/dy(blkLimits(LOW,JAXIS))
        delzinv = 1.0/dz(blkLimits(LOW,KAXIS))

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 sndspd2  = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 dt_ltemp = max((abs(U(VELX_VAR,i,j,k)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                      (abs(U(VELY_VAR,i,j,k)-uygrid(j))+sqrt(sndspd2))*delyinv,&
                      (abs(U(VELZ_VAR,i,j,k)-uzgrid(k))+sqrt(sndspd2))*delzinv)

                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = tileDesc%level
                    temploc(5) = hy_meshMe
                 endif

              enddo
           enddo
        enddo

     elseif (hy_geometry == CYLINDRICAL) then

        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
        delyinv = 1.0/dy(blkLimits(LOW,JAXIS))

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 delzinv = 1.0/(x(i)*dz(k)) ! z is phi

                 sndspd2  = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 dt_ltemp = max((abs(U(VELX_VAR,i,j,k)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                      (abs(U(VELY_VAR,i,j,k)-uygrid(j))+sqrt(sndspd2))*delyinv,&
                      (abs(U(VELZ_VAR,i,j,k)-uzgrid(k))+sqrt(sndspd2))*delzinv)

                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = tileDesc%level
                    temploc(5) = hy_meshMe
                 endif

              enddo
           enddo
        enddo
        
     elseif (hy_geometry == SPHERICAL) then

        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 delyinv = 1.0/(x(i)*dy(j)) ! y is theta
                 delzinv = 1.0/(x(i)*sin(y(j))*dz(k)) ! z is phi

                 sndspd2  = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 dt_ltemp = max((abs(U(VELX_VAR,i,j,k)-uxgrid(i))+sqrt(sndspd2))*delxinv,&
                      (abs(U(VELY_VAR,i,j,k)-uygrid(j))+sqrt(sndspd2))*delyinv,&
                      (abs(U(VELZ_VAR,i,j,k)-uzgrid(k))+sqrt(sndspd2))*delzinv)

                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = tileDesc%level
                    temploc(5) = hy_meshMe
                 endif

              enddo
           enddo
        enddo

     else ! Polar in 3D (that's a no no)
        call Driver_abort("[Hydro_computeDt] ERROR: Polar geometry not supported in 3D")
     endif

  endif


  dt_temp = hy_cfl / dt_temp
  if (dt_temp < dtCheck) then
     dtCheck = dt_temp
     dtMinLoc = temploc
  endif

  !! For the purpose of having screen output for CFL
  if (present(extraInfo)) then
     if (hy_useVaryingCFL) then
        extraInfo = hy_cfl
     else
        extraInfo = 0.
     end if
  end if

  !! Conversion for unitSystem if needed ---------------------------------------
  associate(lo => blkLimits(LOW,  :), &
            hi => blkLimits(HIGH, :))
     do       k = lo(KAXIS), hi(KAXIS)
        do    j = lo(JAXIS), hi(JAXIS)
           do i = lo(IAXIS), hi(IAXIS)
              U(DENS_VAR,i,j,k) = U(DENS_VAR,i,j,k)*hy_dref
              U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k)*hy_eref
              U(PRES_VAR,i,j,k) = U(PRES_VAR,i,j,k)*hy_pref
              U(VELX_VAR:VELZ_VAR,i,j,k) = U(VELX_VAR:VELZ_VAR,i,j,k)*hy_vref
           end do
        end do
     end do
  end associate
  !! ---------------------------------------------------------------------------
  if(dtCheck <= 0.0) call Driver_abort("[Hydro]: Computed dt is not positive! Aborting!")
  return

End Subroutine Hydro_computeDt


