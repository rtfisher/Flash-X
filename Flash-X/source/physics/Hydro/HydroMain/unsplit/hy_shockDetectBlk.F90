!!****if* source/physics/Hydro/HydroMain/unsplit/hy_shockDetectBlk
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
!!  hy_shockDetectBlk
!!
!! SYNOPSIS
!!
!!  hy_shockDetectBlk( integer (IN) :: blockID )
!!
!! DESCRIPTION
!!
!!  This routine detects strongly compressive motions in simulation
!!  by calculating undivided pressure gradients and divergence of
!!  velocity fields. Two parameters beta and delta have been set 
!!  to detect strong shocks. If such shocks exist then the unsplit
!!  scheme applies its robust flux differencings using limited slopes
!!  in data reconstruction step (see hy_dataReconstruct.F90).
!!  Different shock strengths can also be detected by lowering/increasing
!!  beta and delta values.
!!
!! ARGUMENTS
!!
!!  blockID  - local block ID
!!
!! REFERENCE 
!!
!!  Balsara and Spicer, JCP, 149:270--292, 1999.
!!
!! NOTES
!!
!!  The set of variables read from Uin does not overlap with the set
!!  of variables updated in Uout. As long as it stays that way, there
!!  is no problem with violationg Fortran's no-aliasing assumptions,
!!  and thus no need to turn those dummy arguments into POINTERs
!!  or give them TARGET attributes.
!!
!!***

!!REORDER(4): Uin,Uout

Subroutine hy_shockDetectBlk(Uin,loI,blkLimitsGC,Uout,loO,blkLimits,del )


  use Hydro_data,        ONLY : hy_cfl, hy_cfl_original,&
                                hy_RiemannSolver,       &
                                hy_geometry,            &
                                hy_fallbackLowerCFL,    &
                                hy_shockLowerCFL
  use Driver_data,       ONLY : dr_nStep
  use Logfile_interface, ONLY : Logfile_open,Logfile_close
  use Driver_interface,  ONLY : Driver_abort

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer,dimension(LOW:HIGH,MDIM),INTENT(IN) :: blkLimits,blkLimitsGC
  integer, intent(IN)  :: loI(*), loO(*)
  real,    intent(IN)  :: UIN(loI(1):,loI(2):,loI(3):,loI(4):)  !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: UOUT(loO(1):,loO(2):,loO(3):,loO(4):) !CAPITALIZATION INTENTIONAL!
!  real, dimension(:,:,:,:),INTENT(IN) :: Uin
!  real,dimension(:,:,:,:),INTENT(INOUT) :: Uout
  real,dimension(MDIM),INTENT(IN) :: del
  !! -----------------------------------------------------

  integer :: i,j,k,ib,ie,jb,je,kb,ke
  logical :: SW1, SW2
  logical :: doLowerCFL
  integer :: k2,k3


  real :: divv,gradPx,gradPy,gradPz
  real :: minP,minC,beta,delta
  real :: cflMax
  real, dimension(:,:,:), allocatable :: Vc

  ib=blkLimitsGC(LOW,IAXIS)
  ie=blkLimitsGC(HIGH,IAXIS)  
  jb=blkLimitsGC(LOW,JAXIS)
  je=blkLimitsGC(HIGH,JAXIS)  
  kb=blkLimitsGC(LOW,KAXIS)
  ke=blkLimitsGC(HIGH,KAXIS)  

  ! Two parameters that can be adjusted to detect shocks
  ! with different strengths:
  ! (a) The lower the values the weaker shocks detected 
  !     (lower the values to detect more shock regions)
  ! (b) The larger the values the stronger shocks detected
  !     (increase the values to detect less shock regions)
  beta = 0.5 !0.5 !10. ! gradP
  delta= 0.1  !0.1 !2. ! divV
!!$  beta  = 0.1 !0.1
!!$  delta = 0.01

  k2=0
  k3=0
  ! Set dimensional indices
  if (NDIM > 1) k2=1
  if (NDIM > 2) k3=1

#ifdef DEBUG_UHD
!!$     print*,'_shockDetectBlk top: associated(Uin ) is',associated(Uin )
!!$     print*,'_shockDetectBlk top: associated(Uout) is',associated(Uout)
     print*,'_shockDetectBlk top: lbound(Uin ):',lbound(Uin )
     print*,'_shockDetectBlk top: ubound(Uin ):',ubound(Uin )
     print*,'_shockDetectBlk top: lbound(Uout):',lbound(Uout)
     print*,'_shockDetectBlk top: ubound(Uout):',ubound(Uout)
     print*,'_shockDetectBlk top: lbound(Uout(SHOK_VAR,ib:ie,jb:je,kb:ke)):',lbound(Uout(SHOK_VAR,ib:ie,jb:je,kb:ke))
     print*,'_shockDetectBlk top: ubound(Uout(SHOK_VAR,ib:ie,jb:je,kb:ke)):',ubound(Uout(SHOK_VAR,ib:ie,jb:je,kb:ke))
#endif

#ifdef SHOK_VAR
     Uout(SHOK_VAR,ib:ie,jb:je,kb:ke)=0.
#else
  if (hy_RiemannSolver == HYBR) then
     call Driver_abort&
          ("[hy_shockDetectBlk]: SHOK_VAR has not been defined for shock detection")
  endif
#endif


  !! Allocate a temporary cell-centered array for sound speed
  allocate(Vc(blkLimits(LOW,IAXIS)- 1:blkLimits(HIGH,IAXIS)+1,  &
              blkLimits(LOW,JAXIS)-k2:blkLimits(HIGH,JAXIS)+k2, &
              blkLimits(LOW,KAXIS)-k3:blkLimits(HIGH,KAXIS)+k3))

  !! Compute sound speed
  do k=blkLimits(LOW,KAXIS)-K3D,blkLimits(HIGH,KAXIS)+K3D
     do j=blkLimits(LOW,JAXIS)-K2D,blkLimits(HIGH,JAXIS)+K2D
        do i=blkLimits(LOW,IAXIS)-1,blkLimits(HIGH,IAXIS)+1
           Vc(i,j,k) = sqrt(Uin(GAMC_VAR,i,j,k)*Uin(PRES_VAR,i,j,k)/Uin(DENS_VAR,i,j,k))
        enddo
     enddo
  enddo

  ! Note: it is the job of the caller to revert back to the original cfl, if desired
  doLowerCFL = .FALSE.

  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
#ifdef CFL_VAR
           if (hy_shockLowerCFL .OR. .NOT. hy_fallbackLowerCFL) &
                Uout(CFL_VAR,i,j,k) = hy_cfl
#endif

           ! initialize switch values
           SW1 = .false.
           SW2 = .false.

           minP = minval(Uin(PRES_VAR,i-1:i+1,j-k2:j+k2,k-k3:k+k3))
           minC = minval(        Vc(i-1:i+1,j-k2:j+k2,k-k3:k+k3))

           !! We do not need to include non-Cartesian geom factors here.
           !! Undivided divV
           divv =        Uin(VELX_VAR,i+1,j,  k  ) - Uin(VELX_VAR,i-1,j,  k  )
#if NDIM > 1
           divv = divv + Uin(VELY_VAR,i,  j+1,k  ) - Uin(VELY_VAR,i,  j-1,k  )
#if NDIM == 3
           divv = divv + Uin(VELZ_VAR,i,  j,  k+1) - Uin(VELZ_VAR,i,  j,  k-1)
#endif
#endif
           divv = 0.5*divv  

           !! Undivided grad pres
           gradPx = 0.5*(Uin(PRES_VAR,i+1,j,  k  ) - Uin(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.
           gradPz = 0.
#if NDIM > 1
           gradPy = 0.5*(Uin(PRES_VAR,i,  j+1,k  ) - Uin(PRES_VAR,i,  j-1,k  ))
#if NDIM == 3
           gradPz = 0.5*(Uin(PRES_VAR,i,  j,  k+1) - Uin(PRES_VAR,i,  j,  k-1))
#endif
#endif
           if ( abs(gradPx)+abs(gradPy)+abs(gradPz) .ge. beta*minP ) then
              SW1 = .true.
           endif

           if (-delta*minC .ge. divv) then
              SW2 = .true.
           endif


           if (SW1 .and. SW2) then
              ! Set SHOCK_VAR to 1.0 if a shock is detected. 
              ! One use is for a local hybrid method in the Hydro unit which
              ! applies (a diffusive) HLL solver when SHOK_VAR = 1.
#ifdef SHOK_VAR
              Uout(SHOK_VAR,i,j,k) = 1.
#endif

              if (hy_shockLowerCFL) then
                 ! if lowering of the CFL factor within shocks is requested...

                 doLowerCFL = .TRUE.

#ifdef CFL_VAR
                 if (NDIM == 1) then
                    Uout(CFL_VAR,i,j,k) = min(Uout(CFL_VAR,i,j,k),0.60)
                 elseif (NDIM == 2) then
                    Uout(CFL_VAR,i,j,k) = min(Uout(CFL_VAR,i,j,k),0.45)
                 elseif (NDIM == 3) then
                    Uout(CFL_VAR,i,j,k) = min(Uout(CFL_VAR,i,j,k),0.25)
                 endif
#endif
              endif ! endif (hy_shockLowerCFL) then

           endif !endif (SW1 .and. SW2) then

        enddo !enddo i-loop
     enddo !enddo j-loop
  enddo !enddo k-loop

  ! Release block pointer

  ! Deallocate sound speed array
  deallocate(Vc)

  if (doLowerCFL) then
     if (NDIM == 1) then
        cflMax = 0.60
     elseif (NDIM == 2) then
        cflMax = 0.45
     elseif (NDIM == 3) then
        cflMax = 0.25
     endif

!!$     !$omp critical (Update_cfl)
!!$     hy_cfl = min(hy_cfl,cflMax)
!!$     !$omp end critical (Update_cfl)
!!$
     !$omp atomic 
     hy_cfl = min(hy_cfl,cflMax)

  endif

End Subroutine hy_shockDetectBlk
