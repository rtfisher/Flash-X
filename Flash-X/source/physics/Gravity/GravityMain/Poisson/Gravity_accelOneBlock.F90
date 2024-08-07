!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_accelOneBlock
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
!!  Gravity_accelOneBlock
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneBlock(Grid_tile_t, intent(in) :: tileDesc, 
!!                        integer, intent(in) :: ngcellcomp,
!!                        real(:,:,:,:)),intent(out) :: gvec, 
!!                        integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration for this block.  Include ngcell layers outside
!!  block interior.
!!
!!  This routine computes the gravitational acceleration for
!!  zones in a given block. First-order
!!  finite-volume differencing is used everywhere.  It is assumed
!!  here that the requisite number of guard cells have peen appropriately
!!  filled for the variable containting the gravitational potential.
!!
!!  Dean Townsley 2008
!!  Contributed to Flash Center at the University of Chicago 2008
!!
!! ARGUMENTS
!!
!!  block            -  The local block metadata
!!  gvec(:,:,:,:)   -  Array to receive gravitational acceleration
!!                        as as NDIM-dimensional vector.  It is assumed
!!                        the the space provided is the size of the block
!!                        plus all guard cells.  The first index is the vector
!!                        component and the latter are cell indices.
!!  ngcellcomp         -  Number of layers outside of block interior to
!!                        compute gravity
!!  potentialIndex     -  if specified,  Variable # to take as potential.
!!                        Default is GPOT_VAR for the potential stored in the
!!                        gpot slot of unk, which should correspond to the
!!                        potential at the current timestep.
!!
!!
!!***

subroutine Gravity_accelOneBlock ( tileDesc, ngcellcomp, gvec, potentialIndex)


  use Grid_interface, ONLY : Grid_getGeometry
  use Driver_interface, ONLY : Driver_abort
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

#include "Simulation.h"
#include "constants.h"

  type(Grid_tile_t) :: tileDesc
  integer, intent(in)                    :: ngcellcomp
  real, dimension(:,:,:,:), intent(out)  :: gvec
  integer, intent(in),optional           :: potentialIndex

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(LOW:HIGH,MDIM)        :: blkLimits, blkLimitsGC
  real, dimension(MDIM)             :: delta, inv_2delta
  integer         :: sizeI, sizeJ, sizeK
  integer         :: potVar, geom, istat
  real, dimension(:,:,:), allocatable :: gpot(:,:,:)
  real, dimension(:), allocatable     :: inv_r, inv_sintheta
  integer :: i,j,k
  
  call Driver_abort("[Gravity_accelOneBlock] Implement for tiling")

  !==================================================

  ! If a variable index is explicitly specified, assume that as the potential
  ! otherwise use the default current potential GPOT_VAR  
  if(present(potentialIndex)) then
     potVar=potentialIndex
  else
     potVar=GPOT_VAR
  end if

  ! get block data and block info
  call tileDesc%getDataPtr(solnData, CENTER)
  call tileDesc%deltas(delta)
  inv_2delta(1:NDIM) = 0.5/delta(1:NDIM)
  blkLimits = tileDesc%limits
  blkLimitsGC = tileDesc%blkLimitsGC
  call Grid_getGeometry(geom)

  ! DEV: FIXME Should this be tile and halo or grown tile?
  call Driver_abort("[Gravity_accelOneBlock] Implement for tiling")
  sizeI = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeJ = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeK = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  ! pull out potential for better cache locality
  allocate( gpot(sizeI,sizeJ,sizeK), STAT=istat)
  if (istat/=0) call Driver_abort("unable to allocate gpot in Gravity_accelOneBlock")
  ! DEV: TODO For tiling, we might need to write this as an explicit
  !           loop nest
  gpot(:,:,:) = solnData(potVar, :,:,:)

  call tileDesc%releaseDataPtr(solnData, CENTER)

  ! now calculate gravity, depending on geometry
  select case (geom)
  case (CARTESIAN)

     do k = blkLimits(LOW,KAXIS)-K3D*ngcellcomp, blkLimits(HIGH,KAXIS)+K3D*ngcellcomp
        do j = blkLimits(LOW,JAXIS)-K2D*ngcellcomp, blkLimits(HIGH,JAXIS)+K2D*ngcellcomp
           do i = blkLimits(LOW,IAXIS)-ngcellcomp, blkLimits(HIGH,IAXIS)+ngcellcomp

              gvec(IAXIS,i,j,k) = (gpot(i+1,j,k)-gpot(i-1,j,k))*inv_2delta(IAXIS)
              if (NDIM>=2) gvec(JAXIS,i,j,k) = (gpot(i,j+1,k)-gpot(i,j-1,k))*inv_2delta(JAXIS)
              if (NDIM==3) gvec(KAXIS,i,j,k) = (gpot(i,j,k+1)-gpot(i,j,k-1))*inv_2delta(KAXIS)

           enddo
        enddo
     enddo

  case (CYLINDRICAL)

     ! need radius if theta axis is included
     if (NDIM==3) then
        allocate(inv_r(sizeI), STAT=istat)
        if (istat/=0) call Driver_abort("unable to allocate inv_r in Gravity_accelOneBlock")
        inv_r(:) = 1.0/inv_r(:)
     endif
     do k = blkLimits(LOW,KAXIS)-K3D*ngcellcomp, blkLimits(HIGH,KAXIS)+K3D*ngcellcomp
        do j = blkLimits(LOW,JAXIS)-K2D*ngcellcomp, blkLimits(HIGH,JAXIS)+K2D*ngcellcomp
           do i = blkLimits(LOW,IAXIS)-ngcellcomp, blkLimits(HIGH,IAXIS)+ngcellcomp
           
              gvec(IAXIS,i,j,k) = (gpot(i+1,j,k)-gpot(i-1,j,k))*inv_2delta(IAXIS)
              if (NDIM>=2) gvec(JAXIS,i,j,k) = (gpot(i,j+1,k)-gpot(i,j-1,k))*inv_2delta(JAXIS)
              if (NDIM==3) gvec(KAXIS,i,j,k) = inv_r(i)*(gpot(i,j,k+1)-gpot(i,j,k-1))*inv_2delta(KAXIS)

           enddo
        enddo
     enddo
     if (NDIM==3) deallocate(inv_r)

  case (SPHERICAL)

     ! need radius if 2d or 3d
     if (NDIM>=2) then
        allocate(inv_r(sizeI), STAT=istat)
        if (istat/=0) call Driver_abort("unable to allocate inv_r in Gravity_accelOneBlock")
        inv_r(:) = 1.0/inv_r(:)
     endif
     ! need sin_theta if 3d
     if (NDIM==3) then
        allocate(inv_sintheta(sizeJ), STAT=istat)
        if (istat/=0) call Driver_abort("unable to allocate inv_sintheta in Gravity_accelOneBlock")
        inv_sintheta(:) = 1.0/inv_sintheta(:)
     endif
     do k = blkLimits(LOW,KAXIS)-K3D*ngcellcomp, blkLimits(HIGH,KAXIS)+K3D*ngcellcomp
        do j = blkLimits(LOW,JAXIS)-K2D*ngcellcomp, blkLimits(HIGH,JAXIS)+K2D*ngcellcomp
           do i = blkLimits(LOW,IAXIS)-ngcellcomp, blkLimits(HIGH,IAXIS)+ngcellcomp
           
              gvec(IAXIS,i,j,k) = (gpot(i+1,j,k)-gpot(i-1,j,k))*inv_2delta(IAXIS)
              if (NDIM>=2) gvec(JAXIS,i,j,k) = inv_r(i)*(gpot(i,j+1,k)-gpot(i,j-1,k))*inv_2delta(JAXIS)
              if (NDIM==3) gvec(KAXIS,i,j,k) = inv_r(i)*inv_sintheta(j)* &
                                                    (gpot(i,j,k+1)-gpot(i,j,k-1))*inv_2delta(KAXIS)

           enddo
        enddo
     enddo
     if (NDIM==3) deallocate(inv_sintheta)
     if (NDIM>=2) deallocate(inv_r)

  case default
     call Driver_abort("unhandled geometry in Gravity_accelOneBlock")
  end select

  deallocate(gpot)
  
  return
   
end subroutine Gravity_accelOneBlock
