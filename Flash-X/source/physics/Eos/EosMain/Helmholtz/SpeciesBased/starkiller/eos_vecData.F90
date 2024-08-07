!!****if* source/physics/Eos/EosMain/Helmholtz/starkiller/eos_vecData
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
!!  eos_vecData
!! 
!! SYNOPSIS
!!
!!  use eos_vecData
!!
!! DESCRIPTION
!!
!!   Contains arrays for EOS Helmholtz computation.
!!   Incorporates both fixed block size and non-fixed blocksize computations
!!
!! ARGUMENTS
!!
!!
!!*** 

module eos_vecData

#include "Simulation.h"


#ifdef FIXEDBLOCKSIZE
  integer, parameter :: NROWMAX = GRID_IHI_GC*GRID_JHI_GC*GRID_KHI_GC*MAXBLOCKS
  !..thermodynamic and composition inputs
  real, dimension(NROWMAX), save ::  tempRow,denRow,                    &
       &                             abarRow,zbarRow

  !..totals and their derivatives
  real, dimension(NROWMAX), save ::  ptotRow,dptRow,             &
       &                             etotRow,detRow, stotRow            
! UNUSED       &                             dpaRow,dpzRow,                      &
! UNUSED       &                             deaRow,dezRow,                      &
! UNUSED       &                             dsaRow,dszRow
  real, dimension(NROWMAX), save ::  dedRow, dstRow,dsdRow,dpdRow            
  real, dimension(NROWMAX), save ::  deaRow, dezRow  !Calhoun            

  !..radiation contributions 
  !UNUSED  real, dimension(NROWMAX), save ::  pradRow, eradRow, sradRow

  !..ion contributions 
  !UNUSED  real, dimension(NROWMAX), save ::  pionRow,  eionRow,sionRow, xniRow


  !..electron-positron contributions -- most UNUSED and REMOVED
  real, dimension(NROWMAX), save :: pelRow, neRow, etaRow, detatRow
  !REMOVED real, dimension(NROWMAX), save:: etaeleRow, xneRow,                &
  !REMOVED                                  peleRow, seleRow, eeleRow,        &
  !REMOVED                                  pposRow, sposRow, eposRow,        &
  !REMOVED                                  dpeptRow, dseptRow, deeptRow,     &
  !REMOVED                                  dpepdRow, dsepdRow, deepdRow
  !REMOVED                                  dpepaRow,dpepzRow,                   &
  !REMOVED                                  deepaRow,deepzRow,                   &
  !REMOVED                                  dsepaRow,dsepzRow,                   &
  !REMOVED                                  dxnetRow,dxnedRow,           &
  !REMOVED                                  dxneaRow,dxnezRow,         &
  !REMOVED                                  detatRow,detadRow,                   &
  !REMOVED                                  detaaRow,detazRow

  !..coulomb contributions
  !REMOVED   real, dimension(NROWMAX), save ::   pcouRow, ecouRow,                    &
  !REMOVED                                       scouRow, plasgRow


  !..thermodynamic consistency checks; maxwell relations 
  !REMOVED   real, dimension(NROWMAX), save ::   dseRow,dpeRow,dspRow

  !..derivative based quantities
  real, dimension(NROWMAX), save ::    gamcRow
  real, dimension(NROWMAX), save ::    cpRow,cvRow 

#else
  !! NONFIXEDBLOCKSIZE begins here
  !
  !..thermodynamic and composition inputs
  real, allocatable, dimension(:), save ::  tempRow,denRow,                    &
       &                             abarRow,zbarRow

  !..totals and their derivatives
  real, allocatable, dimension(:), save ::  ptotRow,dptRow,             &
       &                             etotRow,detRow, stotRow
! UNUSED       &                             dpaRow,dpzRow,                      &
! UNUSED       &                             deaRow,dezRow,                      &
! UNUSED       &                             dsaRow,dszRow
  real, allocatable, dimension(:), save ::  dsdRow, dstRow, dedRow,dpdRow
  real, allocatable, dimension(:), save ::  deaRow, dezRow ! DL following Calhoun
  !..radiation contributions 
  !UNUSED  real, allocatable, dimension(:), save ::  pradRow, eradRow, sradRow

  !..ion contributions 
  !UNUSED  real, allocatable, dimension(:), save ::  pionRow,  eionRow, sionRow, xniRow


  !..electron-positron contributions -- all UNUSED and REMOVED from allocation
  real, allocatable, dimension(:), save :: pelRow, neRow, etaRow, detatRow
  !REMOVED  real, allocatable, dimension(:), save ::   etaeleRow, xneRow                     &
  !REMOVED                                             peleRow, seleRow, eeleRow,            &
  !REMOVED                                             pposRow, sposRow, eposRow,            &
  !REMOVED                                             dpeptRow, dseptRow, deeptRow,         &
  !REMOVED                                             dpepdRow, dsepdRow, deepdRow
  !REMOVED                                             dpepaRow,dpepzRow,                   &
  !REMOVED                                             deepaRow,deepzRow,                   &
  !REMOVED                                             dsepaRow,dsepzRow,                   &
  !REMOVED                                             dxnetRow,dxnedRow,           &
  !REMOVED                                             dxneaRow,dxnezRow,           &
  !REMOVED                                             detatRow,detadRow,                   &
  !REMOVED                                             detaaRow,detazRow


  !..coulomb contributions -- all UNUSED and REMOVED from allocation
  !REMOVED   real, allocatable, dimension(:), save ::   pcouRow, ecouRow,                    &
  !REMOVED        scouRow, plasgRow


  !..thermodynamic consistency checks; maxwell relations -- all UNUSED and REMOVED from allocation
  !REMOVED   real, allocatable, dimension(:), save ::   dseRow,dpeRow,dspRow

  !..derivative based quantities
  real, allocatable, dimension(:), save ::   gamcRow
  real, allocatable, dimension(:), save ::    cpRow,cvRow

#endif

!All variables in this module must be threadprivate!!!
!$omp threadprivate(tempRow, denRow, etotRow, abarRow, zbarRow, gamcRow, ptotRow, &
!$omp deaRow, dezRow, detRow, dptRow, dpdRow, dedRow, pelRow, neRow, etaRow, detatRow, &
!$omp cvRow, cpRow, dstRow, dsdRow, stotRow)

end module eos_vecData

