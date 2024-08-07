subroutine gr_markRefineDerefineCallback(lev, tags, time, tagval, clearval) bind(c)
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
   use iso_c_binding
   use amrex_fort_module,      ONLY : wp => amrex_real
   use amrex_box_module,       ONLY : amrex_box
   use amrex_tagbox_module,    ONLY : amrex_tagboxarray
   use amrex_multifab_module,  ONLY : amrex_mfiter, &
                                      amrex_mfiter_build, &
                                      amrex_mfiter_destroy, &
                                      amrex_multifab_build
 
   use Grid_tile,         ONLY : Grid_tile_t
   use gr_physicalMultifabs,   ONLY : unk

   implicit none
 
#include "Simulation.h"
#include "constants.h"

   integer,           intent(IN), value :: lev
   type(c_ptr),       intent(in), value :: tags 
   real(wp),          intent(in), value :: time
   character(c_char), intent(in), value :: tagval
   character(c_char), intent(in), value :: clearval

   type(amrex_tagboxarray) :: tag
   type(amrex_mfiter)      :: mfi                                                             
   type(amrex_box)         :: box

   character(kind=c_char), contiguous, pointer :: tagData(:,:,:,:)

   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Started on level ", lev + 1
   
   tag = tags

   call amrex_mfiter_build(mfi, unk(lev), tiling=.FALSE.)                             
   do while(mfi%next())
      box   = mfi%tilebox()

      associate(lo => box%lo)
        tagData(lo(1):, lo(2):, lo(3):, 1:) => tag%dataptr(mfi)
      end associate

      tagData(:, :, :, :) = clearval
   end do
   call amrex_mfiter_destroy(mfi)

   write(*,'(A,A,I2)') "[gr_markRefineDerefineCallback]", &
                       "      Finished on level ", lev + 1
end subroutine gr_markRefineDerefineCallback

