!!****if* source/physics/sourceTerms/Burn/BurnIntegrate/Ma28Stub/bn_ma28Stub
!!
!! NAME
!!
!!  bn_ma28Stub
!!
!! DESCRIPTION
!!
!!  Abort stubs for the three Harwell MA28 sparse-solver entry points that the
!!  Bader-Deuflhard and Rosenbrock network integrators reference: ma28ad,
!!  ma28bd and ma28cd.
!!
!!  Those integrators keep the MA28 and GIFT linear-algebra paths in a single
!!  source file each (bn_baderMa28.F90, bn_rosenMa28.F90).  Both paths are
!!  therefore compiled in and both sets of symbols must resolve at link time,
!!  even though bn_burner selects between them at runtime via the "algebra"
!!  runtime parameter (1 = MA28, 2 = GIFT).  MA28 is Harwell (HSL) code that
!!  Flash-X does not distribute -- lib/ma28 holds only a TYPE EXTERNAL
!!  declaration -- so a build without an external MA28 library fails to link
!!  purely because of that dead branch.
!!
!!  These stubs satisfy the linker.  Reaching one at runtime means algebra=1
!!  was requested without a real MA28 library, so they abort with a message
!!  pointing at algebra=2 rather than returning silently and producing
!!  meaningless abundances.
!!
!!  Compiled only when the setup variable Gift is set.
!!
!! NOTES
!!
!!  Dummy-argument declarations mirror the MA28 interfaces as called from the
!!  integrators.  The routines abort before touching any argument.
!!
!!***

subroutine ma28ad(n, nz, a, licn, irn, lirn, icn, u, ikeep, iw, w, iflag)

   use Driver_interface, ONLY: Driver_abort

   implicit none

   integer, intent(IN)    :: n, nz, licn, lirn
   real,    intent(INOUT) :: a(licn)
   integer, intent(INOUT) :: irn(lirn), icn(licn)
   real,    intent(IN)    :: u
   integer, intent(INOUT) :: ikeep(n, 5), iw(n, 8)
   real,    intent(INOUT) :: w(n)
   integer, intent(INOUT) :: iflag

   call Driver_abort("[ma28ad] MA28 is not available in this build. &
        &Set the runtime parameter algebra = 2 to use the GIFT solver, &
        &or rebuild against an external MA28 library.")

end subroutine ma28ad


subroutine ma28bd(n, nz, a, licn, ivect, jvect, icn, ikeep, iw, w, iflag)

   use Driver_interface, ONLY: Driver_abort

   implicit none

   integer, intent(IN)    :: n, nz, licn
   real,    intent(INOUT) :: a(licn)
   integer, intent(IN)    :: ivect(nz), jvect(nz), icn(licn)
   integer, intent(INOUT) :: ikeep(n, 5), iw(n, 5)
   real,    intent(INOUT) :: w(n)
   integer, intent(INOUT) :: iflag

   call Driver_abort("[ma28bd] MA28 is not available in this build. &
        &Set the runtime parameter algebra = 2 to use the GIFT solver, &
        &or rebuild against an external MA28 library.")

end subroutine ma28bd


subroutine ma28cd(n, a, licn, icn, ikeep, rhs, w, mtype)

   use Driver_interface, ONLY: Driver_abort

   implicit none

   integer, intent(IN)    :: n, licn, mtype
   real,    intent(IN)    :: a(licn)
   integer, intent(IN)    :: icn(licn), ikeep(n, 5)
   real,    intent(INOUT) :: rhs(n)
   real,    intent(INOUT) :: w(n)

   call Driver_abort("[ma28cd] MA28 is not available in this build. &
        &Set the runtime parameter algebra = 2 to use the GIFT solver, &
        &or rebuild against an external MA28 library.")

end subroutine ma28cd
