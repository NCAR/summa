! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module summa_modelRun
! calls the model physics

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing double precision number

! named variables
USE globalData,only:yes,no           ! .true. and .false.
USE var_lookup,only:iLookTIME        ! named variables for time data structure
USE var_lookup,only:iLookDIAG        ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookINDEX       ! look-up values for local column index variables
USE summa_util,only:handle_err

! safety: set private unless specified otherwise
implicit none
private
public::summa_runPhysics
contains

 ! calls the model physics
 subroutine summa_runPhysics(modelTimeStep, summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nrtype                                                     ! variable types, etc.
 USE summa_type, only:summa1_type_dec                           ! master summa data type
 ! subroutines and functions
 USE nr_utility_module,only:indexx                              ! sort vectors in ascending order
 USE vegPhenlgy_module,only:vegPhenlgy                          ! module to compute vegetation phenology
 USE run_oneGRU_module,only:run_oneGRU                          ! module to run for one GRU
 USE time_utils_module,only:elapsedSec                          ! calculate the elapsed time
 ! global data
 USE globalData,only:gru_struc                                  ! gru-hru mapping structures
 USE globalData,only:model_decisions                            ! model decision structure
 USE globalData,only:startPhysics,endPhysics                    ! date/time for the start and end of the initialization
 USE globalData,only:elapsedPhysics                             ! elapsed time for the initialization
 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 integer(i4b),intent(in)               :: modelTimeStep         ! time step index
 type(summa1_type_dec),intent(inout)   :: summa1_struc          ! master summa data structure
 integer(i4b),intent(out)              :: err                   ! error code
 character(*),intent(out)              :: message               ! error message
 ! ---------------------------------------------------------------------------------------
 ! local variables: general
 character(LEN=256)                    :: cmessage              ! error message of downwind routine
 integer(i4b)                          :: iHRU                  ! HRU index
 integer(i4b)                          :: iGRU,jGRU,kGRU        ! GRU indices
 ! local variables: veg phenology
 logical(lgt)                          :: computeVegFluxFlag    ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(rkind)                              :: notUsed_canopyDepth   ! NOT USED: canopy depth (m)
 real(rkind)                              :: notUsed_exposedVAI    ! NOT USED: exposed vegetation area index (m2 m-2)
 ! local variables: parallelize the model run
 integer(i4b), allocatable             :: ixExpense(:)          ! ranked index GRU w.r.t. computational expense
 integer(i4b), allocatable             :: totalFluxCalls(:)     ! total number of flux calls for each GRU
 ! local variables: timing information
 integer*8                             :: openMPstart,openMPend ! time for the start of the parallelization section
 integer*8, allocatable                :: timeGRUstart(:)       ! time GRUs start
 real(rkind),  allocatable                :: timeGRUcompleted(:)   ! time required to complete each GRU
 real(rkind),  allocatable                :: timeGRU(:)            ! time spent on each GRU
 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&

  ! primary data structures (scalars)
  timeStruct           => summa1_struc%timeStruct          , & ! x%var(:)                   -- model time data
  forcStruct           => summa1_struc%forcStruct          , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
  attrStruct           => summa1_struc%attrStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
  typeStruct           => summa1_struc%typeStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
  idStruct             => summa1_struc%idStruct            , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU

  ! primary data structures (variable length vectors)
  indxStruct           => summa1_struc%indxStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  mparStruct           => summa1_struc%mparStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
  progStruct           => summa1_struc%progStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes

  ! basin-average structures
  bparStruct           => summa1_struc%bparStruct          , & ! x%gru(:)%var(:)            -- basin-average parameters
  bvarStruct           => summa1_struc%bvarStruct          , & ! x%gru(:)%var(:)%dat        -- basin-average variables

  ! run time variables
  greenVegFrac_monthly => summa1_struc%greenVegFrac_monthly, & ! fraction of green vegetation in each month (0-1)
  computeVegFlux       => summa1_struc%computeVegFlux      , & ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  dt_init              => summa1_struc%dt_init             , & ! used to initialize the length of the sub-step for each HRU
  nGRU                 => summa1_struc%nGRU                  & ! number of grouped response units

 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_runPhysics/'

 ! *******************************************************************************************
 ! *** initialize computeVegFlux (flag to indicate if we are computing fluxes over vegetation)
 ! *******************************************************************************************

 ! if computeVegFlux changes, then the number of state variables changes, and we need to reoranize the data structures
 if(modelTimeStep==1)then
  do iGRU=1,nGRU
   do iHRU=1,gru_struc(iGRU)%hruCount

    ! get vegetation phenology
    ! (compute the exposed LAI and SAI and whether veg is buried by snow)
    call vegPhenlgy(&
                    ! input/output: data structures
                    model_decisions,                & ! intent(in):    model decisions
                    typeStruct%gru(iGRU)%hru(iHRU), & ! intent(in):    type of vegetation and soil
                    attrStruct%gru(iGRU)%hru(iHRU), & ! intent(in):    spatial attributes
                    mparStruct%gru(iGRU)%hru(iHRU), & ! intent(in):    model parameters
                    progStruct%gru(iGRU)%hru(iHRU), & ! intent(in):    model prognostic variables for a local HRU
                    diagStruct%gru(iGRU)%hru(iHRU), & ! intent(inout): model diagnostic variables for a local HRU
                    ! output
                    computeVegFluxFlag,             & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                    notUsed_canopyDepth,            & ! intent(out): NOT USED: canopy depth (m)
                    notUsed_exposedVAI,             & ! intent(out): NOT USED: exposed vegetation area index (m2 m-2)
                    err,cmessage)                     ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! save the flag for computing the vegetation fluxes
    if(computeVegFluxFlag)      computeVegFlux%gru(iGRU)%hru(iHRU) = yes
    if(.not.computeVegFluxFlag) computeVegFlux%gru(iGRU)%hru(iHRU) = no

    ! define the green vegetation fraction of the grid box (used to compute LAI)
    diagStruct%gru(iGRU)%hru(iHRU)%var(iLookDIAG%scalarGreenVegFraction)%dat(1) = greenVegFrac_monthly(timeStruct%var(iLookTIME%im))

   end do  ! looping through HRUs
  end do  ! looping through GRUs
 end if  ! if the first time step

 ! ****************************************************************************
 ! *** model simulation
 ! ****************************************************************************

 ! initialize the start of the physics
 call date_and_time(values=startPhysics)

 ! ----- rank the GRUs in terms of their anticipated computational expense -----

 ! estimate computational expense based on persistence
 !  -- assume that that expensive GRUs from a previous time step are also expensive in the current time step

 ! allocate space for GRU timing
 allocate(totalFluxCalls(nGRU), timeGRU(nGRU), timeGRUstart(nGRU), timeGRUcompleted(nGRU), ixExpense(nGRU), stat=err)
 if(err/=0)then; message=trim(message)//'unable to allocate space for GRU timing'; return; endif
 timeGRU(:) = realMissing ! initialize because used for ranking

 ! compute the total number of flux calls from the previous time step
 do jGRU=1,nGRU
  totalFluxCalls(jGRU) = 0._rkind
  do iHRU=1,gru_struc(jGRU)%hruCount
   totalFluxCalls(jGRU) = totalFluxCalls(jGRU) + indxStruct%gru(jGRU)%hru(iHRU)%var(iLookINDEX%numberFluxCalc)%dat(1)
  end do
 end do

 ! get the indices that can rank the computational expense
 call indexx(timeGRU, ixExpense) ! ranking of each GRU w.r.t. computational expense
 ixExpense=ixExpense(nGRU:1:-1)  ! reverse ranking: now largest to smallest

 ! initialize the GRU count
 ! NOTE: this needs to be outside the parallel section so it is not reinitialized by different threads
 kGRU=0
 end associate summaVars

 ! initialize the time that the openMP section starts
 call system_clock(openMPstart)

 ! ----- use openMP directives to run GRUs in parallel -------------------------
 ! start of parallel section: define shared and private structure elements
  !$omp parallel  default(none) &
  !$omp          private(iGRU, jGRU) &  ! GRU indices are private for a given thread
  !$omp          shared(openMPstart, openMPend)   & ! access constant variables
  !$omp          shared(timeGRUstart, timeGRUcompleted, timeGRU, ixExpense, kGRU)  & ! time variables shared
  !$omp          shared(summa1_struc, gru_struc) &
  !$omp          private(err, cmessage)
 ! associate to elements in the data structur, gru_struce
 ! need to associate again for the parallelism to work
 summaVars2: associate(&

  ! primary data structures (scalars)
  timeStruct           => summa1_struc%timeStruct          , & ! x%var(:)                   -- model time data
  forcStruct           => summa1_struc%forcStruct          , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
  attrStruct           => summa1_struc%attrStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
  typeStruct           => summa1_struc%typeStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
  idStruct             => summa1_struc%idStruct            , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU

  ! primary data structures (variable length vectors)
  indxStruct           => summa1_struc%indxStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  mparStruct           => summa1_struc%mparStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
  progStruct           => summa1_struc%progStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes

  ! basin-average structures
  bparStruct           => summa1_struc%bparStruct          , & ! x%gru(:)%var(:)            -- basin-average parameters
  bvarStruct           => summa1_struc%bvarStruct          , & ! x%gru(:)%var(:)%dat        -- basin-average variables

  ! lookup table structure
  lookupStruct         => summa1_struc%lookupStruct        , & ! x%gru(:)%hru(:)%z(:)%var(:)%lookup    -- lookup-tables
  
  ! run time variables
  greenVegFrac_monthly => summa1_struc%greenVegFrac_monthly, & ! fraction of green vegetation in each month (0-1)
  computeVegFlux       => summa1_struc%computeVegFlux      , & ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  dt_init              => summa1_struc%dt_init             , & ! used to initialize the length of the sub-step for each HRU
  nGRU                 => summa1_struc%nGRU                  & ! number of grouped response units

 ) ! assignment to variables in the data structures

 !$omp do schedule(dynamic, 1)
 do jGRU=1,nGRU  ! loop through GRUs

  !----- process GRUs in order of computational expense -------------------------
  !$omp critical(setGRU)
  ! assign expensive GRUs to threads that enter first
  kGRU = kGRU+1
  iGRU = ixExpense(kGRU)
  ! get the time that the GRU started
  call system_clock( timeGRUstart(iGRU) )
  !$omp end critical(setGRU)

  !----- run simulation for a single GRU ----------------------------------------
  call run_oneGRU(&
                  ! model control
                  gru_struc(iGRU),              & ! intent(inout): HRU information for given GRU (# HRUs, #snow+soil layers)
                  dt_init%gru(iGRU)%hru,        & ! intent(inout): used to initialize the length of the sub-step for each HRU
                  computeVegFlux%gru(iGRU)%hru, & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                  ! data structures (input)
                  timeStruct%var,               & ! intent(in):    model time data
                  typeStruct%gru(iGRU),         & ! intent(in):    local classification of soil veg etc. for each HRU
                  idStruct%gru(iGRU),           & ! intent(in):    local classification of soil veg etc. for each HRU
                  attrStruct%gru(iGRU),         & ! intent(in):    local attributes for each HRU
                  lookupStruct%gru(iGRU),       & ! intent(in):    local lookup tables for each HRU
                  ! data structures (input-output)
                  mparStruct%gru(iGRU),         & ! intent(inout): local model parameters
                  indxStruct%gru(iGRU),         & ! intent(inout): model indices
                  forcStruct%gru(iGRU),         & ! intent(inout): model forcing data
                  progStruct%gru(iGRU),         & ! intent(inout): prognostic variables for a local HRU
                  diagStruct%gru(iGRU),         & ! intent(inout): diagnostic variables for a local HRU
                  fluxStruct%gru(iGRU),         & ! intent(inout): model fluxes for a local HRU
                  bvarStruct%gru(iGRU),         & ! intent(inout): basin-average variables
                  ! error control
                  err,cmessage)                   ! intent(out):   error control

  ! check errors
  call handle_err(err, cmessage)

  !----- save timing information ------------------------------------------------
  !$omp critical(saveTiming)
  ! save timing information
  call system_clock(openMPend)
  timeGRU(iGRU)          = real(openMPend - timeGRUstart(iGRU), kind(rkind))
  timeGRUcompleted(iGRU) = real(openMPend - openMPstart       , kind(rkind))
  !$omp end critical(saveTiming)

 end do  ! (looping through GRUs)
 !$omp end do
 end associate summaVars2
 !$omp end parallel

 ! identify the end of the physics
 call date_and_time(values=endPhysics)

 ! aggregate the elapsed time for the physics
 elapsedPhysics = elapsedPhysics + elapsedSec(startPhysics, endPhysics)

 ! deallocate space used to determine the GRU computational expense
 deallocate(totalFluxCalls, ixExpense, timeGRU, stat=err)
 if(err/=0)then; message=trim(message)//'unable to deallocate space for GRU timing'; return; endif

 ! end associate statements

 end subroutine summa_runPhysics

end module summa_modelRun
