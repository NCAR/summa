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

MODULE summa_type
! used to define master summa data structure
! *****************************************************************************
! * higher-level derived data types
! *****************************************************************************
USE nrtype          ! variable types, etc.
USE data_types,only:&
                    ! no spatial dimension
                    var_i,               & ! x%var(:)            (i4b)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    ! no variable dimension
                    hru_i,               & ! x%hru(:)            (i4b)
                    hru_d,               & ! x%hru(:)            (dp)
                    gru_i,               & ! x%gru(:)%hru(:)     (i4b)
                    gru_d,               & ! x%gru(:)%hru(:)     (dp)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)     (i4b)
                    gru_double,          & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_int8,        & ! x%gru(:)%hru(:)%var(:)     (i8b)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec,   & ! x%gru(:)%hru(:)%var(:)%dat (dp)
                    ! gru+hru+z dimension
                    gru_hru_z_vLookup      ! x%gru(:)%hru(:)%z(:)%var(:)%lookup(:)  (dp)
implicit none
private

! ************************************************************************
! * master summa data type
! *****************************************************************************
type, public :: summa1_type_dec
    ! define the lookup tables
    type(gru_hru_z_vLookup)          :: lookupStruct               ! x%gru(:)%hru(:)%z(:)%var(:)%lookup(:) -- lookup tables
    
    ! define the statistics structures
    type(gru_hru_doubleVec)          :: forcStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
    type(gru_hru_doubleVec)          :: progStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
    type(gru_hru_doubleVec)          :: diagStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
    type(gru_hru_doubleVec)          :: fluxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
    type(gru_hru_doubleVec)          :: indxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model indices
    type(gru_doubleVec)              :: bvarStat                   ! x%gru(:)%var(:)%dat        -- basin-average variabl

    ! define the primary data structures (scalars)
    type(var_i)                      :: timeStruct                 ! x%var(:)                   -- model time data
    type(gru_hru_double)             :: forcStruct                 ! x%gru(:)%hru(:)%var(:)     -- model forcing data
    type(gru_hru_double)             :: attrStruct                 ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
    type(gru_hru_int)                :: typeStruct                 ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
    type(gru_hru_int8)               :: idStruct                   ! x%gru(:)%hru(:)%var(:)     --

    ! define the primary data structures (variable length vectors)
    type(gru_hru_intVec)             :: indxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model indices
    type(gru_hru_doubleVec)          :: mparStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
    type(gru_hru_doubleVec)          :: progStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
    type(gru_hru_doubleVec)          :: diagStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
    type(gru_hru_doubleVec)          :: fluxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes

    ! define the basin-average structures
    type(gru_double)                 :: bparStruct                 ! x%gru(:)%var(:)            -- basin-average parameters
    type(gru_doubleVec)              :: bvarStruct                 ! x%gru(:)%var(:)%dat        -- basin-average variables

    ! define the ancillary data structures
    type(gru_hru_double)             :: dparStruct                 ! x%gru(:)%hru(:)%var(:)     -- default model parameters

    ! define the run-time variables
    type(gru_i)                      :: computeVegFlux             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
    type(gru_d)                      :: dt_init                    ! used to initialize the length of the sub-step for each HRU
    type(gru_d)                      :: upArea                     ! area upslope of each HRU

    ! define miscellaneous variables
    integer(i4b)                     :: summa1open                 ! flag to define if the summa file is open??
    integer(i4b)                     :: numout                     ! number of output variables??
    real(rkind)                         :: ts                         ! model time step ??
    integer(i4b)                     :: nGRU                       ! number of grouped response units
    integer(i4b)                     :: nHRU                       ! number of global hydrologic response units
    integer(i4b)                     :: hruCount                   ! number of local hydrologic response units
    real(rkind),dimension(12)           :: greenVegFrac_monthly       ! fraction of green vegetation in each month (0-1)
    character(len=256)               :: summaFileManagerFile       ! path/name of file defining directories and files

end type summa1_type_dec

END MODULE summa_type
