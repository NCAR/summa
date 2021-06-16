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

program summa_driver
! driver program for summa simulations
! *****************************************************************************
! * use desired modules
! *****************************************************************************
! data types
USE nrtype                                                  ! variable types, etc.
USE summa_type, only: summa1_type_dec                       ! master summa data type
! subroutines and functions: model setup
USE summa_init, only: summa_initialize                      ! used to allocate/initialize summa data structures
USE summa_setup, only: summa_paramSetup                     ! used to initialize parameter data structures (e.g. vegetation and soil parameters)
USE summa_restart, only: summa_readRestart                  ! used to read restart data and reset the model state
! subroutines and functions: model simulation
USE summa_forcing, only: summa_readForcing                  ! used to read forcing data
USE summa_modelRun, only: summa_runPhysics                  ! used to run the summa physics for one time step
USE summa_writeOutput, only: summa_writeOutputFiles         ! used to write the summa output files
! utility functions
USE summa_util, only: stop_program                          ! used to stop the summa program (with errors)
USE summa_util, only: handle_err                            ! used to process errors
! global data
USE globalData, only: numtim                                ! number of model time steps
implicit none

! *****************************************************************************
! * variable definitions
! *****************************************************************************
! define the master summa data structure
type(summa1_type_dec), allocatable :: summa1_struc(:)
! define parameters for the model simulation
integer(i4b), parameter            :: n=1                        ! number of instantiations
! define timing information
integer(i4b)                       :: modelTimeStep              ! index of model time step
! error control
integer(i4b)                       :: err=0                      ! error code
character(len=1024)                :: message=''                 ! error message

! *****************************************************************************
! * preliminaries
! *****************************************************************************

! allocate space for the master summa structure
allocate(summa1_struc(n), stat=err)
if(err/=0) call stop_program(1, 'problem allocating master summa structure')

! *****************************************************************************
! * model setup/initialization
! *****************************************************************************

! declare and allocate summa data structures and initialize model state to known values
call summa_initialize(summa1_struc(n), err, message)
call handle_err(err, message)

! initialize parameter data structures (e.g. vegetation and soil parameters)
call summa_paramSetup(summa1_struc(n), err, message)
call handle_err(err, message)

! read restart data and reset the model state
call summa_readRestart(summa1_struc(n), err, message)
call handle_err(err, message)

! *****************************************************************************
! * model simulation
! *****************************************************************************
! loop through time
do modelTimeStep=1,numtim

 ! read model forcing data
 call summa_readForcing(modelTimeStep, summa1_struc(n), err, message)
 call handle_err(err, message)
 ! run the summa physics for one time step
 call summa_runPhysics(modelTimeStep, summa1_struc(n), err, message)
 call handle_err(err, message)

 ! write the model output
 call summa_writeOutputFiles(modelTimeStep, summa1_struc(n), err, message)
 call handle_err(err, message)
end do  ! looping through time

! successful end
call stop_program(0, 'finished simulation successfully.')
end program summa_driver
