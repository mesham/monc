!> Performs the actual time stepping over groups of components. Each group can be the whole (which is one call per
!! component per timestep) or column, which calls components for each column of the timestep. Groups are executed
!! sequentially in the order that they have been configured (which is already set up in the registry)
module timestepper_mod
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX
  use registry_mod, only : GROUP_TYPE_WHOLE, GROUP_TYPE_COLUMN, group_descriptor_type, get_ordered_groups, &
       execute_timestep_callbacks
  use optionsdatabase_mod, only : options_get_integer
  use registry_mod, only : is_component_enabled

  implicit none

#ifndef TEST_MODE
  private
#endif

  type(group_descriptor_type), dimension(:), allocatable :: group_descriptors !< Prefetched ordered group descriptors

  integer :: radiation_interval
  logical :: socrates_enabled

  public init_timestepper, timestep, finalise_timestepper
contains

  !> Initialises the timestepper by prefetching the groups in the order that they will be executed, this is for optimised
  !! execution in the timestep calls
  !! @param current_state The current model state
  subroutine init_timestepper(current_state)
    type(model_state_type), intent(inout) :: current_state

    call get_ordered_groups(group_descriptors)
    radiation_interval=options_get_integer(current_state%options_database, "rad_interval")
    socrates_enabled=is_component_enabled(current_state%options_database, "socrates_couple")
  end subroutine init_timestepper  

  !> Performs a timestep, which is comprised of executing each group of components in the order that they have been configured
  !! in. The components in a group can be called, depending on the type, just once per timestep (WHOLE) or per column (COLUMN).
  !! @param current_state The current model state
  subroutine timestep(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i

    call handle_sampling(current_state)

    do i=1,size(group_descriptors)      
      if (group_descriptors(i)%type == GROUP_TYPE_WHOLE) then
        call timestep_whole(current_state, group_descriptors(i))
      else if (group_descriptors(i)%type == GROUP_TYPE_COLUMN) then
        call timestep_column(current_state, group_descriptors(i))
      end if
    end do    
  end subroutine timestep

  !> Finalises the timestepper by cleaning up allocated memory
  subroutine finalise_timestepper()
    deallocate(group_descriptors)
  end subroutine finalise_timestepper  

  !> Performs timestepping for a group of components on a per column basis. Each component in the group is executed 
  !! for every column.
  !! @param current_state The current model state
  !! @param group_descriptor Description of the group of components to execute
  subroutine timestep_column(current_state, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    type(group_descriptor_type), intent(in) :: group_descriptor

    current_state%column_global_x=current_state%local_grid%start(X_INDEX) - current_state%local_grid%halo_size(X_INDEX)
    current_state%column_local_x=1
    do while (current_state%column_global_x .le. &
         current_state%local_grid%end(X_INDEX)+current_state%local_grid%halo_size(X_INDEX))
      current_state%column_global_y = current_state%local_grid%start(Y_INDEX) - current_state%local_grid%halo_size(Y_INDEX)
      current_state%column_local_y=1
      do while (current_state%column_global_y .le. &
           current_state%local_grid%end(Y_INDEX)+current_state%local_grid%halo_size(Y_INDEX))
        call update_state_sitation_flags(current_state)
        call execute_timestep_callbacks(current_state, group_descriptor%id)
        current_state%column_global_y = current_state%column_global_y + 1
        current_state%column_local_y = current_state%column_local_y + 1
      end do
      current_state%column_global_x = current_state%column_global_x + 1
      current_state%column_local_x = current_state%column_local_x + 1
    end do    
  end subroutine timestep_column

  !> Executes a timestep for components in a group which are designed to be executed once per timestep
  !! @param current_state The current model state
  !! @param group_descriptor Description of the group of components to execute
  subroutine timestep_whole(current_state, group_descriptor)
    type(model_state_type), intent(inout) :: current_state
    type(group_descriptor_type), intent(in) :: group_descriptor

    ! For print_debug_data, the column_global fields must match the requested coordinate.
    ! This is already handled for the timestep_column, but needs to be specially set for timestep_whole.
    !  This should not affect update_state_sitation_flags, as that is only used for timestep_column.
    if (current_state%print_debug_data) then
      current_state%column_global_x = current_state%pdd_x
      current_state%column_global_y = current_state%pdd_y
    end if
    call execute_timestep_callbacks(current_state, group_descriptor%id)
  end subroutine timestep_whole

  !> Updates the states situation flags for easy retrieval in the components that are
  !! run per timestep
  !! @param state The current model state
  subroutine update_state_sitation_flags(current_state)
    type(model_state_type), intent(inout) :: current_state

    current_state%first_timestep_column = (current_state%column_local_x == 1 .and. current_state%column_local_y == 1)
    current_state%last_timestep_column = (current_state%column_global_x == &
         current_state%local_grid%end(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) .and. &
         current_state%column_global_y == current_state%local_grid%end(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX))   

    current_state%first_nonhalo_timestep_column = (current_state%column_local_x == current_state%local_grid%halo_size(X_INDEX)+1 &
       .and. current_state%column_local_y == current_state%local_grid%halo_size(Y_INDEX)+1)

    current_state%halo_column = current_state%column_local_y .lt. current_state%local_grid%local_domain_start_index(Y_INDEX) .or.&
         current_state%column_local_x .lt. current_state%local_grid%local_domain_start_index(X_INDEX) .or.&
         current_state%column_local_y .gt. current_state%local_grid%local_domain_end_index(Y_INDEX) .or.&
         current_state%column_local_x .gt. current_state%local_grid%local_domain_end_index(X_INDEX)
  end subroutine update_state_sitation_flags

  !> Updates the diagnostic sampling flag for the new timestep
  !! @param state The current model state
  subroutine handle_sampling(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: i

    current_state%diagnostic_sample_timestep = .false.
    current_state%sampling(:)%active = .false.
    current_state%radiation_timestep = .false.  ! for computation timing under time_basis

    if (.not. current_state%only_compute_on_sample_timestep) then
      ! always compute the diagnostic in this case
      current_state%diagnostic_sample_timestep = .true.
    end if

    ! The following three cases will only compute diangnostics at requested intervals.
    ! However, it does ALL diagnostics regardless of specific request.
    ! MONC isn't STATSH-smart...though radiation diagnostics come close
    if (current_state%time_basis) then
      ! enable calculations and sampling at specified step only 
      ! (at sampling time interval, which is also an output or write interval)
      do i=1, size(current_state%sampling)
        if (current_state%timestep .eq. current_state%sampling(i)%next_step) then
          if (current_state%sampling(i)%radiation) then
            ! Only possible when socrates_enabled and radiation_interval .gt. 0 (iobridge)
            ! Permits radiation without needing to do all diagnostics
            ! Never set %active for the %radiation case - does not denote a iob data_definition
            current_state%radiation_timestep = .true. 
          else
            current_state%diagnostic_sample_timestep = .true.
            current_state%sampling(i)%active = .true.
          end if
        end if
      end do
    else ! timestep basis or force_output_on_interval
      ! enable radiation calculation
      if (socrates_enabled .and. radiation_interval .gt. 0) then
        if (mod(current_state%timestep, radiation_interval) == 0) &
             current_state%radiation_timestep = .true.
      end if 
      ! enable diagnostic calculation and sampling on the sampling timestep interval.
      do i=1,size(current_state%sampling)
        if (mod(current_state%timestep, current_state%sampling(i)%interval) == 0) then
          current_state%diagnostic_sample_timestep = .true.
          current_state%sampling(i)%active = .true.
        end if
      end do
    end if
  end subroutine handle_sampling

end module timestepper_mod
