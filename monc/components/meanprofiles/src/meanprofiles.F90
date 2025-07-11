!> Calculates the mean profiles of prognostic variables which are then used in smoothing and other areas
module meanprofiles_mod
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use mpi, only : MPI_SUM, MPI_IN_PLACE
  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: start_x, end_x, start_y, end_y, bar_fields

  real(kind=DEFAULT_PRECISION) :: rnhpts
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: bartmp
  ! Immersed boundary variables
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: rnspts_u, rnspts_v, rnspts_s

  public meanprofiles_get_descriptor
contains

  !> Returns the component descriptor of the mean profiles module
  !! @returns Component descriptor of mean profiles
  type(component_descriptor_type) function meanprofiles_get_descriptor()
    meanprofiles_get_descriptor%name="mean_profiles"
    meanprofiles_get_descriptor%version=0.1
    meanprofiles_get_descriptor%initialisation=>init_callback
    meanprofiles_get_descriptor%timestep=>timestep_callback
    meanprofiles_get_descriptor%finalisation=>finalisation_callback
  end function meanprofiles_get_descriptor

  !> Called on MONC initialisation, will allocate appropriate data structures
  !! @param current_state The current model state
  subroutine init_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION) :: tmp
    integer :: k, ierr, myproc

    myproc = current_state%parallel%my_rank
    

    bar_fields=0

    rnhpts=1.0_DEFAULT_PRECISION/real(current_state%global_grid%size(X_INDEX)*current_state%global_grid%size(Y_INDEX))

    start_x=current_state%local_grid%local_domain_start_index(X_INDEX)
    end_x=current_state%local_grid%local_domain_end_index(X_INDEX)
    start_y=current_state%local_grid%local_domain_start_index(Y_INDEX)
    end_y=current_state%local_grid%local_domain_end_index(Y_INDEX)

#ifdef U_ACTIVE
    if (.not. current_state%continuation_run) then
      allocate(current_state%global_grid%configuration%vertical%olubar(current_state%local_grid%size(Z_INDEX)),&
           current_state%global_grid%configuration%vertical%olzubar(current_state%local_grid%size(Z_INDEX)))
    end if
    bar_fields=bar_fields+2
#endif
#ifdef V_ACTIVE
    if (.not. current_state%continuation_run) then
      allocate(current_state%global_grid%configuration%vertical%olvbar(current_state%local_grid%size(Z_INDEX)),&
           current_state%global_grid%configuration%vertical%olzvbar(current_state%local_grid%size(Z_INDEX)))
    end if
    bar_fields=bar_fields+2
#endif
    if (current_state%th%active) then
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olthbar(current_state%local_grid%size(Z_INDEX)),&
             current_state%global_grid%configuration%vertical%olzthbar(current_state%local_grid%size(Z_INDEX)))
      end if
      bar_fields=bar_fields+2
    end if
    if (current_state%number_q_fields .gt. 0) then
      bar_fields=bar_fields+(current_state%number_q_fields*2)
      if (.not. current_state%continuation_run) then
        allocate(current_state%global_grid%configuration%vertical%olqbar(current_state%local_grid%size(Z_INDEX), &
             current_state%number_q_fields), current_state%global_grid%configuration%vertical%olzqbar(&
             current_state%local_grid%size(Z_INDEX), current_state%number_q_fields))
      end if
    end if
    allocate(bartmp(current_state%local_grid%size(Z_INDEX), bar_fields))


    ! Immersed boundary
    if(current_state%immersed%ib_enabled) then

      tmp = real(current_state%local_grid%size(X_INDEX)*current_state%local_grid%size(Y_INDEX))

#ifdef U_ACTIVE
      allocate(rnspts_u(current_state%local_grid%size(Z_INDEX)))
      rnspts_u = 1.0_DEFAULT_PRECISION
      do k=1,current_state%local_grid%size(Z_INDEX)
        rnspts_u(k) = tmp-sum(current_state%immersed%indic_u(k, start_y:end_y, start_x:end_x))
      end do
      call mpi_allreduce(MPI_IN_PLACE, rnspts_u, current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
                         current_state%parallel%monc_communicator, ierr)
      do k=1,current_state%local_grid%size(Z_INDEX)
        if(rnspts_u(k).gt.0.0)rnspts_u(k) = 1.0_DEFAULT_PRECISION/rnspts_u(k)
      end do
#endif

#ifdef V_ACTIVE
      allocate(rnspts_v(current_state%local_grid%size(Z_INDEX)))
      rnspts_v = 1.0_DEFAULT_PRECISION
      do k=1,current_state%local_grid%size(Z_INDEX)
        rnspts_v(k) = tmp-sum(current_state%immersed%indic_v(k, start_y:end_y, start_x:end_x))
      end do
      call mpi_allreduce(MPI_IN_PLACE, rnspts_v, current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
                         current_state%parallel%monc_communicator, ierr)
      do k=1,current_state%local_grid%size(Z_INDEX)
        if(rnspts_v(k).gt.0.0)rnspts_v(k) = 1.0_DEFAULT_PRECISION/rnspts_v(k)
      end do
#endif

      if(current_state%th%active.or.(current_state%number_q_fields.gt.0)) then
        allocate(rnspts_s(current_state%local_grid%size(Z_INDEX)))
        rnspts_s = 1.0_DEFAULT_PRECISION
        do k=1,current_state%local_grid%size(Z_INDEX)
          rnspts_s(k) = tmp-sum(current_state%immersed%indic_s(k, start_y:end_y, start_x:end_x))
        end do
        call mpi_allreduce(MPI_IN_PLACE, rnspts_s, current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
                           current_state%parallel%monc_communicator, ierr)
        do k=1,current_state%local_grid%size(Z_INDEX)
          if(rnspts_s(k).gt.0.0)rnspts_s(k) = 1.0_DEFAULT_PRECISION/rnspts_s(k)
        end do
      end if

    end if


    ! Do the initial calculation for the first timestep
    if (.not. current_state%continuation_run) call calculate_mean_profiles(current_state)
  end subroutine init_callback

  !> Will recalculate the mean profiles of each prognostic when called (for the entire local domain)
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call calculate_mean_profiles(current_state)
    
  end subroutine timestep_callback

  !> Frees up the temporary data for the bars
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (allocated(bartmp)) deallocate(bartmp)    
  end subroutine finalisation_callback

  !> Calculates the global mean profiles and stores these in the ol bar arrays
  !! @param current_state The current model state
  subroutine calculate_mean_profiles(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: bar_index, i, k
    
    call calculate_sum_profiles(current_state)
    
    if(current_state%immersed%ib_enabled)then

      bar_index=1
#ifdef U_ACTIVE
      current_state%global_grid%configuration%vertical%olubar(:)=bartmp(:, bar_index)*rnspts_u(:)
      current_state%global_grid%configuration%vertical%olzubar(:)=bartmp(:, bar_index+1)*rnspts_u(:)
      bar_index=bar_index+2
#endif
#ifdef V_ACTIVE
      current_state%global_grid%configuration%vertical%olvbar(:)=bartmp(:, bar_index)*rnspts_v(:)
      current_state%global_grid%configuration%vertical%olzvbar(:)=bartmp(:, bar_index+1)*rnspts_v(:)
      bar_index=bar_index+2
#endif
      if (current_state%th%active) then
        current_state%global_grid%configuration%vertical%olthbar(:)=bartmp(:, bar_index)*rnspts_s(:)
        current_state%global_grid%configuration%vertical%olzthbar(:)=bartmp(:, bar_index+1)*rnspts_s(:)
        bar_index=bar_index+2
      end if
      do i=1,current_state%number_q_fields
        if (current_state%q(i)%active) then
          current_state%global_grid%configuration%vertical%olqbar(:, i)=bartmp(:, bar_index)*rnspts_s(:)
          current_state%global_grid%configuration%vertical%olzqbar(:, i)=bartmp(:, bar_index+1)*rnspts_s(:)
          bar_index=bar_index+2
        end if
      end do

    else

      bar_index=1
#ifdef U_ACTIVE
      current_state%global_grid%configuration%vertical%olubar(:)=bartmp(:, bar_index)*rnhpts
      current_state%global_grid%configuration%vertical%olzubar(:)=bartmp(:, bar_index+1)*rnhpts
      bar_index=bar_index+2
#endif
#ifdef V_ACTIVE
      current_state%global_grid%configuration%vertical%olvbar(:)=bartmp(:, bar_index)*rnhpts
      current_state%global_grid%configuration%vertical%olzvbar(:)=bartmp(:, bar_index+1)*rnhpts
      bar_index=bar_index+2
#endif
      if (current_state%th%active) then
        current_state%global_grid%configuration%vertical%olthbar(:)=bartmp(:, bar_index)*rnhpts
        current_state%global_grid%configuration%vertical%olzthbar(:)=bartmp(:, bar_index+1)*rnhpts
        bar_index=bar_index+2
      end if
      do i=1,current_state%number_q_fields
        if (current_state%q(i)%active) then
          current_state%global_grid%configuration%vertical%olqbar(:, i)=bartmp(:, bar_index)*rnhpts
          current_state%global_grid%configuration%vertical%olzqbar(:, i)=bartmp(:, bar_index+1)*rnhpts
          bar_index=bar_index+2
        end if
      end do


    end if

  end subroutine calculate_mean_profiles    

  !> Calculates the sum profiles for the bars for each level globally
  !! @param current_state The current model state_mod
  subroutine calculate_sum_profiles(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: k, n, bar_index, ierr    

    if(current_state%immersed%ib_enabled)then
      do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
        bar_index=1
#ifdef U_ACTIVE
        bartmp(k, bar_index)=sum(current_state%u%data(k, start_y:end_y, start_x:end_x),&
                                 mask=current_state%immersed%indic_u(k, start_y:end_y,start_x:end_x).eq.0 )
        bartmp(k, bar_index+1)=sum(current_state%zu%data(k, start_y:end_y, start_x:end_x),&
                                   mask=current_state%immersed%indic_u(k, start_y:end_y,start_x:end_x).eq.0 )
        bar_index=bar_index+2
#endif
#ifdef V_ACTIVE
        bartmp(k, bar_index)=sum(current_state%v%data(k, start_y:end_y, start_x:end_x),&
                                 mask=current_state%immersed%indic_v(k, start_y:end_y,start_x:end_x).eq.0 )
        bartmp(k, bar_index+1)=sum(current_state%zv%data(k, start_y:end_y, start_x:end_x),&      
                                 mask=current_state%immersed%indic_v(k, start_y:end_y,start_x:end_x).eq.0 )
        bar_index=bar_index+2
#endif
        if (current_state%th%active) then
          bartmp(k, bar_index)=sum(current_state%th%data(k, start_y:end_y, start_x:end_x),&
                                   mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
          bartmp(k, bar_index+1)=sum(current_state%zth%data(k, start_y:end_y, start_x:end_x),&
                                     mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
          bar_index=bar_index+2
        end if
        do n=1,current_state%number_q_fields
          if (current_state%q(n)%active) then
            bartmp(k, bar_index)=sum(current_state%q(n)%data(k, start_y:end_y, start_x:end_x),&
                                     mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
            bartmp(k, bar_index+1)=sum(current_state%zq(n)%data(k, start_y:end_y, start_x:end_x),&
                                       mask=current_state%immersed%indic_s(k, start_y:end_y,start_x:end_x).eq.0 )
            bar_index=bar_index+2
          end if
        end do
      end do

    else

     do k=current_state%local_grid%local_domain_start_index(Z_INDEX), current_state%local_grid%local_domain_end_index(Z_INDEX)
        bar_index=1
#ifdef U_ACTIVE
        bartmp(k, bar_index)=sum(current_state%u%data(k, start_y:end_y, start_x:end_x))
        bartmp(k, bar_index+1)=sum(current_state%zu%data(k, start_y:end_y, start_x:end_x))
        bar_index=bar_index+2
#endif
#ifdef V_ACTIVE
        bartmp(k, bar_index)=sum(current_state%v%data(k, start_y:end_y, start_x:end_x))
        bartmp(k, bar_index+1)=sum(current_state%zv%data(k, start_y:end_y, start_x:end_x))      
        bar_index=bar_index+2
#endif
        if (current_state%th%active) then
          bartmp(k, bar_index)=sum(current_state%th%data(k, start_y:end_y, start_x:end_x))
          bartmp(k, bar_index+1)=sum(current_state%zth%data(k, start_y:end_y, start_x:end_x))
          bar_index=bar_index+2
        end if
        do n=1,current_state%number_q_fields
          if (current_state%q(n)%active) then
            bartmp(k, bar_index)=sum(current_state%q(n)%data(k, start_y:end_y, start_x:end_x))
            bartmp(k, bar_index+1)=sum(current_state%zq(n)%data(k, start_y:end_y, start_x:end_x))
            bar_index=bar_index+2
          end if
        end do
      end do

    end if

    call mpi_allreduce(MPI_IN_PLACE, bartmp, bar_fields*current_state%local_grid%size(Z_INDEX), PRECISION_TYPE, MPI_SUM, &
         current_state%parallel%monc_communicator, ierr)
  end subroutine calculate_sum_profiles
end module meanprofiles_mod
