!> Add random noise into the fields
module randomnoise_mod
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : FORWARD_STEPPING, model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
       options_get_logical_array, options_get_real_array, options_get_string_array, options_get_array_size, &
       options_compare_profile_arrays
  use interpolation_mod, only: piecewise_linear_1d
  use q_indices_mod, only: get_q_index, standard_q_names
  use logging_mod, only: LOG_ERROR, log_master_log
  use conversions_mod, only : conv_to_string

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer, parameter :: MAX_SIZE_SEED_ARRAY=256,  &  ! large value to work on multiple systems
                        I_SEED     = 7,           &  ! initial seed value,        non-reproducible case
                        THETA_SEED = -1731191804, &  ! initial seed for theta,    reproducible case
                        Q_SEED     =  2011234875, &  ! initial seed for q-fields, reproducible case
                        W_SEED     =  -163411914     ! initial seed for w,        reproducible case

  public randomnoise_get_descriptor
contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The termination check component descriptor
  type(component_descriptor_type) function randomnoise_get_descriptor()
    randomnoise_get_descriptor%name="randomnoise"
    randomnoise_get_descriptor%version=0.1
    randomnoise_get_descriptor%initialisation=>initialisation_callback
  end function randomnoise_get_descriptor

  !> The initialisation callback sets up the buoyancy coefficient
  !! @param current_state The current model state
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer, dimension(MAX_SIZE_SEED_ARRAY) :: iranseed
    real(kind=DEFAULT_PRECISION) :: random_num

    integer :: nq_rand   ! The number of q fields to add noise to
    integer :: nzq       ! The number of input levels for noise
    integer :: i,j,k,n,s ! loop counters
    integer :: iq        ! temporary q varible index
    integer :: iloc, jloc

    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: f_rand_pl_q   ! Random Noise node amplitude for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_rand_pl_q     ! Random Noise node height values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_theta ! Random Noise node amplitude for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_rand_pl_theta ! Random Noise node height values for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_rand_pl_w     ! Random Noise node amplitude for theta variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_rand_pl_w     ! Random Noise node height values for theta variable

    logical :: l_rand_pl_theta ! if .true. then random noise added to theta field
    logical :: l_rand_pl_q     ! if .true. then random noise added to q fields
    logical :: l_rand_pl_w     ! if .true. then random noise added to w field
    logical :: l_rand_bit_reproducible ! if .true. then is bit reproducible between runs

    character(len=STRING_LENGTH), dimension(:), allocatable :: names_rand_pl_q ! names of q variables to add random noise to

    real(kind=DEFAULT_PRECISION), allocatable :: f_rand_pl_q_tmp(:) !temporary 1D storage of random noise for q field
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid to use in interpolation

    if (current_state%continuation_run) return

    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))

    l_rand_pl_theta=options_get_logical(current_state%options_database, "l_rand_pl_theta")
    l_rand_pl_w=options_get_logical(current_state%options_database, "l_rand_pl_w")
    l_rand_pl_q=options_get_logical(current_state%options_database, "l_rand_pl_q")
    l_rand_bit_reproducible=options_get_logical(current_state%options_database, "l_rand_bit_reproducible")

    ! Initialise random seed to be different in every MONC
    !  This is overwritten if l_rand_bit_reproducible=.true.
    iranseed = I_SEED + current_state%parallel%my_rank
    call random_seed(put = iranseed)

    if (l_rand_pl_q) then

      allocate(names_rand_pl_q(options_get_array_size(current_state%options_database, "names_rand_pl_q")))
      if (size(names_rand_pl_q) .eq. 0) then
        call log_master_log(LOG_ERROR, "Model configured with l_rand_pl_q=.true., but no names_rand_pl_q "//&
                                       "have been provided to specify the species names to perturb.")
      end if

      call options_get_string_array(current_state%options_database, "names_rand_pl_q", names_rand_pl_q)
    end if


    !---------------------------------------------------
    ! Handle randomnoise for theta
    if (l_rand_pl_theta)then

      ! Get amplitude profiles
      allocate(z_rand_pl_theta(options_get_array_size(current_state%options_database, "z_rand_pl_theta")), &
           f_rand_pl_theta(options_get_array_size(current_state%options_database, "f_rand_pl_theta")))
      call options_compare_profile_arrays(current_state%options_database, &
                               "z_rand_pl_theta", "f_rand_pl_theta", "theta perturbation")
      call options_get_real_array(current_state%options_database, "z_rand_pl_theta", z_rand_pl_theta)
      call options_get_real_array(current_state%options_database, "f_rand_pl_theta", f_rand_pl_theta)
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_rand_pl_theta(1:size(z_rand_pl_theta)), f_rand_pl_theta(1:size(f_rand_pl_theta)), zgrid, &
           current_state%global_grid%configuration%vertical%theta_rand)

      ! Cycle through all points in the global grid, record if on the local grid
      if (l_rand_bit_reproducible) then
        iranseed(:) = THETA_SEED
        call random_seed(put = iranseed)
        do i=1, current_state%global_grid%size(X_INDEX)
          iloc = i - current_state%local_grid%start(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) + 1
          do j=1, current_state%global_grid%size(Y_INDEX)
            jloc = j - current_state%local_grid%start(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) + 1
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
              call random_number(random_num)

              if (iloc .ge. current_state%local_grid%local_domain_start_index(X_INDEX) .and. &
                  iloc .le. current_state%local_grid%local_domain_end_index(X_INDEX)   .and. &
                  jloc .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) .and. &
                  jloc .le. current_state%local_grid%local_domain_end_index(Y_INDEX) ) then
                  current_state%th%data(k,jloc,iloc) = current_state%th%data(k,jloc,iloc) + &
                     current_state%global_grid%configuration%vertical%theta_rand(k) * 2.0 * (random_num-0.5)

              end if ! in local space?
            end do ! k
          end do ! j
        end do ! i 

      else ! not l_rand_bit_reproducible
            
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
  
              ! Apply random number for this grid point
              !   if not l_rand_bit_reproducible, then random_num is determined by local rank
              call random_number(random_num)
              current_state%th%data(k,j,i) = current_state%th%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%theta_rand(k) * 2.0 * (random_num-0.5)
             end do
          end do
        end do
      end if !check l_rand_bit_reproducible
      deallocate(z_rand_pl_theta, f_rand_pl_theta)
    end if ! l_rand_pl_theta



    !---------------------------------------------------
    ! Handle randomnoise for q-fields
    if (l_rand_pl_q)then
      nq_rand=size(names_rand_pl_q)
      allocate(z_rand_pl_q(options_get_array_size(current_state%options_database, "z_rand_pl_q")))
      call options_get_real_array(current_state%options_database, "z_rand_pl_q", z_rand_pl_q)
      nzq=size(z_rand_pl_q)
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      allocate(f_rand_pl_q_tmp(options_get_array_size(current_state%options_database, "f_rand_pl_q")))
      if (nq_rand*nzq .ne. size(f_rand_pl_q_tmp)) then
        call log_master_log(LOG_ERROR, "There is a mismatch between the number of moisture perturbation heights, "//&
                                       "size(z_rand_pl_q)="//trim(conv_to_string(nzq))//                            &
                                       ", and the perturbation values, "//                                               &
                                       "size(f_rand_pl_q)="//trim(conv_to_string(size(f_rand_pl_q_tmp)))//          &
                                       ".  The length of f_rand_pl_q should equal the length of z_rand_pl_q "//     &
                                       "multiplied by the number of names_rand_pl_q.")
      end if
      call options_get_real_array(current_state%options_database, "f_rand_pl_q", f_rand_pl_q_tmp)
      allocate(f_rand_pl_q(nzq, nq_rand))
      f_rand_pl_q(1:nzq, 1:nq_rand)=reshape(f_rand_pl_q_tmp, (/nzq, nq_rand/))
      do n=1,nq_rand

        iq=get_q_index(trim(names_rand_pl_q(n)), 'random noise')
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        call piecewise_linear_1d(z_rand_pl_q(1:size(z_rand_pl_q)), f_rand_pl_q(1:nzq,n), zgrid, &
             current_state%global_grid%configuration%vertical%q_rand(:,iq))

        ! Cycle through all points in the global grid, record if on the local grid
        if (l_rand_bit_reproducible) then
          iranseed(:) = Q_SEED/n
          call random_seed(put = iranseed)
          do i=1, current_state%global_grid%size(X_INDEX)
            iloc = i - current_state%local_grid%start(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) + 1
            do j=1, current_state%global_grid%size(Y_INDEX)
              jloc = j - current_state%local_grid%start(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) + 1
              do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
                call random_number(random_num)
  
                if (iloc .ge. current_state%local_grid%local_domain_start_index(X_INDEX) .and. &
                    iloc .le. current_state%local_grid%local_domain_end_index(X_INDEX)   .and. &
                    jloc .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) .and. &
                    jloc .le. current_state%local_grid%local_domain_end_index(Y_INDEX) ) then
  
                  current_state%q(iq)%data(k,jloc,iloc) = current_state%q(iq)%data(k,jloc,iloc) + &
                       current_state%global_grid%configuration%vertical%q_rand(k,iq) * 2.0 * (random_num-0.5)
                end if ! in local space?
              end do ! k
            end do ! j
          end do ! i
  
        else ! not l_rand_bit_reproducible
  
          do i=current_state%local_grid%local_domain_start_index(X_INDEX), &
               current_state%local_grid%local_domain_end_index(X_INDEX)
            do j=current_state%local_grid%local_domain_start_index(Y_INDEX), &
                 current_state%local_grid%local_domain_end_index(Y_INDEX)
              do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
  
                ! Apply random number for this grid point
                !   if not l_rand_bit_reproducible, then random_num is determined by local rank
                call random_number(random_num)
                current_state%q(iq)%data(k,j,i) = current_state%q(iq)%data(k,j,i) + &
                       current_state%global_grid%configuration%vertical%q_rand(k,iq) * 2.0 * (random_num-0.5)
              end do ! k
            end do ! j
          end do ! i
        end if ! l_rand_bit_reproducible
      end do ! n
      deallocate(z_rand_pl_q, f_rand_pl_q_tmp, f_rand_pl_q, names_rand_pl_q)
    end if ! l_rand_pl_q


    !---------------------------------------------------
    ! Handle randomnoise for w
    if (l_rand_pl_w)then

      ! Get amplitude profiles
      allocate(z_rand_pl_w(options_get_array_size(current_state%options_database, "z_rand_pl_w")), &
           f_rand_pl_w(options_get_array_size(current_state%options_database, "f_rand_pl_w")))
      call options_compare_profile_arrays(current_state%options_database, &
                                "z_rand_pl_w", "f_rand_pl_w", "w perturbation")
      call options_get_real_array(current_state%options_database, "z_rand_pl_w", z_rand_pl_w)
      call options_get_real_array(current_state%options_database, "f_rand_pl_w", f_rand_pl_w)

      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_rand_pl_w(1:size(z_rand_pl_w)), f_rand_pl_w(1:size(f_rand_pl_w)), zgrid, &
           current_state%global_grid%configuration%vertical%w_rand)

      ! Cycle through all points in the global grid, record if on the local grid
      if (l_rand_bit_reproducible) then
        iranseed(:) = W_SEED
        call random_seed(put = iranseed)
        do i=1, current_state%global_grid%size(X_INDEX)
          iloc = i - current_state%local_grid%start(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) + 1
          do j=1, current_state%global_grid%size(Y_INDEX)
            jloc = j - current_state%local_grid%start(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) + 1
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
              call random_number(random_num)

              if (iloc .ge. current_state%local_grid%local_domain_start_index(X_INDEX) .and. &
                  iloc .le. current_state%local_grid%local_domain_end_index(X_INDEX)   .and. &
                  jloc .ge. current_state%local_grid%local_domain_start_index(Y_INDEX) .and. &
                  jloc .le. current_state%local_grid%local_domain_end_index(Y_INDEX) ) then
                current_state%w%data(k,jloc,iloc) = current_state%w%data(k,jloc,iloc) + &
                   current_state%global_grid%configuration%vertical%w_rand(k) * (random_num-0.5)

              end if ! in local space?
            end do ! k
          end do ! j
        end do ! i

      else ! not l_rand_bit_reproducible

        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)

              ! Apply random number for this grid point
              !   if not l_rand_bit_reproducible, then random_num is determined by local rank
              call random_number(random_num)
              current_state%w%data(k,j,i) = current_state%w%data(k,j,i) + &
                     current_state%global_grid%configuration%vertical%w_rand(k) * (random_num-0.5)
            end do
#ifdef W_ACTIVE
            current_state%w%data(current_state%local_grid%local_domain_end_index(Z_INDEX),j,i)=0.0_DEFAULT_PRECISION
            current_state%w%data(1,j,i)=0.0_DEFAULT_PRECISION
#endif
            if (current_state%use_viscosity_and_diffusion) then
#ifdef U_ACTIVE
              current_state%u%data(1,j,i)=-current_state%u%data(2,j,i)
#endif
#ifdef V_ACTIVE
              current_state%v%data(1,j,i)=-current_state%v%data(2,j,i)
#endif
            else
#ifdef U_ACTIVE
              current_state%u%data(1,j,i)=current_state%u%data(2,j,i)
#endif
#ifdef V_ACTIVE
              current_state%v%data(1,j,i)=current_state%v%data(2,j,i)
#endif
            end if
          end do ! j
        end do ! i
      end if ! l_rand_bit_reproducible
      deallocate(z_rand_pl_w, f_rand_pl_w)
    end if ! l_rand_pl_w
    deallocate(zgrid)
  end subroutine initialisation_callback
end module randomnoise_mod
