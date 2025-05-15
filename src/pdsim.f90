program simulate_inception
  use m_config
  use m_pdsim
  use m_cross_sec
  use m_particle_core
  use m_units_constants
  use m_gas
  use m_photoi
  use m_lookup_table
  use m_interp_unstructured

  implicit none

  type(CFG_t)   :: cfg
  type(pdsim_t) :: pd
  type(PC_part_t) :: initial_electron

  integer              :: iteration
  real(dp)             :: dt, time
  logical              :: output_particles_now
  integer              :: max_iterations = huge(1)
  integer              :: n_output, i_start
  integer, allocatable :: i_start_cells(:)
  integer, allocatable :: n_ionizations(:)

  call pdsim_create_config(cfg)

  call CFG_update_from_arguments(cfg)
  call CFG_check(cfg)

  call pdsim_initialize(cfg, pd)

  call pdsim_get_gas_cells(pdsim_ug, pd%n_initial_positions, &
       pd%initial_rmin, pd%initial_rmax, i_start_cells)
  allocate(n_ionizations(pd%n_initial_positions))
  n_ionizations(:) = 0

  ! Set properties for initial electron
  initial_electron%v(:)   = 0.0_dp
  initial_electron%w      = 1.0_dp
  initial_electron%t_left = 0.0_dp
  initial_electron%a(:)   = 0.0_dp

  do i_start = 1, pd%n_initial_positions
     call pd%pc%remove_particles()

     ! Create electron at center of cell
     initial_electron%x = iu_get_cell_center(pdsim_ug, i_start_cells(i_start))
     call pd%pc%add_part(initial_electron)

     time     = 0.0_dp
     dt       = pd%max_dt
     n_output = 0

     do iteration = 1, max_iterations
        ! Set acceleration
        call pd%pc%set_accel()

        output_particles_now = (pd%output_particles .and. &
             time >= n_output * pd%output_particles_dt)

        if (output_particles_now) then
           n_output = n_output + 1
           call pdsim_write_particles(pd, n_output, i_start)
        end if

        call pd%pc%advance(dt)
        call pd%pc%after_mover(dt)
        call pd%pc%clean_up()
        time = time + dt

        call pdsim_handle_events(pd%pc, n_ionizations(i_start))

        if (pd%pc%n_part == 0) exit
        if (pd%pc%n_part >= pd%num_electrons_inception) exit
     end do

     print *, i_start, iteration, n_ionizations(i_start)
  end do

  call write_ionization_count(i_start_cells, n_ionizations)

contains

  subroutine write_ionization_count(i_start_cells, n_ionizations)
    integer, intent(in) :: i_start_cells(:)
    integer, intent(in) :: n_ionizations(:)

    character(len=200)        :: fname
    integer                   :: my_unit, n

    write(fname, "(A,I0.5,A,I0.5,A)") trim(pd%output_name) // &
         "_ionization_count.3D"

    open(newunit=my_unit, file=trim(fname), action="write")
    write(my_unit, *) "X Y Z n_ionizations"
    do n = 1, size(i_start_cells)
       write(my_unit, *) iu_get_cell_center(pdsim_ug, i_start_cells(n)), &
            n_ionizations(n)
    end do
    close(my_unit)

    print *, "Wrote ", trim(fname)
  end subroutine write_ionization_count

end program simulate_inception
