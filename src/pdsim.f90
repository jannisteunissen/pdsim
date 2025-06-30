program pdsim
  use m_config
  use m_interp_unstructured

  use m_pdsim
  use m_photoi
  use m_integral
  use m_avalanche
  use m_particles

  implicit none

  type(CFG_t) :: cfg

  character(len=40) :: what_to_simulate

  call CFG_add(cfg, "simulate", undefined_str, &
         "What to simulate (integral, avalanches, particles)", &
         required=.true.)

  call pdsim_create_config(cfg)
  call photoi_create_cfg(cfg)
  call integral_create_config(cfg)
  call avalanche_create_config(cfg)
  call particles_create_config(cfg)

  call CFG_update_from_arguments(cfg)
  call CFG_check(cfg)

  call pdsim_initialize(cfg)
  call photoi_initialize(cfg)
  call CFG_get(cfg, "simulate", what_to_simulate)

  select case (what_to_simulate)
     case ("integral")
        print *, "Computing ionization integrals"
        call integral_compute(cfg)
     case ("avalanche")
        print *, "Computing ionization integrals"
        call integral_compute(cfg)
        print *, "Starting avalanche simulation"
        call avalanche_simulate(cfg)
     case ("particles")
        print *, "Starting particle simulation"
        call particles_simulate(cfg)
     case default
        error stop "simulate option can be integral, avalanche, particles"
  end select

  call iu_write_vtk(pdsim_ug, trim(pdsim_output_name) // ".vtu")
  print *, "Wrote ", trim(pdsim_output_name) // ".vtu"

end program pdsim
