! Re-exports the shared capacities and their established packing and dimension aliases.
      module grasp_parameters

      use compile_time_capacities, only: orbital_index_radix
      use compile_time_capacities, only: max_radial_points
      use compile_time_capacities, only: radial_workspace_capacity
      use compile_time_capacities, only: max_orbitals
      use compile_time_capacities, only: max_orbitals_minus_one
      use compile_time_capacities, only: max_orbitals_minus_two
      ! The legacy names refer to the same constants.
      use compile_time_capacities, only: keyorb => orbital_index_radix
      use compile_time_capacities, only: nnnp => max_radial_points
      use compile_time_capacities, only: nnn1 => radial_workspace_capacity
      use compile_time_capacities, only: nnnw => max_orbitals
      use compile_time_capacities, only: nnnwm1 => max_orbitals_minus_one
      use compile_time_capacities, only: nnnwm2 => max_orbitals_minus_two
      implicit none (type, external)
      end module grasp_parameters
