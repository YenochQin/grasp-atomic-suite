"""Check production wavefunction procedure interfaces at alternate suite capacities."""

from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys
import tempfile


STATE_SOURCES = (
    "lib/grasp_shared/compile_time_capacities.f90",
    "lib/grasp_shared/grasp_parameters.f90",
    "lib/grasp_shared/radial_grid_state.f90",
    "lib/grasp_shared/scf_radial_potentials.f90",
    "lib/grasp_shared/radial_quadrature_workspace.f90",
)

# Compile the production module procedures and their USE dependencies in order.
# The probe only needs their .mod files for procedure-pointer association.
PROCEDURE_SOURCES = (
    "lib/grasp_shared/atomic_calculation_state.f90",
    "lib/grasp_shared/orbital_csf_basis.f90",
    "apps/rwfnestimate/integrate_radial_tail.f90",
    "lib/grasp_shared/radial_wavefunctions.f90",
    "lib/grasp_shared/standard_io_units.f90",
    "lib/grasp_shared/radial_dirac_workspace.f90",
    "lib/core/radial_dirac_coefficients.f90",
    "lib/grasp_shared/radial_difference_coefficients.f90",
    "lib/grasp_shared/nuclear_model_parameters_state.f90",
    "lib/grasp_shared/scf_iteration_state.f90",
    "lib/core/radial_origin_solution.f90",
    "lib/grasp_shared/radial_stepping_coefficients.f90",
    "apps/rwfnestimate/integrate_homogeneous_radial_segment.f90",
    "lib/grasp_shared/calculation_file_units.f90",
    "lib/grasp_shared/radial_node_count_settings.f90",
    "lib/core/radial_node_count.f90",
    "lib/grasp_shared/radial_quadrature_weights.f90",
    "lib/core/radial_quadrature.f90",
    "apps/rwfnestimate/solve_homogeneous_radial_equation.f90",
    "lib/grasp_shared/initial_orbital_selection.f90",
    "lib/grasp_shared/initial_orbital_source.f90",
    "apps/rwfnestimate/estimate_thomas_fermi_orbitals.f90",
    "lib/grasp_shared/debug_output_settings.f90",
    "lib/grasp_shared/array_memory_management.f90",
    "lib/core/file_open.f90",
    "lib/core/weighted_radial_overlap.f90",
    "lib/core/radial_orbital_interpolation.f90",
    "apps/rwfnestimate/estimate_from_saved_orbitals.f90",
    "apps/rwfnestimate/write_initial_wavefunction_summary.f90",
)


def main() -> None:
    compiler = sys.argv[1]
    suite = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="standard-grasp-capacity-") as temporary:
        work = Path(temporary)
        sources = [suite / "src" / source for source in STATE_SOURCES]
        capacities = work / sources[0].name
        text = sources[0].read_text()
        for name, original, replacement in (("max_radial_points", 590, 790), ("max_orbitals", 127, 149)):
            pattern = rf"^(\s*integer\s*,\s*parameter\s*::\s*{name}\s*=\s*){original}\b"
            text, count = re.subn(pattern, rf"\g<1>{replacement}", text, flags=re.IGNORECASE | re.MULTILINE)
            assert count == 1, f"alternate-capacity test expects one pristine {name} declaration"
        capacities.write_text(text)
        sources[0] = capacities

        flags = [compiler, "-std=f2018", "-fcheck=all"]
        result = subprocess.run(
            [*flags, "-fsyntax-only", *sources, *(suite / "src" / source for source in PROCEDURE_SOURCES)],
            cwd=work, capture_output=True, text=True, timeout=60, check=False,
        )
        assert result.returncode == 0, result.stderr
        for fixture in ["shared_grid_capacity.f90", "standard_grasp/interface_capacity.f90"]:
            binary = work / Path(fixture).stem
            result = subprocess.run(
                [*flags, *sources, suite / "test" / fixture, "-o", binary],
                cwd=work, capture_output=True, text=True, timeout=60, check=False,
            )
            assert result.returncode == 0, result.stderr
            subprocess.run([str(binary)], cwd=work, check=True, timeout=10)
        print("Ordinary GRASP radial/orbital interfaces: 790 / 800 / 149")


if __name__ == "__main__":
    main()
