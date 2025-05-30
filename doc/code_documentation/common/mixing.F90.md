## Overview

The `esl_mixing_m` module provides functionalities for mixing quantities, typically electron densities, during Self-Consistent Field (SCF) calculations. It defines a derived type `mixing_t` that stores the mixing parameter and offers a `linear` mixing procedure. Linear mixing is a common technique used to improve the stability and convergence of SCF iterations by combining the density (or another quantity) from the previous iteration with the newly computed density.

## Key Components

- **Module:** `esl_mixing_m`
    - **Description:** Implements linear mixing, primarily for use in SCF charge density mixing.
- **Type:** `mixing_t`
    - **Description:** A derived type designed to hold parameters and methods related to a mixing scheme.
    - **Components:**
        - `alpha (real(dp))`: The linear mixing parameter. This value determines the weight of the newly calculated quantity in the mixed output. A typical range is (0, 1].
    - **Procedures (public, type-bound):**
        - `init`: Initializes the `mixing_t` object, primarily by reading the `alpha` parameter from an input file.
        - `linear`: Performs the element-wise linear mixing of two arrays.

### Type-Bound Procedures

#### `init(this)`
- **Description:** Initializes the `mixing_t` object. Its main task is to set the mixing parameter `this%alpha`. The value of `alpha` is read from an FDF (Flexible Data Format) input file using the key `SCF.Mix.alpha`. If this key is not found in the input file, `alpha` defaults to `0.1_dp`.
- **Arguments:**
    - `this (class(mixing_t), intent(inout))`: The `mixing_t` object to be initialized.
- **Returns:** Not applicable (Subroutine).

#### `linear(this, np, in, out, next)`
- **Description:** Performs an element-wise linear mixing of two input arrays, `in` and `out`, and stores the result in the `next` array. The mixing formula applied to each element `ip` from 1 to `np` is:
  `next(ip) = in(ip) * beta + out(ip) * this%alpha`
  where `beta = 1.0_dp - this%alpha`.
  This means the `next` quantity is a weighted average of the `in` (e.g., old density) and `out` (e.g., new density) quantities.
- **Arguments:**
    - `this (class(mixing_t), intent(in))`: The `mixing_t` object containing the mixing parameter `alpha`.
    - `np (integer, intent(in))`: The number of elements in the arrays `in`, `out`, and `next`.
    - `in (real(dp), intent(in) :: in(:))`: The array representing the quantity from the previous iteration (e.g., input density $\rho_{in}$).
    - `out (real(dp), intent(in) :: out(:))`: The array representing the quantity computed in the current iteration (e.g., output density $\rho_{out}$).
    - `next (real(dp), intent(inout) :: next(:))`: The output array where the mixed quantity (e.g., density for the next iteration $\rho_{next}$) is stored. This array can be the same as the `in` array if `in` is also `intent(inout)`.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants
- **`alpha (real(dp))`**: The primary parameter within the `mixing_t` type, controlling the degree of mixing between the input and output quantities.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within an SCF loop:
!
! module scf_cycle_m
!   use esl_mixing_m, only: mixing_t
!   use prec, only: dp
!   implicit none
!
!   subroutine perform_scf(mixer_settings, initial_density, num_points)
!     type(mixing_t), intent(inout) :: mixer_settings ! Initialized elsewhere
!     real(dp), dimension(:), allocatable, intent(inout) :: initial_density ! rho_in, becomes rho_next
!     integer, intent(in) :: num_points
!
!     real(dp), dimension(size(initial_density)) :: current_output_density ! rho_out
!     integer :: iter
!
!     ! Call mixer_settings%init() would have been done before this routine
!
!     do iter = 1, max_iterations
!       ! 1. Compute current_output_density based on initial_density (representing rho_in)
!       !    (e.g., solve KS equations, then compute density)
!       call calculate_new_density(initial_density, current_output_density)
!
!       ! 2. Check for convergence (not shown)
!       !    if converged, exit loop
!
!       ! 3. Mix to get density for the next iteration
!       !    initial_density (rho_in) and current_output_density (rho_out) are mixed
!       !    into initial_density (which becomes rho_next for the subsequent iteration)
!       call mixer_settings%linear(num_points, initial_density, current_output_density, initial_density)
!     end do
!
!   contains
!     subroutine calculate_new_density(rho_in_scf, rho_out_scf)
!       real(dp), dimension(:), intent(in) :: rho_in_scf
!       real(dp), dimension(:), intent(out) :: rho_out_scf
!       ! ... actual calculation ...
!       rho_out_scf = rho_in_scf * 0.5_dp ! Dummy calculation
!     end subroutine calculate_new_density
!   end subroutine perform_scf
!
! end module scf_cycle_m
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only : dp, ip`: Imports precision kinds `dp` (double precision) and `ip` (integer). `ip` is not directly used in the public interface shown but might be used by private components or be a general import.
    - `fdf, only: fdf_get`: Used in the `init` procedure to read the `SCF.Mix.alpha` mixing parameter from an FDF input file.
- **External Libraries:**
    - This module has no direct dependencies on external libraries.
- **Interactions with other components:**
    - **SCF Engine (e.g., `esl_scf_m`):** The `esl_mixing_m` module is a core utility for the SCF machinery. Typically, an SCF module (like `esl_scf_m`) would contain an instance of `mixing_t`. In each iteration of the SCF loop, after computing a new output quantity (e.g., electron density), the `linear` mixing procedure is called to combine this new quantity with the quantity from the previous iteration. The result of this mixing then serves as the input for the subsequent SCF iteration.
    - **Convergence Control:** The primary role of this mixing process is to aid in achieving stable and efficient convergence of the SCF cycle. The mixing parameter `alpha` allows for tuning how aggressively the new information is incorporated, which can be crucial for preventing oscillations or divergence in challenging systems.
    - **Input Configuration (`fdf`):** The mixing behavior is user-configurable via the `SCF.Mix.alpha` parameter in the FDF input file, providing flexibility in controlling the SCF process.
    - **Data Types:** While often used for electron densities (which are `real(dp)` arrays), the `linear` mixing routine is generic enough to be applied to any set of `real(dp)` arrays of the same size.
