## Overview

The `esl_basis_pw_m` module defines and manages plane-wave (PW) basis sets for electronic structure calculations. It introduces the `basis_pw_t` derived type, which extends `basis_base_t` (from `esl_basis_base_m`). A PW basis set is characterized by an energy cutoff (`ecut`) that determines the maximum kinetic energy (and thus the shortest wavelength or highest G-vector magnitude) of the plane waves included.

This module is responsible for:
-   Initializing the PW basis based on this energy cutoff and the system's geometry (cell vectors). This includes determining appropriate dimensions for the Fast Fourier Transform (FFT) grid.
-   Storing key properties of the PW basis, such as the energy cutoff, the FFT grid dimensions (`ndims`), the reciprocal space metric tensor (`gmet`), the squared magnitudes of the included G-vectors (`gmod2`), and a map from a PW index to the Miller indices of the G-vectors (`gmap`).
-   The total number of plane waves in the basis set is stored in the inherited `size` component.
-   It also manages an auxiliary real-space grid (inherited as `this%grid`) consistent with the PW cutoff.

## Key Components

- **Module:** `esl_basis_pw_m`
    - **Description:** Defines and manages plane-wave (PW) basis sets, including setup based on an energy cutoff and cell geometry.

- **Type:** `basis_pw_t`, extends `basis_base_t`
    - **Description:** Represents a plane-wave basis set.
    - **Inherited Components (from `basis_base_t`):**
        - `grid (type(grid_t))`: An auxiliary real-space grid compatible with the PW basis.
        - `size (integer)`: The total number of plane waves (G-vectors + k) included in the basis set, satisfying the energy cutoff condition.
    - **Specific Components:**
        - `ecut (real(dp))`: The kinetic energy cutoff for plane waves, given in Hartree units. Only plane waves $\mathbf{G}$ such that $\frac{1}{2} |\mathbf{G}+\mathbf{k}|^2 \le E_{cut}$ are included (for a given k-point $\mathbf{k}$).
        - `ndims(3) (integer)`: The dimensions (nx, ny, nz) of the underlying FFT grid used for representing wavefunctions and densities in real and reciprocal space. These are determined based on `ecut`.
        - `gmet(3,3) (real(dp))`: The reciprocal space metric tensor, defined as $gmet_{ij} = \mathbf{b}_i \cdot \mathbf{b}_j$, where $\mathbf{b}_i$ are the reciprocal lattice vectors.
        - `gmod2 (real(dp), allocatable :: gmod2(:))`: A 1D array of size `this%size`, storing the squared magnitude of each G-vector (plus k-vector, i.e., $|\mathbf{G}+\mathbf{k}|^2$) included in the basis set.
        - `gmap (integer, allocatable :: gmap(:,:))`: A 2D array of shape `(3, this%size)`, mapping each plane-wave index (from 1 to `this%size`) to its corresponding Miller indices (h, k, l) which define the G-vector $\mathbf{G} = h\mathbf{b}_1 + k\mathbf{b}_2 + l\mathbf{b}_3$.
    - **Procedures (Public, Type-Bound):**
        - `init(this, geo)`: Initializes the plane-wave basis set.
        - `summary(this)`: Prints a summary of the PW basis properties.
    - **Finalizer (Type-Bound):**
        - `finalizer(this)`: Deallocates the `gmod2` and `gmap` arrays.

### Type-Bound Procedures

#### `init(this, geo)`
- **Description:** Initializes the `basis_pw_t` object using the system geometry (`geo` of type `geometry_t`).
    1.  Reads the plane-wave energy cutoff `this%ecut` from an FDF input file (key: `cut-off`, default: 10 Ha).
    2.  Calls the private module subroutine `ndims_from_ecut` to determine suitable FFT grid dimensions `this%ndims` based on `this%ecut`, the inverse cell vectors `geo%icell`, and a k-point (currently hardcoded to $\mathbf{k}=[0,0,0]$ for this determination step).
    3.  Initializes the inherited auxiliary real-space grid `this%grid` using these `this%ndims` and the direct cell vectors `geo%cell`.
    4.  Calculates the reciprocal space metric tensor `this%gmet` from `geo%icell`.
    5.  Calls `get_number_of_pw` (from `esl_utils_pw_m`) to determine `this%size`, which is the count of G-vectors satisfying the energy cutoff for the calculated `this%ndims` and `this%gmet` (again, for $\mathbf{k}=[0,0,0]$).
    6.  Allocates `this%gmod2` (size `this%size`) and `this%gmap` (shape `(3, this%size)`).
    7.  Calls `construct_mod_map_tables` (from `esl_utils_pw_m`) to populate `this%gmod2` with $|\mathbf{G}|^2$ values and `this%gmap` with the Miller indices for each G-vector in the basis (for $\mathbf{k}=[0,0,0]$).
- **Arguments:**
    - `this (class(basis_pw_t), intent(inout))`: The PW basis object to initialize.
    - `geo (type(geometry_t), intent(in))`: The system geometry, providing cell and inverse cell matrices.

#### `summary(this)`
- **Description:** Prints a summary of the plane-wave basis set configuration to standard output using YAML formatting. The summary includes:
    -   Type ("Plane waves").
    -   Energy cut-off (`this%ecut`) in Hartree.
    -   Total number of plane waves (`this%size`).
    It also calls `this%grid%summary()` to print details of the associated auxiliary real-space grid.
- **Arguments:**
    - `this (class(basis_pw_t), intent(in))`: The PW basis object to summarize.

#### `finalizer(this)`
- **Description:** This is a type-bound `final` subroutine that is automatically invoked when a `basis_pw_t` object is deallocated. It deallocates the `this%gmod2` and `this%gmap` arrays if they were previously allocated.
- **Arguments:**
    - `this (type(basis_pw_t))`: The PW basis object being finalized.

### Private Module Subroutine

#### `ndims_from_ecut(ndims_out, ecut_in, gcell_in, kpt_in)`
- **Description:** Calculates suitable FFT grid dimensions (`ndims_out`) that can accurately represent plane waves up to a given kinetic energy cutoff `ecut_in`. It considers the reciprocal cell vectors (`gcell_in`, which is `geo%icell`) and a specific k-point (`kpt_in`). The routine iteratively increases grid dimensions (using `fourier_dim` for FFT-friendliness) until the Nyquist frequency along each direction is sufficient for the given `ecut_in`. It uses internal helper functions `smallest` (to find the limiting G-vector component for a given grid size) and `dsq` (to calculate $|\mathbf{G}+\mathbf{k}|^2$).
- **Arguments:**
    - `ndims_out (integer, intent(out) :: ndims(3))`: The calculated FFT grid dimensions.
    - `ecut_in (real(dp), intent(in))`: The target energy cutoff.
    - `gcell_in (real(dp), intent(in) :: gcell(3,3))`: The reciprocal lattice vectors (inverse of the direct cell matrix).
    - `kpt_in (real(dp), intent(in) :: kpt(3))`: The k-point vector (in fractional reciprocal coordinates).

## Important Variables/Constants
- **`ecut (real(dp))`**: The kinetic energy cutoff defining the size of the PW basis.
- **`ndims(3) (integer)`**: Dimensions of the FFT grid.
- **`gmet(3,3) (real(dp))`**: Reciprocal space metric tensor.
- **`gmod2(:)`, `gmap(:,:)`**: Arrays storing the squared magnitudes and Miller indices of the G-vectors in the basis.
- **`size (integer)`**: (Inherited) Total number of plane waves.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_basis_pw_m
! use esl_geometry_m
! type(basis_pw_t) :: pw_basis
! type(geometry_t) :: system_geo
!
! ! ... (Initialize system_geo, e.g., by reading from input) ...
!
! ! Initialize the plane-wave basis
! call pw_basis%init(system_geo)
!
! ! Print a summary of the basis
! call pw_basis%summary()
!
! ! ... (pw_basis can now be used by Hamiltonian and density modules) ...
!
! ! pw_basis will be finalized automatically when it goes out of scope.
```

## Dependencies and Interactions

- **Base Type:** `esl_basis_base_m`: `basis_pw_t` extends `basis_base_t`, inheriting `grid` and `size`.
- **Core ESL Modules:**
    - `prec`: For precision kinds (`dp`).
    - `yaml_output`: Used by `summary` for formatted output.
    - `fdf`: Used by `init` to read the `cut-off` parameter.
    - `esl_constants_m`: Provides `PI`, used in `ndims_from_ecut`.
    - `esl_geometry_m`: Provides `geometry_t` (cell and inverse cell matrices) needed for `init` and `ndims_from_ecut`.
    - `esl_grid_m`: The `grid_t` type (inherited `this%grid`) is initialized and managed.
- **PW Specific Utility Modules:**
    - `esl_utils_pw_m`: Provides `get_number_of_pw` and `construct_mod_map_tables` which are essential for determining the basis size and populating G-vector information.
- **FFT Utility Modules:**
    - `module_fft_sg` (via `esl_grid_m` or directly in `ndims_from_ecut`): Provides `fourier_dim` for ensuring FFT-friendly grid dimensions.
- **Interactions with other components:**
    - **System Setup (`esl_system_m`, `esl_main_m`):** If a plane-wave basis is chosen for a calculation, an instance of `basis_pw_t` would be created and initialized, typically as part of an overarching `system_t` or `basis_t` object.
    - **Hamiltonian Construction (`esl_hamiltonian_pw_m`):** This module would heavily rely on `basis_pw_t` to:
        - Know the set of G-vectors (`gmap`, `gmod2`, `size`).
        - Access the FFT grid (`grid`) for transformations between real and reciprocal space.
        - Use `ecut` for consistency.
    - **Density Calculation (`esl_density_pw_m`):** Represents the electron density using PW coefficients corresponding to the G-vectors in `basis_pw_t`. It uses the FFT grid for real/reciprocal space transforms.
    - **Input Configuration (`fdf`):** The primary parameter `ecut` is read from an FDF input file, allowing user control over the basis set size and accuracy.
</tbody>
</table>
