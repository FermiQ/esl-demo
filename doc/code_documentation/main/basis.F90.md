## Overview

The module `esl_basis_m` defines and manages basis set information for electronic structure calculations. It allows the user to choose between Plane-waves (PW) and Atomic-centered (AC) basis sets via an input option. The module is responsible for initializing the chosen basis set based on the system's geometry and providing a summary of its configuration.

## Key Components

- **Module:** `esl_basis_m`
    - **Description:** Manages the selection, initialization, and summary of the electronic structure basis set (either Plane-waves or Atomic-centered orbitals).
- **Type:** `basis_t`
    - **Description:** A derived type that encapsulates all information about the chosen basis set. It contains an indicator for the basis type (PW or AC) and holds the actual basis data in nested types (`basis_pw_t` or `basis_ac_t`).
    - **Components:**
        - `type (integer)`: Stores whether the basis is `PLANEWAVES` or `ATOMCENTERED`.
        - `pw (type(basis_pw_t))`: Contains the plane-wave basis information if `type` is `PLANEWAVES`.
        - `ac (type(basis_ac_t))`: Contains the atomic-centered basis information if `type` is `ATOMCENTERED`.
    - **Procedures:**
        - `init`: Initializes the basis set.
        - `summary`: Prints a summary of the basis set configuration.

### Subroutine/Function: `init(this, geo)`

- **Description:** Initializes the `basis_t` object. It reads the `BasisSet` parameter from an FDF input file (defaulting to 'Planewaves'). Based on this choice, it initializes either the plane-wave (`this%pw`) or atomic-centered (`this%ac`) basis components using their respective `init` procedures. If the `WITH_FLOOK` preprocessor flag is enabled, it also registers the plane-wave cutoff energy with a dictionary.
- **Arguments:**
    - `this (class(basis_t))`: The `basis_t` object to be initialized.
    - `geo (type(geometry_t), intent(in))`: Input geometry information, necessary for basis initialization.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `summary(this)`

- **Description:** Prints a summary of the initialized basis set. It calls the `summary` procedure of either the `pw` or `ac` component depending on the selected `this%type`.
- **Arguments:**
    - `this (class(basis_t), intent(in))`: The `basis_t` object whose summary is to be printed.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `ndims_from_spacing(ndims, spacing, cell)`

- **Description:** Calculates the required dimensions for a Fast Fourier Transform (FFT) grid (`ndims`) based on a desired real-space `spacing` and the simulation `cell` vectors. For each dimension, it determines the number of grid points by dividing the cell vector length by the spacing and then calls `fourier_dim` to get a suitable FFT dimension.
- **Arguments:**
    - `ndims (integer, intent(out) :: ndims(3))`: Output array containing the calculated FFT grid dimensions for each of the 3 spatial directions.
    - `spacing (real(dp), intent(in))`: The desired real-space grid spacing.
    - `cell (real(dp), intent(in) :: cell(3,3))`: The simulation cell lattice vectors. `cell(idim, idim)` is assumed to be the length of the cell in the idim-th direction for this calculation.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants

- `PLANEWAVES (integer, public, parameter)`: A constant with value `1`, used to identify the plane-wave basis set type.
- `ATOMCENTERED (integer, public, parameter)`: A constant with value `2`, used to identify the atomic-centered basis set type.

## Usage Examples

```fortran
! TODO: Add usage example
! Example for init:
! use esl_basis_m
! use esl_geometry_m
! type(basis_t) :: basis
! type(geometry_t) :: geo
! ! ... initialize geo ...
! call basis%init(geo)
!
! Example for summary:
! call basis%summary()
!
! Example for ndims_from_spacing:
! integer :: fft_dims(3)
! real(dp) :: grid_spacing, unit_cell(3,3)
! ! ... define grid_spacing and unit_cell ...
! call ndims_from_spacing(fft_dims, grid_spacing, unit_cell)
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec`: Used for defining double precision (`dp`) real kinds.
    - `fdf, only : fdf_get, leqi`: Used to read the `BasisSet` choice from an FDF (Flexible Data Format) input file.
    - `esl_message_m`: Used for error handling, specifically `message_error` is called for unknown basis sets.
    - `esl_basis_ac_m`: Provides the `basis_ac_t` type and its associated procedures (e.g., `init`, `summary`) for the atomic-centered basis.
    - `esl_basis_pw_m`: Provides the `basis_pw_t` type and its associated procedures (e.g., `init`, `summary`) for the plane-wave basis.
    - `esl_constants_m`: [How it's used - Likely provides physical or mathematical constants, though not directly shown in use in the provided snippet.]
    - `esl_geometry_m`: Provides the `geometry_t` derived type, which is an input to the `init` subroutine.
    - `esl_grid_m`: Likely provides the `fourier_dim` subroutine used in `ndims_from_spacing` for determining FFT-friendly dimensions.
    - `yaml_output`: [How it's used - Potentially for structured data output, but not explicitly used in the visible code.]
    - `dictionary` (conditional, if `WITH_FLOOK` is defined): Used if the FLOOK library is enabled, for dictionary functionalities.
    - `esl_dict_m` (conditional, if `WITH_FLOOK` is defined): Provides procedures like `esl_dict_var_add` for adding variables to a dictionary when FLOOK is enabled.
- **External Libraries:**
    - `fdf` could be considered an external input parsing library.
- **Interactions with other components:**
    - The `esl_basis_m` module is fundamental for defining the representation of wavefunctions and densities. It would be directly used by modules that set up or operate on these quantities, such as Hamiltonian construction modules (`hamiltonian_pw_m`, `hamiltonian_ac_m`) and density calculation modules. The choice of basis also dictates which computational kernels are subsequently called.
