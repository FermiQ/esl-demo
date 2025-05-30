## Overview

The `esl_numeric_m` module is a utility library providing a collection of common numerical functions. These functions are likely used in various parts of the electronic structure application for tasks such as matrix operations, random number generation, calculation of spherical harmonics, and geometric measurements.

The current public offerings include:
-   Calculating the determinant of a 3x3 matrix (`determinant`).
-   Calculating the inverse of a 3x3 matrix (`matr3inv`).
-   Initializing the pseudo-random number generator (`init_random`).
-   Computing real spherical harmonics Y_lm and their gradients (`grylmr`).
-   Calculating the Euclidean distance between two 3D points (`distance`).

## Key Components

- **Module:** `esl_numeric_m`
    - **Description:** A collection of numerical utility routines.

### Public Functions and Subroutines

#### `determinant(a) result(d)`
- **Type:** Pure Function
- **Description:** Computes the determinant of a given 3x3 real matrix `a`.
- **Arguments:**
    - `a (real(dp), intent(in) :: a(3,3))`: The input 3x3 matrix.
- **Returns:** `d (real(dp))`: The determinant of matrix `a`.

#### `matr3inv(a) result(b)`
- **Type:** Pure Function
- **Description:** Calculates the inverse of a given 3x3 real matrix `a`. The resulting inverse matrix is stored in `b`. This function internally uses the `determinant` function.
- **Arguments:**
    - `a (real(dp), intent(in) :: a(3,3))`: The input 3x3 matrix to be inverted.
- **Returns:** `b (real(dp) :: b(3,3))`: The inverse of matrix `a`.

#### `init_random()`
- **Type:** Subroutine
- **Description:** Initializes the Fortran intrinsic pseudo-random number generator using `random_seed`. It attempts to read an integer seed value from an FDF (Flexible Data Format) input file using the key `seed`. If this key is not found, a default seed (related to the value 13) is used. The seed is processed to generate a seed array which is then supplied to `random_seed(put=seed)`.
- **Arguments:** None.
- **Returns:** Not applicable (Subroutine).

#### `grylmr(x, y, z, li, mi, ylm, grylm)`
- **Type:** Subroutine
- **Description:** Computes the value of a real spherical harmonic, Y_lm(θ, φ), at a point specified by Cartesian coordinates (x, y, z). It can also optionally compute the Cartesian gradient of the spherical harmonic at that point. The real spherical harmonics are defined as:
    -   Y_lm = c * P_lm(cos(θ)) * sin(mφ) for m < 0
    -   Y_lm = c * P_lm(cos(θ)) * cos(mφ) for m ≥ 0
    where `c` is a normalization constant and P_lm is an associated Legendre polynomial.
    The routine uses hardcoded explicit formulas for l=0, 1, and 2 for efficiency. For higher l values (up to `lmaxd=20`), it employs a general algorithm based on the `plgndr` routine from "Numerical Recipes" to calculate associated Legendre polynomials. Normalization constants are pre-calculated and stored for reuse to optimize performance.
- **Arguments:**
    - `x (real(dp), intent(in))`: The x-coordinate of the point.
    - `y (real(dp), intent(in))`: The y-coordinate of the point.
    - `z (real(dp), intent(in))`: The z-coordinate of the point.
    - `li (integer, intent(in))`: The angular momentum quantum number, l.
    - `mi (integer, intent(in))`: The magnetic quantum number, m.
    - `ylm (real(dp), intent(out))`: The calculated value of the real spherical harmonic Y_lm at (x,y,z).
    - `grylm(3) (real(dp), optional, intent(out))`: If present, this array is filled with the Cartesian components of the gradient of Y_lm (∂Y_lm/∂x, ∂Y_lm/∂y, ∂Y_lm/∂z) at (x,y,z).
- **Returns:** Not applicable (Subroutine).

#### `distance(r1, r2) result(d)`
- **Type:** Pure Function
- **Description:** Calculates the Euclidean distance between two points, `r1` and `r2`, in 3D Cartesian space. The distance is computed as `sqrt((r1(1)-r2(1))^2 + (r1(2)-r2(2))^2 + (r1(3)-r2(3))^2)`.
- **Arguments:**
    - `r1(3) (real(dp), intent(in))`: A 3-element array representing the Cartesian coordinates of the first point.
    - `r2(3) (real(dp), intent(in))`: A 3-element array representing the Cartesian coordinates of the second point.
- **Returns:** `d (real(dp))`: The Euclidean distance between points `r1` and `r2`.

## Important Variables/Constants

- **Within `grylmr` subroutine:**
    - `lmaxd (integer, parameter)`: Defined as `20`. This sets the maximum l value for which normalization constants `c` are pre-calculated and stored.
    - `tiny (real(dp), parameter)`: Defined as `1.e-30_dp`. A small floating-point number used to prevent division by zero or numerical issues when a point is very close to the origin (r=0).
    - `c(0:(lmaxd+1)*(lmaxd+1)) (real(dp), save)`: A saved array containing pre-calculated normalization constants for the real spherical harmonics, indexed by l and m.

## Usage Examples
```fortran
! TODO: Add usage example
! Example for matr3inv:
! use esl_numeric_m, only: matr3inv
! use prec, only: dp
! real(dp) :: matrix_A(3,3), matrix_B(3,3)
! matrix_A(1,:) = [1.0_dp, 2.0_dp, 3.0_dp]
! matrix_A(2,:) = [0.0_dp, 1.0_dp, 4.0_dp]
! matrix_A(3,:) = [5.0_dp, 6.0_dp, 0.0_dp]
! matrix_B = matr3inv(matrix_A)
! print *, "Inverse of A:", matrix_B

! Example for init_random and random_number:
! use esl_numeric_m, only: init_random
! use prec, only: dp
! real(dp) :: random_val
! call init_random() ! Initialize (seed might be read from FDF)
! call random_number(random_val)
! print *, "Random number:", random_val

! Example for grylmr:
! use esl_numeric_m, only: grylmr
! use prec, only: dp
! real(dp) :: x, y, z, ylm_val, grad_ylm(3)
! x = 1.0_dp; y = 0.0_dp; z = 0.0_dp ! Point on x-axis
! call grylmr(x, y, z, 1, 1, ylm_val, grad_ylm) ! Y_1,1 (p_x like)
! print *, "Y_1,1 at (1,0,0):", ylm_val, " Gradient:", grad_ylm

! Example for distance:
! use esl_numeric_m, only: distance
! use prec, only: dp
! real(dp) :: point1(3), point2(3), dist
! point1 = [0.0_dp, 0.0_dp, 0.0_dp]
! point2 = [1.0_dp, 1.0_dp, 1.0_dp]
! dist = distance(point1, point2)
! print *, "Distance:", dist
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only : dp, ip`: Imports double precision (`dp`) and integer precision (`ip`) kind specifiers from the `prec` module.
    - `fdf, only : fdf_get`: Used by `init_random` to attempt to read a random seed value from an FDF input file.
    - `esl_constants_m, only : PI`: Used by `grylmr` for the value of π, which is required for normalizing spherical harmonics.
- **External Libraries:**
    - The `grylmr` subroutine's algorithm is noted as being based on routines from "Numerical Recipes," but it does not appear to directly call any external library functions from Numerical Recipes. It's an internal implementation.
- **Interactions with other components:**
    - **`esl_geometry_m`:** The `matr3inv` function is used by `geometry_t%init` (in `esl_geometry_m`) for calculating the inverse of the unit cell matrix. The `distance` function can be used for various geometric calculations within that module or others.
    - **`esl_grid_m`:** The `grylmr` function is a key component for `grid_t%radial_function_ylm` and `grid_t%radial_function_ylm_gradient` (in `esl_grid_m`), which evaluate functions involving spherical harmonics on grid points.
    - **Basis Set Construction & Pseudopotentials:** The `grylmr` routine is essential for any part of the code that deals with angular dependencies, such as constructing atomic orbitals, defining projectors for pseudopotentials, or analyzing charge density in terms of multipole moments.
    - **Stochastic Methods:** If any part of the application uses Monte Carlo or other stochastic methods, `init_random` would be called to initialize the random number stream, followed by calls to the intrinsic `random_number()` Fortran subroutine.
    - **General Linear Algebra:** `matr3inv` and `determinant` are general-purpose 3x3 matrix utilities that could be used in various contexts requiring small matrix linear algebra (e.g., coordinate transformations, tensor manipulations).
