## Overview

The `esl_mulliken_ac_m` module provides a subroutine, `mulliken_ac_summary`, for performing Mulliken population analysis and printing the results. Mulliken analysis is a method used to assign partial atomic charges to atoms within a molecule or solid, based on the distribution of electrons among the basis orbitals. This routine is specifically designed for systems described by atomic-centered (AC) basis sets and utilizes sparse matrix representations for the overlap (S) and density (DM) matrices.

The calculated "orbital populations" and summed "atomic charges" are then output in a structured YAML format.

*Note: The formula used in the provided source code for calculating orbital populations ($F(io) = \sum_{jo} S_{io,jo} \cdot DM_{io,jo}$, an element-wise product summed over connected elements in a row) deviates from the standard Mulliken population analysis formula, which typically involves matrix products like $(P \cdot S)_{ii}$ or related expressions. The documentation below describes the implemented behavior.*

## Key Components

- **Module:** `esl_mulliken_ac_m`
    - **Description:** Calculates and summarizes Mulliken-like atomic and orbital charges for systems using atomic-centered basis sets.

### Public Subroutine

#### `mulliken_ac_summary(basis, S, DM)`
- **Description:** This subroutine calculates and prints a Mulliken-like population analysis.
    1.  **Initialization:**
        *   It retrieves the sparsity pattern (`sp`) from the input overlap matrix `S`. It is assumed that the density matrix `DM` shares a compatible sparsity for the element-wise operations performed.
        *   Allocates two real arrays: `M` of size `basis%n_site` (to store summed charge per atomic site) and `F` of size `basis%n_orbital` (to store population per orbital).
    2.  **Orbital Population Calculation:**
        *   It iterates through each atomic site `ia` and then through each orbital `io` centered on that site.
        *   For each orbital `io`, it calculates a quantity `F(io)`. This is done by iterating through the non-zero elements in the row `io` of the sparsity pattern `sp`. If `ind` is the index of such a non-zero element (representing matrix element at `(io, col(ind))`), then `F(io)` accumulates the sum `S%M(ind) * DM%M(ind)`. Effectively, $F(io) = \sum_{j} S_{io,j} \cdot DM_{io,j}$, where the sum is over columns `j` for which $S_{io,j}$ (and $DM_{io,j}$) are non-zero according to the pattern.
    3.  **Site Population Calculation:**
        *   After calculating `F(io)` for all orbitals `io` on site `ia`, it sums these `F(io)` values to get the total population `M(ia)` for that site.
    4.  **Output Generation:**
        *   The results are printed in YAML format under a main mapping "Mulliken".
        *   The total sum of `M(ia)` over all sites is printed as 'Total'.
        *   For each atomic site `ia`:
            *   Its site index is printed.
            *   The summed population `M(ia)` for that site is printed as 'Sum'.
            *   The individual `F(io)` values for all orbitals `io` belonging to site `ia` are printed as 'Individual'.
        *   (A commented-out line in the source suggests an intent to print the species label for each site, but this is not currently active.)
    5.  **Cleanup:** Deallocates the temporary arrays `M` and `F`.
- **Arguments:**
    - `basis (type(basis_ac_t), intent(in))`: The atomic-centered basis set definition. This is used to access information about atomic sites, the number of orbitals, and the mapping between orbitals and sites (e.g., `basis%n_site`, `basis%n_orbital`, `basis%site_orbital_start`).
    - `S (type(sparse_matrix_t), intent(in))`: The sparse overlap matrix $S_{\mu\nu} = \langle \psi_\mu | \psi_\nu \rangle$.
    - `DM (type(sparse_matrix_t), intent(in))`: The sparse density matrix $P_{\mu\nu}$.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants
- This module does not define any public Fortran `parameter` constants. The key inputs are the `basis`, `S` matrix, and `DM`.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage after an SCF calculation:
!
! use esl_mulliken_ac_m, only: mulliken_ac_summary
! use esl_basis_ac_m, only: basis_ac_t
! use esl_sparse_matrix_m, only: sparse_matrix_t
! use esl_density_ac_m, only: density_ac_t ! Assuming DM is part of density_ac_t
!
! type(basis_ac_t) :: my_ac_basis
! type(density_ac_t) :: my_ac_density ! Contains my_ac_density%DM
! ! my_ac_basis%S would hold the overlap matrix
!
! ! ... (Initialize my_ac_basis, and perform SCF to get my_ac_density%DM and my_ac_basis%S) ...
!
! ! Perform and print Mulliken analysis
! call mulliken_ac_summary(my_ac_basis, my_ac_basis%S, my_ac_density%DM)
!
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec`: Provides precision kinds (e.g., `dp`).
    - `esl_basis_ac_m`: Provides the `basis_ac_t` type, which defines the atomic-centered basis set and gives access to site/orbital mappings.
    - `esl_sparse_pattern_m`: Provides `sparse_pattern_t`, used via `S%sp` to iterate through non-zero elements.
    - `esl_sparse_matrix_m`: Provides `sparse_matrix_t` for the overlap matrix `S` and density matrix `DM`.
    - `yaml_output`: Used for formatting the output of the Mulliken analysis as YAML.
- **External Libraries:**
    - This module does not have direct dependencies on external libraries.
- **Interactions with other components:**
    - **`esl_density_ac_m`:** The `density_ac_t` type defined in `esl_density_ac_m` also contains a `mulliken_summary` subroutine with an identical signature. It is highly probable that one implementation calls the other or that they are intended to be the same piece of logic. The routine in `esl_density_ac_m` might be the primary interface point.
    - **SCF Solvers (e.g., `esl_scf_m`):** The `mulliken_ac_summary` routine would typically be called after an SCF calculation has converged. The converged density matrix (`DM`) and the system's overlap matrix (`S`) in the given `basis` are required inputs.
    - **Output/Analysis Stage:** This routine serves as an analysis tool to provide insights into the charge distribution within the simulated system by assigning partial charges to atoms. The YAML output is suitable for logging and parsing by other tools.
    - **Basis Set (`basis_ac_t`) and Matrices (`sparse_matrix_t`):** The accuracy and interpretation of Mulliken charges are highly dependent on the choice of basis set. The routine requires a fully defined `basis_ac_t` object, the corresponding overlap matrix `S`, and the density matrix `DM`.
</tbody>
</table>
