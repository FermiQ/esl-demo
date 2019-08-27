module esl_system_m

  use prec, only : dp, ip

#ifdef WITH_FLOOK
  use dictionary
  use esl_dict_m
#endif

  ! Sparse pattern for LO
  use esl_sparse_pattern_m, only: sparse_pattern_t
  ! Sparse matrix for LO (only real(dp))
  use esl_sparse_matrix_m, only: sparse_matrix_t

  use esl_basis_m
  use esl_energy_m
  use esl_force_m
  use esl_geometry_m
  use esl_ion_interaction_m
  use esl_smear_m
  use esl_states_m
  use esl_species_m
  
  implicit none
  
  private

  public :: system_t

  !< System data structure
  !<
  !< The system data-structure is a basic container that encompass
  !< everything related to a single DFT calculation.
  !< It contains the basis, energy, geometry and force information.
  !< This is the highest level type which is necessary for performing
  !< any DFT calculation.
  type system_t
    
    !< Container for the basis used (container for PW and AC)
    type(basis_t) :: basis
    !< All system related energies
    type(energy_t) :: energy
    !< Geometry information
    type(geometry_t) :: geo
    !< Force inforamtion (connected to geo)
    type(force_t) :: force
    !< Ion-ion interaction (connected to geo)
    type(ion_interaction_t) :: ion_inter

    ! TODO decide whether the system should be inherited for the LO/PW
    ! case. It may make the system a lot easier to figure out.
    ! However, it will prohibit switching from PW/LO -> LO/PW within the same
    ! calculation.

    ! LO dependent variables
    type(sparse_pattern_t):: sparse_pattern
    ! Perhaps these should be transferred to the basis_ac type.
    ! However, they could be defined higher level since the basis
    ! does not necessarily use the overlap matrix
    type(sparse_matrix_t) :: S ! always 1D

    !< Total number of electrons
    real(dp) :: nElectrons = 0._dp
  contains
    private
    procedure, public :: init
    procedure, public :: update
    procedure, public :: summary

    final  :: cleanup
  end type system_t

contains

  !< Initialize all contained types.
  !<
  !< If your compilation is performed with Lua, it is only *after*
  !< this call are the listed items exposed to Lua:
  !<   - Geometry.xyz
  !<   - Geometry.cell
  !<   - Geometry.Force.Total
  !<   - Geometry.Force.Local
  !<   - Geometry.Force.NonLocal
  !<   - Geometry.Force.IonIon
  !<   - IonIon.Ewald.Alpha
  !<   - IonIon.Ewald.Alpha
  !<   - E.Total
  !<   - E.Hartree
  !<   - E.Fermi
  !<   - E.IonIon
  !<   - E.External
  !<   - E.Exchange
  !<   - E.Correlation
  !<   - E.Kinetic
  subroutine init(this)
    class(system_t), intent(inout) :: this

    call this%energy%init()
    this%geo = geometry_t()
    call this%basis%init(this%geo)
    call this%force%init(this%geo%n_atoms)
    call this%ion_inter%init()
    call this%energy%init()

    call this%summary()

#ifdef WITH_FLOOK
    ! Add variables
    call esl_dict_var_add('Geometry.xyz', this%geo%xyz)
    call esl_dict_var_add('Geometry.cell', this%geo%cell)
    call esl_dict_var_add('Geometry.Force.Total', this%force%total)
    call esl_dict_var_add('Geometry.Force.Local', this%force%loc)
    call esl_dict_var_add('Geometry.Force.NonLocal', this%force%nl)
    call esl_dict_var_add('Geometry.Force.IonIon', this%force%ionion)
    call esl_dict_var_add('IonIon.Ewald.Alpha', this%ion_inter%alpha)
    
    call esl_dict_var_add('E.Total', this%energy%total)
    call esl_dict_var_add('E.Hartree', this%energy%hartree)
    call esl_dict_var_add('E.Fermi', this%energy%fermi)
    call esl_dict_var_add('E.IonIon', this%energy%ionion)
    call esl_dict_var_add('E.External', this%energy%extern)
    call esl_dict_var_add('E.Exchange', this%energy%exchange)
    call esl_dict_var_add('E.Correlation', this%energy%correlation)
    call esl_dict_var_add('E.Kinetic', this%energy%kinetic)

#endif

  end subroutine init

  !< Update the system after an MD step
  !< Currently this re-calculates ion-specific quantities that are
  !< independent on the SCF.
  subroutine update(this, periodic)
    class(system_t), intent(inout) :: this
    logical, intent(in) :: periodic

    if ( periodic ) then
      call this%ion_inter%calculate_periodic(this%geo, this%force%ionion, this%energy%ionion)
    else
      call this%ion_inter%calculate_isolated(this%geo, this%force%ionion, this%energy%ionion)
    end if
    
  end subroutine update


  !< Clean up all the contained types.
  subroutine cleanup(sys)
    type(system_t) :: sys

  end subroutine cleanup

  !< Show information about the type as YAML output
  subroutine summary(sys)
    use yaml_output
    class(system_t) :: sys

    call yaml_mapping_open("System")
    call sys%geo%summary()
    call sys%basis%summary()
    call yaml_mapping_close()

  end subroutine summary

end module esl_system_m
