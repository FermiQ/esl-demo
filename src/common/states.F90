module esl_states_m

  use prec, only : dp,ip
  use fdf
  use yaml_output

  use esl_geometry_m

  implicit none
  private

  public :: states_t

  type wfn_t
     real(dp),    allocatable :: dcoef(:) !<Coefficients of the wavefunction in the basis
     complex(dp), allocatable :: zcoef(:)
  end type wfn_t


  !Data structure for the states
  type states_t
     integer :: nstates
     integer :: nspin
     integer :: nkpt
     real(dp) :: nel   !< Number of electrons

     logical :: complex_states
     integer :: ncoef !< Number of coefficients

     type(wfn_t), allocatable :: states(:,:,:)  !nstates, nspin, nkpt
     real(dp), allocatable :: occ_numbers(:,:,:)
     real(dp), allocatable :: k_weights(:)
     real(dp), allocatable :: eigenvalues(:,:,:)
   contains
     private
     procedure, public :: init
     procedure, public :: summary
     procedure, public :: randomize
     final :: cleanup
  end type states_t

contains

  !Initialize the states
  !----------------------------------------------------
  subroutine init(this, geo, nspin, nkpt, complex, ncoef)
    class(states_t)  :: this
    type(geometry_t), intent(inout) :: geo
    integer,          intent(in)    :: nspin
    integer,          intent(in)    :: nkpt
    logical,          intent(in)    :: complex !< Should the wavefunctions be complex?
    integer,          intent(in)    :: ncoef    !< Size of wavefunctions (number of coefficients)
    
    integer :: ist, isp, ik, extra_states
    real(dp) :: nel, occ
    
    extra_states = fdf_get('ExtraStates', 0)

    !TODO: the charge hould be a real number, not an integer
    this%nspin = nspin
    this%nel = geo%electronic_charge()
    this%nstates = int(geo%electronic_charge()/2)
    if ( this%nstates*2 < geo%electronic_charge() ) this%nstates = this%nstates + 1
    this%nstates = this%nstates + extra_states
    this%nkpt = nkpt
    this%complex_states = complex
    this%ncoef = ncoef

    allocate(this%states(1:this%nstates, 1:nspin, 1:nkpt))
    allocate(this%occ_numbers(1:this%nstates, 1:nspin, 1:nkpt))
    this%occ_numbers(1:this%nstates, 1:nspin, 1:nkpt) = 0._dp
    nel = this%nel
    do ist = 1, this%nstates
      if ( nel == 0.0_dp) exit
      occ = 2.0/real(this%nspin, dp)
      if ( nel < occ ) occ = nel
      this%occ_numbers(ist,  1:nspin, 1:nkpt) = occ
      nel = nel - occ
    end do

    allocate(this%eigenvalues(1:this%nstates, 1:nspin, 1:nkpt))
    this%eigenvalues(1:this%nstates, 1:nspin, 1:nkpt) = 0._dp    

    allocate(this%k_weights(1:nkpt))
    this%k_weights(:) = 1.d0/this%nkpt

    if (this%complex_states) then
      do ik = 1, nkpt
        do isp = 1, nspin
          do ist = 1, this%nstates
            allocate(this%states(ist, isp, ik)%zcoef(1:this%ncoef))
          end do
        end do
      end do
    else
      do ik = 1, nkpt
        do isp = 1, nspin
          do ist = 1, this%nstates
            allocate(this%states(ist, isp, ik)%dcoef(1:this%ncoef))
          end do
        end do
      end do
    end if

  end subroutine init


  !Release the states
  !----------------------------------------------------
  subroutine cleanup(this)
    type(states_t):: this

    integer :: ist, isp, ik

    if (allocated(this%states)) then
      do ik = 1, this%nkpt
        do isp = 1, this%nspin
          do ist = 1, this%nstates
            if(allocated(this%states(ist, isp, ik)%dcoef)) &
              deallocate(this%states(ist, isp, ik)%dcoef)
            if(allocated(this%states(ist, isp, ik)%zcoef)) &
              deallocate(this%states(ist, isp, ik)%zcoef)
          end do
        end do
      end do
      deallocate(this%states)
    end if

    if(allocated(this%occ_numbers)) deallocate(this%occ_numbers)
    if(allocated(this%eigenvalues)) deallocate(this%eigenvalues)
    if(allocated(this%k_weights)) deallocate(this%k_weights)

  end subroutine cleanup


  !Randomize the states
  !----------------------------------------------------
  subroutine randomize(this)
    class(states_t):: this

    integer :: ist, isp, ik
    real(dp), allocatable :: tmp_re(:), tmp_im(:)
    real(dp) :: norm

    if(this%complex_states) then
      allocate(tmp_re(1:this%ncoef))
      allocate(tmp_im(1:this%ncoef))
      do ik = 1, this%nkpt
        do isp = 1, this%nspin
          do ist = 1, this%nstates
            call random_number(tmp_re(1:this%ncoef))
            call random_number(tmp_im(1:this%ncoef))
            this%states(ist, isp, ik)%zcoef(1:this%ncoef) = tmp_re(1:this%ncoef) &
                   + cmplx(0.d0,1.d0,kind=dp)*tmp_im(1:this%ncoef) 
            norm = sum(abs(this%states(ist, isp, ik)%zcoef(1:this%ncoef))**2)
            this%states(ist, isp, ik)%zcoef(1:this%ncoef) = this%states(ist, isp, ik)%zcoef(1:this%ncoef)/sqrt(norm)
          end do
        end do
      end do
      deallocate(tmp_re, tmp_im)
    else 
      do ik = 1, this%nkpt
        do isp = 1, this%nspin
          do ist = 1, this%nstates
            call random_number(this%states(ist, isp, ik)%dcoef(1:this%ncoef))
            norm = sum(abs(this%states(ist, isp, ik)%dcoef(1:this%ncoef))**2)
            this%states(ist, isp, ik)%dcoef(1:this%ncoef) = this%states(ist, isp, ik)%dcoef(1:this%ncoef)/sqrt(norm)
          end do
        end do
      end do
    end if

  end subroutine randomize

  !Summary
  !----------------------------------------------------
  subroutine summary(this)
    class(states_t) :: this

    call yaml_mapping_open("States")
    call yaml_map("Number of states", this%nstates)
    call yaml_map("Spin components", this%nspin)
    call yaml_map("Number of k-points", this%nkpt)
    call yaml_mapping_close()

  end subroutine summary

end module esl_states_m
