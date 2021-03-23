!===================================================================================================
! FatesVegetationManagementMod.F90
! Implemented by: Joshua M. Rady
! Started: 8/26/2020
!
! This module contains subroutines to simulate aspects of human management of vegetation.  While
! written to implement forest management the code contains generic behaviors that could be used to
! implement other activities, e.g. agriculture.
!
! Note: This code should be revised to conform to FATES coding style, to the extent that it exists.
! The coding style is a bit vague with a mix of 2 and 3 space identing.  I will use 2 since they
! indent starting at the module scope.
!
!===================================================================================================

module FatesVegetationManagementMod
  
  use EDTypesMod, only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDPftvarcon, only : EDPftvarcon_inst
  use FatesAllometryMod, only : h2d_allom, h_allom
  use FatesConstantsMod, only : pi_const
  use FatesConstantsMod, only : r8 => fates_r8
  ! Using generic integers for PFT numbers should be fine but we follow FatesAllometryMod in
  ! explicitly sizing them.
  use FatesConstantsMod, only : i4 => fates_int
  use FatesConstantsMod, only : itrue, ifalse
  use FatesGlobals, only : endrun => fates_endrun 
  use FatesInterfaceTypesMod, only : bc_in_type
  use PRTGenericMod, only : prt_vartypes
  use PRTParametersMod, only : prt_params
  
  ! Log and error reporting:
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use FatesGlobals, only : fates_log
  
  ! Enforce explicit type declarations:
  implicit none
  private
  
  ! Entry Points:
  public :: vegetation_management_init
  public :: managed_mortality
  public :: anthro_mortality_rate
  public :: managed_fecundity
  public :: management_fluxes
  public :: spawn_anthro_disturbed_cohorts
  ! Destructor!!!!!
  
  ! Primitives:
  public :: plant
  public :: init_temporary_cohort ! Currently doesn't need to be public.
  
  private :: kill
  private :: kill_disturbed
  private :: kill_patch
  
  ! Management Operations:
  private :: understory_control
  private :: thin_row_low
  private :: thinnable_patch
  !private :: thin_low
  private :: thin_proportional
  private :: thin_patch_low_perfect
  private :: thin_low_perfect
  private :: thin_patch_low_probabilistic
  private :: thin_low_probabilistic
  private :: harvest_mass_min_area
  private :: plant_harvestable_biomass
  private :: cohort_harvestable_biomass
  private :: clearcut_patch
  private :: clearcut
  
  ! Utilities:
  private :: cohort_effective_basal_area
  private :: patch_effective_basal_area
  private :: cohort_effective_n
  private :: patch_effective_n
  private :: cohort_disturbed_n
  private :: patch_disturbed_n
  private :: cohort_disturbed_basal_area
  public :: patch_disturbed_basal_area ! Temporarily public
  private :: get_flux_profile
  private :: validate_size_specifications

  ! Prescribed Event Driver File:
  private :: load_prescribed_events
  private :: field_pop
  private :: field_pop_int
  private :: field_pop_real
  private :: is_now
  private :: is_here

  ! Interfaces:
  ! JMR_NOTE: These interfaces have given me a bit of trouble and are a work in progress:
  interface effective_basal_area
    module procedure patch_effective_basal_area
    module procedure cohort_effective_basal_area
  end interface
  
  interface effective_n
    module procedure patch_effective_n
    module procedure cohort_effective_n
  end interface
  
  interface effective_stem_density
    module procedure patch_effective_n
    module procedure cohort_effective_n
  end interface
  
  interface disturbed_n
    module procedure cohort_disturbed_n
    module procedure patch_disturbed_n
  end interface
  
!   interface disturbed_stem_density
!     module procedure cohort_disturbed_n
!     module procedure patch_disturbed_n
!   end interface
  
  interface disturbed_basal_area
    module procedure cohort_disturbed_basal_area
    module procedure patch_disturbed_basal_area
  end interface
  
  ! Globals:
  
  ! Debugging flag / switch for module:
  logical, parameter, private :: debug = .true.
  
  ! PFTs:
  ! These class definitions should be determined dynamically. The following definitions assume the
  ! default 16 pft parameter set.  We need a better way to get the classes of PFTs. The woody flag
  ! (EDPftvarcon_inst%woody) includes shrubs.
  !integer(i4), allocatable, dimension(:), target, private :: woody_pfts
  !integer(i4), parameter, private :: woody_pfts(9) = [1,2,3,4,5,6,7,8,9]
  !integer(i4), parameter, private :: tree_pfts(6) = [1,2,3,4,5,6]
  integer(i4), target, private :: woody_pfts(9) = [1,2,3,4,5,6,7,8,9]
  integer(i4), target, public :: tree_pfts(6) = [1,2,3,4,5,6] ! Temporarily public
  ! Shrubs?
  integer(i4), target, private :: understory_pfts(6) = [7,8,9,10,11,12]
  
  ! Understory control modes:
  ! Note: The order an value of these have no significance and are subject to change.  Use the
  ! names.  May change to an enumerated type. Names may change.
  integer, parameter, private :: method_mow = 1 ! or method_cut
  integer, parameter, private :: method_herbicide = 2
  integer, parameter, private :: method_burn = 3
  
  ! Mode / flux profiles:
  ! There are may be more practices than flux or harvest profiles.  logging_traditional -> bole_harvest
  integer, parameter, private :: null_profile = 0
  integer, parameter, private :: logging_traditional = 1 ! logging_module, logging_legacy, logginng_classic
  integer, parameter, private :: bole_harvest = 2 ! harvest_bole ! better namespacing
  integer, parameter, private :: in_place = 3
  
  ! String length specifier: (After FatesInventoryInitMod)
  ! This is longer than what Fortran will likely allow for a line from read(). Consider shortening.
  integer, parameter :: line_strlen = 512
  
  ! VM Events:--------------------------------------------------------------------------------------
  type, private :: vm_event
      integer(i4) :: code ! The event code.
      integer, dimension(1) :: pfts ! Change to array (16 in length?)
      real(r8) :: density ! planting_density?
      real(r8) :: dbh
      real(r8) :: height
      real(r8) :: row_fraction
      real(r8) :: final_basal_area
      real(r8) :: thin_fraction
      real(r8) :: dbh_min
      real(r8) :: ht_min
      real(r8) :: patch_fraction
    contains
      procedure :: zero
      procedure :: load
      procedure :: is_generative ! or is_mortality()?????
      procedure :: dump
  end type vm_event
  
  ! Used to store events ingested from the the prescribed event driver file:
  ! May hold events from other sources in the future as well.
  type(vm_event), private :: vm_generative_event, vm_mortality_event

  ! VM event codes:
  integer, parameter, private :: vm_event_null = 0 ! No event
  integer, parameter, private :: vm_event_plant = 1
  integer, parameter, private :: vm_event_thin_test1 = 2 ! In progress, currently wraps thin_row_low() to do perfect low thinning.
  integer, parameter, private :: vm_event_thin_proportional = 3
  integer, parameter, private :: vm_event_thin_low_perfect = 4
  integer, parameter, private :: vm_event_thin_low_probabilistic = 5
  integer, parameter, private :: vm_event_clearcut = 6
  !integer, parameter, private :: vm_event_XXXXX = X

  integer, parameter, private :: vm_event_generative_max = vm_event_plant
  integer, parameter, private :: vm_event_mortality_max = vm_event_clearcut
  
  ! These values are used to indicate an optional parameter was not provided by the driver call:
  ! Can we safely use only vm_empty_integer?
  integer, parameter, private :: vm_empty_integer = -99 ! Formerly -1
  real, parameter, private :: vm_empty_real = -99.0_r8 ! Formerly -1.0_r8
  !real, parameter, private :: vm_empty_array = vm_empty_integer ! Temporary until vm_event%pfts is a real array.

  !=================================================================================================
  
contains
  
  !=================================================================================================
  ! Main module event loop entry points:
  !
  ! Event loop entry points:
  ! ed_ecosystem_dynamics() [EDMainMod]
  ! |--> vegetation_management_init()
  ! |
  ! |-->  EDPatchDynamicsMod: disturbance_rates()
  ! |     |--> managed_mortality()
  ! |
  ! |-->  ed_integrate_state_variables()
  ! |     |--> EDMortalityFunctionsMod: Mortality_Derivative()
  ! |          |--> anthro_mortality_rate()
  ! |
  ! |--> managed_fecundity()
  ! |--> EDPatchDynamicsMod: spawn_patches()
  !      |--> management_fluxes()
  !      |--> spawn_anthro_disturbed_cohorts()
  !
  ! |--> [Destructor to be added!!!!!]
  !=================================================================================================

! Template:
!   subroutine kill()
!     ! ----------------------------------------------------------------------------------------------
!     ! 
!     ! ----------------------------------------------------------------------------------------------
!     
!     ! Uses:
!     
!     ! Arguments:
!     
!     ! Locals:
!     
!     ! ----------------------------------------------------------------------------------------------
!     
!   end subroutine kill
!
!   !=================================================================================================

  subroutine vegetation_management_init(is_master_processor, site) ! REVIEW!
    ! ----------------------------------------------------------------------------------------------
    ! Perform initialization of module globals and determine what vegetation management is due.
    !
    ! Determining what vegetation management is due at the current time step requires integrating
    ! several sources.  These may include:
    ! Implemented:
    ! - Logging module event codes. [IsItLoggingTime()]
    ! - Activities scheduled via the prescribed event driver input file.
    ! Not Implemented:
    ! - Harvest demand from HLM land use input streams. [hlm_use_lu_harvest == itrue]
    !   Currently only area based harvest is implemented in the main branch:
    !     hlm_harvest_units == hlm_harvest_area_fraction
    !   I have started to implement carbon demand based code but it is not yet linked to the HLM
    ! demand.
    ! - Activities that are due, based on patch flag or regime heuristic.
    !
    ! This routine will likely need a paired destructor!!!!!
    ! ----------------------------------------------------------------------------------------------
    !
    ! Uses:
    use FatesInterfaceTypesMod, only : numpft
    use EDLoggingMortalityMod, only : IsItLoggingTime, logging_time
    
    ! Arguments:
    integer, intent(in) :: is_master_processor ! Too long?
    type(ed_site_type), intent(inout), target :: site
    
    ! Locals: NA
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'vegetation_management_init() entering.'
    
    ! Initialize globals:
    
    ! The following should work but I need to add a destructor if I'm going to use it.
    
    ! The woody PFTs can be determined from parameter file flags:
    ! num_woody = count(EDPftvarcon_inst%woody == itrue)
    !allocate(woody_pfts(num_woody))
!     allocate(woody_pfts(count(EDPftvarcon_inst%woody == itrue)))
!     woody_pfts = pack((/(I, I = 1, numpft)/), (EDPftvarcon_inst%woody(pfts) == itrue))
    
    ! Currently there are no other flags so we define tree_pfts explicitly above.
    
    ! Check if a traditional logging module event is due and initialize IsItLoggingTime:
    call IsItLoggingTime(is_master_processor, site)

    ! Check if there are any vegetation management activities specified and load to globals:
    call load_prescribed_events(site)
    
    ! Make sure that traditional logging events do not co-occur with VM harvest events.
    if (logging_time .and. vm_mortality_event%code /= vm_event_null) then
      write(fates_log(),*) 'Traditional logging events can not currently co-occur in the same time step as vegetation management events that induce mortality.' ! Long message!
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    if (debug) write(fates_log(), *) 'vegetation_management_init() exiting.'
  end subroutine vegetation_management_init

  !=================================================================================================

  subroutine managed_mortality(site_in, bc_in, frac_site_primary)
    ! ----------------------------------------------------------------------------------------------
    ! Perform any mortality generating vegetation management that is scheduled or required at the
    ! current time step and calculate the resulting disturbance.
    !
    ! This routine:
    ! Determines which management activities are due for the site from several possible sources.
    !   (Existing: IsItLoggingTime() & LoggingMortality_frac()
    !
    ! Executes management activities that are due, storing the resulting mortalities in the cohorts.
    ! This will result in an update of:
    ! - The site level harvest_carbon_flux value member.  [Not sure why is this calculated here! Not fully implemented?]
    ! - The patch level disturbance_rates(dtype_ilog) member.
    ! - The cohort level lmort_direct, lmort_collateral, lmort_infra, and l_degrad members.
    ! - The cohort level vegetation management members. [more!!!!!]
    !   Note: Mortalities will not be executed until disturbance calculations are completed
    ! subsequently.  Mortalities staged here may or may not subsequently occur depending n whether
    ! they are the dominant disturbance at a patch level. (See ?????)
    !
    ! Mortality and disturbance calculations...
    !
    ! This routine is called from EDPatchDynamicsMod.F90: disturbance_rates() in the FATES event
    ! sequence.  Some of the traditional logging behavior in this routine is based on code extracted
    ! from that routine (noted below).  The behavior has been expanded and to include other
    ! vegetation management functionality.
    !
    ! This routine replaces the functionality anthro_disturbance_rate(), which was a direct
    ! restructuring of the orignal EDPatchDynamicsMod.F90: disturbance_rates() logging code.
    !
    ! ToDo:----------
    ! [Integrate multiple management notes from anthro_disturbance_rate() here or elsewhere.]
    ! Why do we need to estimate the harvest amount here?
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use FatesInterfaceTypesMod, only : hlm_current_year, hlm_current_month, hlm_current_day ! Temp!
    
    ! From anthro_disturbance_rate():
    use EDLoggingMortalityMod, only : get_harvest_rate_area
    use EDLoggingMortalityMod, only : logging_time
    use EDLoggingMortalityMod, only : LoggingMortality_frac
    use EDTypesMod, only : dtype_ilog
    use EDTypesMod, only : dump_patch, dump_cohort
    use FatesConstantsMod, only : fates_tiny
    use FatesConstantsMod, only : nearzero
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: site_in
    type(bc_in_type), intent(in) :: bc_in
    ! I would like to get this with get_frac_site_primary() but that creates a circular dependency:
    real(r8), intent(in) :: frac_site_primary
    
    ! Locals:
    logical :: thinning_needed, harvest_needed, control_needed
    
    type(ed_patch_type), pointer :: current_patch
    type(ed_cohort_type), pointer :: current_cohort
    
    ! Traditional logging module rates:
    real(r8) :: lmort_direct
    real(r8) :: lmort_collateral
    real(r8) :: lmort_infra
    real(r8) :: l_degrad         ! fraction of trees that are not killed but suffer from forest 
                                 ! degradation (i.e. they are moved to newly-anthro-disturbed 
                                 ! secondary forest patch)
    real(r8) :: dist_rate_ldist_notharvested
    real(r8) :: harvest_rate
    
    real(r8) :: patch_disturbance ! Accumulator
    real(r8) :: cohort_disturbance ! Accumulator
    real(r8) :: cohort_mort ! Mortality rate at the cohort level.
    real(r8) :: patch_mort_d ! The patchwide disturbed area fraction resulting directly from mortality.
    
    real(r8) :: c_1st, c_2nd ! Temporary: Used for temporary harvest implementation.
    
    integer(i4) :: pft_int_temp(1) ! Temporary
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'managed_mortality() entering.'
    
    ! An estimate of the harvest flux will be made and stored subsequently:
    site_in%harvest_carbon_flux = 0.0_r8
    
    ! ----------------------------------------------------------------------------------------------
    ! Note: This section has been rendered (mostly?) unnecessary beyond testing.  Most functionality
    ! has been moved into vegetation_management_init().  Review and remove.
    ! ----------------------------------------------------------------------------------------------
    
    thinning_needed = .false.
    harvest_needed = .false.
    control_needed = .false.
    
    ! Manually trigger events for initial testing: TEMPORARY!
    if (hlm_current_year == 2050 .and. hlm_current_month == 1 .and. hlm_current_day == 1) then
      thinning_needed = .true.
    else if (hlm_current_year == 1942 .and. hlm_current_month == 1 .and. hlm_current_day == 1) then
      harvest_needed = .true.
    else if (hlm_current_year == 2050 .and. hlm_current_month == 1 .and. hlm_current_day == 1) then
      control_needed = .true.
    endif
    
    ! ----------------------------------------------------------------------------------------------
    ! Calculate the impact of mortality inducing vegetation management on mortality rates and store
    ! the results in the cohort data structures.
    ! Note: A class with methods for the following would be a great way to implement the regimes.
    
    ! Activities occur in a specific order:
    ! Some management activities happen elsewhere like planting (fertilization).
    
    ! ----------------------------------------------------------------------------------------------
    ! Logging Module Event:
    !   If a traditional logging module event code has occurred we honor it.
    ! While this is important for backward compatibility several assumptions in the traditional
    ! logging module are not maintained in the other vegation management activities, so care should
    ! be taken when using them.  It seems likely isolated events or events that occur with low
    ! frequently may play well with other vegetation management behaviors.  Testing is needed.
    !
    ! The best postion for this step once multiple events are allowed is unclear.  However, it
    ! seems like placing this first increases the chance of the event occurring as expected.  The
    ! downside is it adds another difference between logging module harvests and other harvests.
    ! 
    ! This segment is based on code extracted from EDPatchDynamicsMod.F90: disturbance_rates().
    ! ----------------------------------------------------------------------------------------------
    if (logging_time) then
      if (debug) write(fates_log(), *) 'Logging module event beginning.'
      
      current_patch => site_in%oldest_patch
      do while (associated(current_patch))
        current_cohort => current_patch%shortest
        
        do while(associated(current_cohort))
          current_cohort%patchptr => current_patch ! Why is this necessary?????
          
          call LoggingMortality_frac(current_cohort%pft, current_cohort%dbh, &
                                     current_cohort%canopy_layer, lmort_direct, lmort_collateral, &
                                     lmort_infra, l_degrad, bc_in%hlm_harvest_rates, &
                                     bc_in%hlm_harvest_catnames, bc_in%hlm_harvest_units, &
                                     current_patch%anthro_disturbance_label, &
                                     current_patch%age_since_anthro_disturbance, frac_site_primary)

          current_cohort%lmort_direct     = lmort_direct
          current_cohort%lmort_collateral = lmort_collateral
          current_cohort%lmort_infra      = lmort_infra
          current_cohort%l_degrad         = l_degrad

        ! Estimate the wood product (trunk_product_site):
        ! Note:  Preexisting logging module code replaced by cohort_harvestable_biomass().
        site_in%harvest_carbon_flux = site_in%harvest_carbon_flux + &
                                      cohort_harvestable_biomass(cohort = current_cohort, &
                                                            harvest_profile = logging_traditional, &
                                                            staged = .true.)
        
          current_cohort => current_cohort%taller
        end do ! Cohort loop.
      
        current_patch => current_patch%younger
      end do ! Patch loop.
    endif ! (logging_time)
    
    ! The following sections of code are being replaced.  The following if() hack turns them off:
    if (.false.) then ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    ! ----------------------------------------------------------------------------------------------
    ! Thinning:
    ! Perform thinning before harvest since wood for thinning can be used to reduce harvest demand.
    !
    ! As an initial test criteria search for an appropriate aged patch and thin it:
    ! Note: Thinning may send mass to the harvest pool.  thin_row_low() estimates harvest but we do
    ! not retrieve it in the code below.
    ! ----------------------------------------------------------------------------------------------
    if (thinning_needed) then
!       if (debug) write(fates_log(), *) 'VM thinning event beginning.'
!       
!       ! Making the following work would be a bit of a pain since it returns an variable length array
!       ! of pointers, which is not simple in Fortran:
!       ! thin_patches = find_thinnable_patches(site_in, criteria)
!       ! For now either just return the first or use a comparator function:
!       
!       current_patch => site_in%oldest_patch
!       do while (associated(current_patch))
!         if (thinnable_patch(patch = current_patch, pfts = tree_pfts, goal_basal_area = 20.0_r8)) then
!           call thin_row_low(patch = current_patch, pfts = woody_pfts, &
!                             row_fraction = 0.2_r8, final_basal_area = 25.0_r8)
!         endif
!         current_patch => current_patch%younger
!       end do
      
      ! Estimate the woodproduct (trunk_product_site) if not done already. !!!!!
!       cohort_harvestable_biomass(cohort = current_cohort, harvest_profile = bole_harvest, &
!                                  staged = .true.)
      
    endif
    
    ! ----------------------------------------------------------------------------------------------
    ! Harvest:
    !
    ! For testing do prioritized harvest across both primary and secondary patches:
    ! ----------------------------------------------------------------------------------------------
    if (harvest_needed) then
      if (debug) write(fates_log(), *) 'VM harvest event beginning.'
      
      ! Initial test:
      !call kill(cohort = site_in%oldest_patch%tallest, kill_fraction = 0.5_r8)
      ! Lets make sure it is noticeable by hitting all the cohorts of 2 PFTs:
!       current_patch => site_in%oldest_patch
!       do while (associated(current_patch))
!         current_cohort => current_patch%shortest
!         do while(associated(current_cohort))
!           !if (current_cohort%pft == 4 .or. current_cohort%pft == 12) then
!           if (any(current_cohort%pft == [2,4,12])) then
!             call kill(cohort = current_cohort, flux_profile = bole_harvest, kill_fraction = 0.5_r8, &
!                       area_fraction = 1.0_r8) ! Leaving out the area_fraction right now won't work.  Fix that.
!           endif
!           current_cohort => current_cohort%taller
!         end do ! Cohort loop.
!         current_patch => current_patch%younger
!       end do ! Patch loop.
      
      ! Test 2:
      c_1st = 10000.0_r8 ! 2000.0_r8
      c_2nd = 0.0_r8
      call harvest_mass_min_area(site_in = site_in, harvest_c_primary = c_1st, harvest_c_secondary = c_2nd, & ! REVIEW!
                                 pfts = tree_pfts, dbh_min = 10.0_r8)
      
      ! Report the amount harvested to the log:
      write(fates_log(), *) 'Amount harvested: ', c_1st
      
      ! Estimate the woodproduct (trunk_product_site) if not done already. !!!!!
      
    endif
    
    ! ----------------------------------------------------------------------------------------------
    ! Understory / Competition Control:
    !
    ! Understory control can probably go anywhere but we put it after harvest so we have the option
    ! to do site prep in the same time step. (One multiple events are enabled.)
    ! ----------------------------------------------------------------------------------------------
    if (control_needed) then
      if (debug) write(fates_log(), *) 'VM understory control event beginning.'
      
      !postharvest_control()
      ! Find the patch that was most recently cleared and perform understory control on it:
      
      !current_patch => site_in%youngest_patch ! Start with the youngest patch:
!       do while (associated(current_patch) & patch_is_bare(current_patch) /= .true.)
!         current_patch => current_patch%older
!       end do     NEED TO ADD patch_is_bare()
      
!       if (associated(current_patch)) then
!         call understory_control(current_patch, method_mow)
!       else
!         ! A bare patch was not found.  Note it in the log and proceed:
!         write(fates_log(),*) 'anthro_disturbance_rate_2()?????: No bare patch found.'
!       endif
      
      current_patch => site_in%oldest_patch
      do while (associated(current_patch))
        call understory_control(current_patch, method_mow)
        current_patch => current_patch%younger
      end do ! Patch loop.
      
    endif
    
    ! End 'comment' if().
    endif ! (.false.) ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    ! ----------------------------------------------------------------------------------------------
    ! Vegetation Management Events (Mortality Inducing):
    ! If a VM event has been prescribed by the driver input file execute it:
    !
    ! Note:  I haven't figured out how heuristic / land use time series events will work yet so this
    ! logic will likely be revised again.
    !
    ! Currently we have only implemented a single event type that is used to call thin_row_low().
    ! This is more less a place holder.  More robust set of behaviors will be added in near future.
    ! ----------------------------------------------------------------------------------------------
    
    if (debug) then
      write(fates_log(),*) 'vm_generative_event:'
      call vm_generative_event%dump()
      write(fates_log(),*) 'vm_mortality_event:'
      call vm_mortality_event%dump()
    end if
    
    if (vm_mortality_event%code /= vm_event_null) then
      if (debug) write(fates_log(), *) 'VM mortality event initiating.'
      
      ! I should probably just use this as the master conditional and do nothing for vm_event_null.
      ! That may not work so well if multiple events are allowed in the future however.
      select case (vm_mortality_event%code)
        !case (vm_event_null)
        case (vm_event_thin_test1)
          ! For testing purposes replicate the manually triggered thinning (all patches) from above.
          
          ! Explanation: thin_row_low() takes a number of parameters, some optional.  We will limit
          ! the options here.
          ! Not relevant to driver file: patch & harvest_estimate
          ! Virtual call:
          !thin_event(pfts, row_fraction, [patch_fraction], final_basal_area, [final_stem_density])?
          
          current_patch => site_in%oldest_patch
          do while (associated(current_patch))
            if (thinnable_patch(patch = current_patch, &
                pfts = vm_mortality_event%pfts, &
                goal_basal_area = vm_mortality_event%final_basal_area)) then
              
              call thin_row_low(patch = current_patch, &
                                pfts = vm_mortality_event%pfts, &
                                row_fraction = vm_mortality_event%row_fraction, &
                                final_basal_area = vm_mortality_event%final_basal_area)
            
              ! Add accumulation of harvest here????
            endif
            current_patch => current_patch%younger
          end do
          
        case (vm_event_thin_proportional)
          call thin_proportional(site = site_in, pfts = vm_mortality_event%pfts, &
                                 thin_fraction = vm_mortality_event%thin_fraction)
          
        case (vm_event_thin_low_perfect)
          call thin_low_perfect(site = site_in, pfts = vm_mortality_event%pfts, &
                                thin_fraction = vm_mortality_event%thin_fraction)
          
        case (vm_event_thin_low_probabilistic)
          call thin_low_probabilistic(site = site_in, pfts = vm_mortality_event%pfts, &
                                      thin_fraction = vm_mortality_event%thin_fraction)
          
        case default
          write(fates_log(),*) 'Unrecognized event code:', vm_mortality_event%code
          call endrun(msg = errMsg(__FILE__, __LINE__))
          ! Add checking for generative events?
      end select ! (vm_mortality_event%code)
      
    end if ! (vm_mortality_event%code /= vm_event_null)
    
    ! ----------------------------------------------------------------------------------------------
    ! Planting will occur later in the event loop.
    ! ----------------------------------------------------------------------------------------------
    
    ! ----------------------------------------------------------------------------------------------
    ! Recalculate total canopy area prior to resolving the disturbance:
    !
    ! JMR Note:
    ! I think this is the first time that the canopy area has been updated since growth was applied.
    ! It is needed for logging disturbance calculation but not for other disturbance types so I
    ! removed it from EDPatchDynamicsMod.F90: disturbance_rates().  It needs to happen prior
    ! to disturbance calculation and it happens following mortality calculation in the original code
    ! because cohort level crown area is calculated in the same loop.  We may be changing that?????
    ! ----------------------------------------------------------------------------------------------
    
    current_patch => site_in%oldest_patch
    do while (associated(current_patch))
      current_patch%total_canopy_area = 0.0_r8
      current_cohort => current_patch%shortest
      do while(associated(current_cohort))
        if (current_cohort%canopy_layer == 1) then
          current_patch%total_canopy_area = current_patch%total_canopy_area + current_cohort%c_area
        endif
        current_cohort => current_cohort%taller
      end do ! Cohort loop.
      current_patch => current_patch%younger
    end do ! Patch loop.
    
    ! ----------------------------------------------------------------------------------------------
    ! Calculate the patch level disturbance rates based on the cumulative effect of management
    ! mortality rates:
    ! 
    ! The resulting rate is stored in current_patch%disturbance_rates(dtype_ilog).  The name of the
    ! dtype_ilog index value should be reconsidered and whether more than one class is needed.
    !
    ! Note: Not all mortality caused by management is disturbance generating?????  A challenge is to
    ! have as much of this happen as possible even when anthro disturbance is not dominant.
    ! [More or move!!!!!]
    !
    ! This traditional logging module method code only consideres the mortality in the top layer of
    ! the canopy as disturbing mortality.  We don't make that distinction.  I'm not sure of the
    ! the consequences of this, although the result does effect the result of
    ! anthro_mortality_rate().  This may all about not yielding a disturbance rate over 1 based on
    ! canopy area.  I added reporting in the new code.
    ! ----------------------------------------------------------------------------------------------
    
    if (logging_time) then
      ! Traditional logging module events:
      ! This segment is based on code extracted from EDPatchDynamicsMod.F90: disturbance_rates().
      !
      ! The complexity of the following code is deceptive.  The end result is that the area logged
      ! is always completely disturbed.  LoggingMortality_frac() calculates direct and indirect
      ! mortality rates and labels any remaining canopy area as degraded.  The following code
      ! adds any tree free ground.  The result is that %disturbance_rates(dtype_ilog) = 1 (or very
      ! close to that) when hlm_use_lu_harvest is false and =  hlm_harvest_rates when
      ! hlm_use_lu_harvest is true. This begs the question if this code should be simplified to
      ! reflect this.
      
      current_patch => site_in%oldest_patch
      do while (associated(current_patch))   
      
        current_patch%disturbance_rates(dtype_ilog)  = 0.0_r8
        dist_rate_ldist_notharvested = 0.0_r8
        
        current_cohort => current_patch%shortest
        do while(associated(current_cohort))
          if (current_cohort%canopy_layer == 1) then
            ! Logging Disturbance Rate:
            current_patch%disturbance_rates(dtype_ilog) = &
                                current_patch%disturbance_rates(dtype_ilog) + &
                                min(1.0_r8, current_cohort%lmort_direct + &
                                current_cohort%lmort_collateral + &
                                current_cohort%lmort_infra + &
                                current_cohort%l_degrad) * &
                                current_cohort%c_area/current_patch%area
            
            ! Non-harvested part of the logging disturbance rate:
            dist_rate_ldist_notharvested = dist_rate_ldist_notharvested + current_cohort%l_degrad * &
                                           current_cohort%c_area / current_patch%area
          endif ! (current_cohort%canopy_layer == 1)
        
          current_cohort => current_cohort%taller
        enddo ! Cohort loop.
        
        ! For non-closed-canopy areas subject to logging, add an additional increment of area
        ! disturbed equivalent to the fraction logged to account for transfer of interstitial ground
        ! area to new secondary lands
        if ((current_patch%area - current_patch%total_canopy_area) .gt. fates_tiny) then

          call get_harvest_rate_area(current_patch%anthro_disturbance_label, &
                                     bc_in%hlm_harvest_catnames, bc_in%hlm_harvest_rates, &
                                     frac_site_primary, current_patch%age_since_anthro_disturbance, &
                                     harvest_rate)

          current_patch%disturbance_rates(dtype_ilog) = current_patch%disturbance_rates(dtype_ilog) + &
              (current_patch%area - current_patch%total_canopy_area) * harvest_rate / current_patch%area

          ! Non-harvested part of the logging disturbance rate, i.e. the area to treeless area to
          ! transfer to the newly disturbed patch along with the timbered portion.
          dist_rate_ldist_notharvested = dist_rate_ldist_notharvested + &
                                         (current_patch%area - current_patch%total_canopy_area) * &
                                         harvest_rate / current_patch%area
        endif ! (current_patch%area ...

        ! Fraction of the logging disturbance rate that is non-harvested:
        if (current_patch%disturbance_rates(dtype_ilog) .gt. nearzero) then
            current_patch%fract_ldist_not_harvested = dist_rate_ldist_notharvested / &
                                                     current_patch%disturbance_rates(dtype_ilog)
        endif

        current_patch => current_patch%younger
      enddo ! Patch loop
      
    else
      ! Note: This block is now being executed whenever logging_time is false, not when other
      ! vegetation management is occurring.  This is almost certainly wrong.  It doesn't seem to
      ! be breaking anything.  Consider something like:
      ! else if (thinning_needed .or. harvest_needed .or. control_needed) then
      
      ! Other vegetation management activities:
      !   When mortality for other vegetation management activities is calculated the two values are
      ! recorded in a cohort, a mortality fraction specifying the fraction of plants in the cohort
      ! that die and a patch fraction that tells how much of the cohort (cohort area = patch area)
      ! that mortality comes from.  The patch fraction is most informative for calculating the
      ! disturbance rate.
      !
      ! What is a managed disturbance?:
      !   Given that some management activities happen 'in place' is possible to consider them
      ! non-disturbing.  For example, thinning reduces the density of trees evenly everywhere and no
      ! 'new' patch is seeming needed.  When the patch fraction is 1 this may be more or less true.
      ! However, there are times when we want to thin an area smaller than an existing patch.  In
      ! these case the patch fraction is less than 1 whether or not thinning matches our internal
      ! definition of a disturbance, the event requires a patch splitting and therefore is
      ! disturbing in the FATES sense.
      !   Since disturbing a whole patch produces a new patch with the same structure as one with
      ! the same mortalities that was not split (both result in one patch) we can treat all
      ! activities that change composition through mortally as fully disturbing and try to ignore
      ! the disturbance philosophy.
      !   However, while this decision doesn't effect the patch composition it does effect the patch
      ! age.  It is a reasonable question as to whether the anthro-disturbance flag should be set
      ! for all activities and whether intermediate opperations should not reset the patch age.
      !
      ! Calculating disturbance:
      !   While the mortality fractions are essentially rates (dead plants / total plants / event,
      ! with the per event implied) and are therefore additive, patch fractions are not. To get the
      ! patch level disturbance rate we combine them in a manner discussed in detail in the Managed
      ! Mortality Primitives section.  In short, with one management mode / flux profile we expect
      ! that some cohorts may be unaffected (patch fraction = 0) while others are (0 < patch
      ! fraction <= 1).  The effected cohorts should all have the same patch fraction.  More than
      ! one value implies more that one activities, and while that is potentially resolvable we are
      ! not ready to handle that yet.  A single non-zero patch fraction(s) gives us the disturbance
      ! rate directly without an need for accumulation.
      
      ! Loop over the cohorts in each patch and find the largest patch fraction while checking
      ! values for valid combinations:
      current_patch => site_in%oldest_patch
      do while (associated(current_patch))
        
        current_patch%disturbance_rates(dtype_ilog) = 0.0_r8
        patch_disturbance = 0.0_r8
        patch_mort_d = 0.0_r8
        
        current_cohort => current_patch%shortest
        do while(associated(current_cohort))
          
          ! Checking for more than one mortality type was previous checked in ?????
          
          cohort_disturbance = max(current_cohort%vm_pfrac_in_place, current_cohort%vm_pfrac_bole_harvest)
          
          ! Make sure the disturbance is consistant across cohorts:
          if (cohort_disturbance /= 0.0_r8) then ! 
            if (patch_disturbance == 0.0_r8) then
              patch_disturbance = cohort_disturbance
              
            ! Switch to a constant!!!!! nearzero is too small.  rsnbl_math_prec maybe?:
            else if (abs(patch_disturbance - cohort_disturbance) > 1.0e-13_r8) then
              ! Don't allow multiple (different) management activities to co-occur:
              write(fates_log(),*) 'More than one disturbance rate was detected in this patch.'
              write(fates_log(),*) 'patch_disturbance  = ', patch_disturbance
              write(fates_log(),*) 'cohort_disturbance = ', cohort_disturbance
              call dump_patch(current_patch)
              call dump_cohort(current_cohort)
              call endrun(msg = errMsg(__FILE__, __LINE__))
             endif
           endif
          
          ! Calculate the total patch wide mortality rate:
          ! With more than one activity we would sum them...
          cohort_mort = max(current_cohort%vm_mort_bole_harvest, &
                            current_cohort%vm_mort_in_place) * &
                        current_cohort%c_area / current_patch%area
          
          ! Calculate the patch wide mortality disturbance by weighing the mortality of each cohort
          ! by its canopy area:
          patch_mort_d = patch_mort_d + (cohort_mort * current_cohort%c_area / current_patch%area)
          
          current_cohort => current_cohort%taller
        end do ! Cohort loop.
        
        ! Error checking:
        if (patch_mort_d > patch_disturbance) then ! Add a tolerance and adjust?????
          write(fates_log(),*) 'Patch level mortality is > patch disturbance:'
          write(fates_log(),*) 'Patch level mortality = ', & patch_mort_d
          write(fates_log(),*) 'Patch disturbance = ', & patch_disturbance
          call endrun(msg = errMsg(__FILE__, __LINE__))
        else if (patch_mort_d > 1.0_r8 + nearzero) then
          write(fates_log(),*) 'Patch level mortality is > 1, = ', & patch_mort_d
          call endrun(msg = errMsg(__FILE__, __LINE__))
        else if (patch_mort_d > 1.0_r8) then
          patch_mort_d = 1.0_r8 ! Correct for floating point math.
        endif
        ! Move patch_disturbance error checking up?
        
        ! Store the accumulated disturbance:
        current_patch%disturbance_rates(dtype_ilog) = patch_disturbance
        ! This may be possible since we aren't filtering by canopy layer and would probably be bad.
        if (patch_disturbance > 1.0_r8 + nearzero) then
          write(fates_log(),*) 'current_patch%disturbance_rates(dtype_ilog) > 1, = ', &
                                current_patch%disturbance_rates(dtype_ilog)
        endif
        
        ! Calculate the unharvested fraction:
        ! patch_mort_d is the mortality area at the whole patch level.  Subtracting that from the
        ! disturbance fraction gives us the unharmed fraction of the original patch that will be
        ! transfered (potentially) to a new patch.  We convert that to a fraction of the new patch.
        ! Should this be recorded in all cases or just for harvest types?????
        if (current_patch%disturbance_rates(dtype_ilog) > nearzero) then
          current_patch%fract_ldist_not_harvested = (patch_disturbance - patch_mort_d) / &
                                                   patch_disturbance
          
          ! Checking the resulting value:
          if (current_patch%fract_ldist_not_harvested > 1.0_r8 + nearzero) then
          ! %fract_ldist_not_harvested is not currently in dump_patch():
            write(fates_log(),*) 'patch%fract_ldist_not_harvested is > 1, = ', &
                                  current_patch%fract_ldist_not_harvested
            call endrun(msg = errMsg(__FILE__, __LINE__))
          else if (current_patch%fract_ldist_not_harvested > 1.0_r8) then
            current_patch%fract_ldist_not_harvested = 1.0_r8
          endif ! (current_patch%fract_ldist_not_harvested > 1.0_r8 + nearzero)
          
        else
          current_patch%fract_ldist_not_harvested = 0.0_r8
        endif ! (current_patch%disturbance_rates(dtype_ilog) > nearzero)
        
        if (debug) then ! Temp reporting:
          write(fates_log(), *) 'managed_mortality(): %disturbance_rates(dtype_ilog) = ', &
                                 current_patch%disturbance_rates(dtype_ilog)
          write(fates_log(), *) 'managed_mortality(): %fract_ldist_not_harvested) = ', &
                                 current_patch%fract_ldist_not_harvested
         call dump_patch(current_patch)
        end if
        
        current_patch => current_patch%younger
      end do ! Patch loop.
      
    endif ! (logging_time)
    
    if (debug) write(fates_log(), *) 'managed_mortality() exiting.'
  end subroutine managed_mortality

  !=================================================================================================

  function anthro_mortality_rate(cohort, bc_in, frac_site_primary) result(dndt_logging)
    ! ----------------------------------------------------------------------------------------------
    ! Calculate mortality resulting from human vegetation management at the cohort level.
    ! Mortality is returned as a change in number (density?) per unit time (year).
    ! This routine is called from EDMortalityFunctionsMod: Mortality_Derivative().
    !
    ! This subroutine is designed to encapsulate the management specific logic so the calling code
    ! does not have to be aware of it.
    ! This code currently extracts the logging specific code from EDMortalityFunctionsMod.F90:
    ! Mortality_Derivative() with only formating and name changes.
    !
    ! Human induced mortality includes logging but other actives as well such as thinning,
    ! understory clearing, agricultural harvest, etc.
    ! This function only returns non-disturbing mortality, and should perhaps be renamed.
    ! The logging module code excludes disturbance inducing mortality resulting from logging of trees in
    ! in the top canopy layer.  We have not made that distinction.  As a result in
    ! managed_mortality() we classify all harvest mortality as disturbing and none of it is returned
    ! here.  This is may be problematic philosophically and mathematically and may change.
    !
    ! I think the call to LoggingMortality_frac() here is incorrect because it was called earlier
    ! but for some patches the staged mortality may be overridden. I will consult with the community
    ! to see if they agree.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use EDLoggingMortalityMod, only : LoggingMortality_frac
    use FatesInterfaceTypesMod, only : hlm_freq_day
    
    ! Arguments:
    ! Similar functions take the site as well but we can get that from the cohort.
    type(ed_cohort_type), intent(inout), target :: cohort
    type(bc_in_type), intent(in) :: bc_in
    ! Would like to get this with get_frac_site_primary() but that creates a circular dependency:
    real(r8), intent(in) :: frac_site_primary
    
    ! Locals:
    real(r8) :: dndt_logging      ! Mortality rate (per day) associated with the a logging event
    ! Name will change to dndt_managment!
    integer  :: ipft              ! local copy of the pft index
    
    ! ----------------------------------------------------------------------------------------------
    
    ipft = cohort%pft
    
    call LoggingMortality_frac(ipft, cohort%dbh, cohort%canopy_layer, &
                               cohort%lmort_direct, &
                               cohort%lmort_collateral, &
                               cohort%lmort_infra, &
                               cohort%l_degrad, &
                               bc_in%hlm_harvest_rates, &
                               bc_in%hlm_harvest_catnames, &
                               bc_in%hlm_harvest_units, &
                               cohort%patchptr%anthro_disturbance_label, &
                               cohort%patchptr%age_since_anthro_disturbance, &
                               frac_site_primary)
    
    ! The logging mortality rate are expressed as a fractional rate of the cohort / event.  They
    ! need to be converted to a rate of fraction per year similar to the natural disturbance rates
    ! using hlm_freq_day.
    
    if (cohort%canopy_layer > 1)then 
       ! Include understory logging mortality rates not associated with disturbance:
       dndt_logging = (cohort%lmort_direct + cohort%lmort_collateral +  cohort%lmort_infra) / &
                       hlm_freq_day
                       ! Consider adding vegetation management moralities here!!!!!
    else
       ! Mortality from logging in the canopy is ONLY disturbance generating, don't
       ! update number densities via non-disturbance inducing death
       dndt_logging = 0.0_r8
    endif
    
  end function anthro_mortality_rate

  !=================================================================================================

  subroutine managed_fecundity(site, bc_in)
    ! ----------------------------------------------------------------------------------------------
    ! This is called during the recruitment phase of the main event loop and executes any human
    ! driven changes in plant fecundity in the form of planting, additions of seeds, and possibly
    ! species changes.
    !
    ! Currently this is largely a hack for testing planting but over time logic will be added to
    ! determine what actions are due based on event codes, input files, or heuristics.
    !
    ! ----------------------------------------------------------------------------------------------
    
    use EDLoggingMortalityMod, only : logging_time ! Temporary trigger for planting
    use FatesGlobals, only : endrun => fates_endrun
    use FatesInterfaceTypesMod, only : hlm_current_year, hlm_current_month, hlm_current_day ! For testing.
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: site
    type(bc_in_type), intent(in) :: bc_in
    
    ! Locals:
    type(ed_patch_type), pointer :: thisPatch
    
    ! ----------------------------------------------------------------------------------------------
    
    if (debug) write(fates_log(), *) 'managed_fecundity() entering.'
    
    ! ToDo: 
    ! Determine which patch or patches to plant into:
    ! Considerations include the area that needs to be planted, the composition of the patches: i.e.
    ! what species are there, is it a bare patch, and history: is it secondary or managed land.
    !
    ! For initial testing the patch doesn't matter much but we want to make sure it is big enough
    ! that any new cohorts are not eliminated.
    
    ! For testing we use the logging event code to trigger a planting event:
    !if (logging_time) then
    ! New testing:
    if (hlm_current_year == 2050 .and. hlm_current_month == 1 .and. hlm_current_day == 1) then
      
      if (debug) then
        !write(fates_log(), *) 'Planting triggered by logging event (for testing).'
      end if
      
      ! Temporary patch determination:
      thisPatch => site%oldest_patch
      do while (associated(thisPatch))
        
        ! Assuming we are using a nominal hectare of 10,000 m^2, the 1/10th should be big enough.
        if (thisPatch%area > 1000) then
          exit
        end if
        
        thisPatch => thisPatch%younger
      enddo
      
      ! Check that the loop was not a bad idea:
      if (.not. associated(thisPatch)) then
        call endrun(msg = errMsg(__FILE__, __LINE__))
      end if
      
      ! For testing the PFT and planting conditions are arbitrary and just to demonstrate things
      ! work. For simplicity of interpretation we should set up the PFTs so that those that we will
      ! plant are not present prior to planting.  We can do this several ways:
      ! If we do not include the PFTs in the parameer file we can't add them later.
      ! If we set the initial density to 0 we have to supply a density but we can't  check the
      ! behavior when the density is omitted.  This is fine for our initial tests.
      ! Using inventory initialization is probably the best initialization to test all features.
      
      ! For the first test I plant to initialize with a grass only and then plant two different tree
      ! PFTs at two different later dates.
      
      ! To test more than one use of plant in the same run we use the date:
      ! For the first test we will start in 2001 in brazil and run for at least 15 years.
!       if (hlm_current_year < 2007) then
!         ! Plant broadleaf_evergreen_tropical_tree with default settings passed in explicitly:
!         call plant(site = site, patch = thisPatch, bc_in = bc_in, pft_index = 1, density = 0.2_r8, &
!                    height = 1.3_r8)
!       else
!         ! Plant broadleaf_hydrodecid_tropical_tree without a size specified:
!         ! We have to provide the density because we have changed the initial density to 0 to
!         ! suppress trees growth initially.
!         call plant(site = site, patch = thisPatch, bc_in = bc_in, pft_index = 5, density = 0.2_r8)
!         
!       end if
      ! Test plant() with all optional arguments omitted. Plant needleleaf_evergreen_extratrop_tree:
      ! call plant(site = site, patch = thisPatch, bc_in = bc_in, pft_index = 2)
      
      ! No planting for now.
        
    end if ! if (logging_time)
    
    
    ! Execute any generative events that came in via the driver file:-------------------------------
    if (vm_generative_event%code /= vm_event_null) then
      if (debug) write(fates_log(), *) 'VM generative event initiating.'
      
      select case (vm_generative_event%code)
        !case (vm_event_null)
        case (vm_event_plant)
          call plant_site(site = site, bc_in = bc_in, pfts = vm_generative_event%pfts, &
                          density = vm_generative_event%density, dbh = vm_generative_event%dbh, &
                          height = vm_generative_event%height)
          
        case default
          write(fates_log(),*) 'Unrecognized event code:', vm_generative_event%code
          call endrun(msg = errMsg(__FILE__, __LINE__))
          ! Add checking for generative events?
      end select ! (vm_generative_event%code)
      
    end if ! (vm_generative_event%code /= vm_event_null)
    
    if (debug) write(fates_log(), *) 'managed_fecundity() exiting.'
  end subroutine managed_fecundity
  
  !=================================================================================================

  subroutine management_fluxes(current_site, current_patch, new_patch, patch_site_areadis) ! REVIEW!
    ! ----------------------------------------------------------------------------------------------
    ! Calculate the fluxes resulting from all the management activities performed during this
    ! timestep to litter, soil, product pools, the atmosphere, and partition stocks between the
    ! existing and new patches.
    !
    ! This code is based on EDLoggingMortalityMod: logging_litter_fluxes() with some improvements
    ! from and EDPatchDynamicsMod: mortality_litter_fluxes().
    ! The major change in the code is the addition of multiple management mode flux profiles.  The
    ! code that loads / calculates mortality rates has been updated and the there are some switches
    ! added for fluxes that differ between profiles.  For most of the remaining code mainly
    ! formatting, code order, and variable name changes have been made.
    !
    ! This could easily be expanded to include natural mortality and probably fire too.
    !
    ! Another, perhaps more elegant and extensible, way to handle these different flux profiles is
    ! via a class of flux profile objects.  That would allow inheritance of shared behaviors
    ! (probably starting from natural mortality) and a more modular strucutre that could be shared
    ! with natural and fire mortality.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    ! use EDtypesMod, only : area
    !use EDLoggingMortalityMod, only : harvest_litter_localization
    use EDParamsMod, only : logging_coll_under_frac
    use EDPftvarcon, only : GetDecompyFrac
    use EDtypesMod,   only : ed_site_type
    use EDtypesMod,   only : ed_patch_type
    use EDtypesMod,   only : ed_cohort_type
    use FatesAllometryMod, only : carea_allom
    use FatesAllometryMod, only : set_root_fraction
    use FatesConstantsMod, only : rsnbl_math_prec
    use FatesInterfaceTypesMod, only : hlm_use_planthydro
    use FatesLitterMod, only : ncwd ! The number of coarse woody debris pools.
    use FatesLitterMod, only : ndcmpy
    use FatesLitterMod, only : litter_type
    use FatesPlantHydraulicsMod, only : AccumulateMortalityWaterStorage
    use SFParamsMod, only : SF_val_cwd_frac
    use EDParamsMod, only : logging_export_frac
    use EDTypesMod, only : site_massbal_type
    use EDTypesMod, only : site_fluxdiags_type
    use PRTGenericMod, only : struct_organ, leaf_organ, fnrt_organ, sapw_organ, store_organ, repro_organ
    use PRTGenericMod, only : num_elements, element_list, carbon12_element
    
    ! Arguments:
    ! Possibly unnecessary, we could get the site from the current_patch:
    type(ed_site_type), intent(inout), target :: current_site
    type(ed_patch_type), intent(inout), target :: current_patch ! The extising / donor patch.
    type(ed_patch_type), intent(inout), target :: new_patch ! The new patch from disturbance.
    ! This can be calculated from data in the cohort so may not really need to be passed in:
    real(r8), intent(in) :: patch_site_areadis ! Total area disturbed in m2 per patch per day
    
    ! Locals:
    type(ed_cohort_type), pointer :: current_cohort
    type(site_massbal_type), pointer :: site_mass
    type(site_fluxdiags_type), pointer :: flux_diags
    type(litter_type), pointer :: new_litt
    type(litter_type), pointer :: cur_litt
    
    real(r8) :: remainder_area      ! Current patch's remaining area after donation [m2]
    real(r8) :: retain_frac         ! Fraction of mass to be retained
    real(r8) :: donate_frac         ! Fraction of mass to be donated
    real(r8) :: donate_m2           ! Area normalization for litter mass destined to new patch [m-2]
    real(r8) :: retain_m2           ! Area normalization for litter mass destined to old patch [m-2]
    
    real(r8) :: direct_dead   ! Death count directly from management.
    real(r8) :: indirect_dead ! Indirect death count: impacts, infrastructure and collateral damage
    real(r8) :: all_dead      ! Both of the above.

    integer :: flux_profile ! Vegetation management flux profile.
    logical :: harvest_mode ! Is the vegetaiton managment flux profile a harvest profile?

    integer :: el                          ! Element loop index
    integer :: element_id                  ! Parteh compatible global element index
    integer :: pft                         ! PFT index
    integer :: cwd_pool                    ! CWD index
    integer :: cwd_pool_stop               ! CWD pool to stop on
    integer, parameter :: cwd_trunk = ncwd ! Trunk CWD pool index
    integer :: soil_layer                  ! soil layer loop index
    integer :: nlevsoil                    ! number of soil layers
    integer :: dcmpy                       ! index for decomposability pools
    real(r8) :: dcmpy_frac                 ! fraction going into each decomposability pool

    real(r8) :: leaf_m              ! leaf element mass      [kg]
    real(r8) :: fnrt_m              ! fineroot element mass  [kg]
    real(r8) :: sapw_m              ! sapwood element mass   [kg]
    real(r8) :: store_m             ! storage element mass   [kg]
    real(r8) :: struct_m            ! structure element mass [kg]
    real(r8) :: repro_m             ! reproductive mass      [kg]
    real(r8) :: ag_wood             ! above ground wood mass [kg]
    real(r8) :: bg_wood             ! below ground wood mass [kg]
    real(r8) :: leaf_litter         ! Leafy biomass transferred through mortality      [kgC/site]
    real(r8) :: root_litter         ! Rooty + storage biomass transferred through mort [kgC/site]

    real(r8) :: trunk_product_site  ! Flux of carbon in trunk wood exported off site   [kgC/site] 
                                    ! Note: we accumulate over the patch, but the ultimate scale is
                                    ! of this flux is at the site level.
    real(r8) :: delta_litter_stock  ! Total litter flux carbon                         [kgC/site]
    real(r8) :: delta_biomass_stock ! Total mortality carbon flux (litter + product)   [kgC/site]
    real(r8) :: delta_individual    ! change in plant number through mortality         [plants/site]
    
    ! See EDLoggingMortalityMod about this:
    real(r8), parameter :: harvest_litter_localization = 0.0_r8
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'management_fluxes() entering.'
    
    direct_dead   = 0.0_r8
    indirect_dead = 0.0_r8
    
    ! The fraction of the patch that has been disturbed and will form the new patch is provided by
    ! patch_site_areadis.  Calculate the fraction remaining in the original donor patch:
    remainder_area = current_patch%area - patch_site_areadis
    
    ! Calculate the fraction of litter to be retained in the donor patch and passed to the new one:
    retain_frac = (1.0_r8 - harvest_litter_localization) * &
                  remainder_area / (new_patch%area + remainder_area)
    donate_frac = 1.0_r8 - retain_frac
    
    ! Convert the factions to areas:
    if(remainder_area > rsnbl_math_prec) then
      retain_m2 = retain_frac / remainder_area
      donate_m2 = (1.0_r8 - retain_frac) / new_patch%area ! donate_frac / new_patch%area
    else
       retain_m2 = 0.0_r8
       donate_m2 = 1.0_r8 / new_patch%area
    end if
    
    nlevsoil = current_site%nlevsoil
    
    ! For each element process the fluxes:
    do el = 1, num_elements
      
      element_id = element_list(el)
      site_mass => current_site%mass_balance(el)
      flux_diags=> current_site%flux_diags(el)
      cur_litt  => current_patch%litter(el) ! Litter pool of "current" patch
      new_litt  => new_patch%litter(el) ! Litter pool of "new" patch
      
      ! Zero some site level accumulator diagnostics:
      trunk_product_site  = 0.0_r8
      delta_litter_stock  = 0.0_r8
      delta_biomass_stock = 0.0_r8
      delta_individual    = 0.0_r8
      
      ! --------------------------------------------------------------------------------------------
      ! Determine and execute the fluxes for each cohort:
      ! --------------------------------------------------------------------------------------------
      current_cohort => current_patch%shortest
      do while(associated(current_cohort))
        pft = current_cohort%pft
        
        sapw_m   = current_cohort%prt%GetState(sapw_organ, element_id)
        struct_m = current_cohort%prt%GetState(struct_organ, element_id)
        leaf_m   = current_cohort%prt%GetState(leaf_organ, element_id)
        fnrt_m   = current_cohort%prt%GetState(fnrt_organ, element_id)
        store_m  = current_cohort%prt%GetState(store_organ, element_id)
        repro_m  = current_cohort%prt%GetState(repro_organ, element_id)
        
        ! ------------------------------------------------------------------------------------------
        ! Determine the flux profile corresponding to the managed mortality:
        ! ------------------------------------------------------------------------------------------
        flux_profile = get_flux_profile(current_cohort)
        
        ! ------------------------------------------------------------------------------------------
        ! Get the number of plants that have died:
        ! 
        ! Both logging_litter_fluxes() & mortality_litter_fluxes():
        ! - Assume there will be no mortality of grasses.
        ! - Only calculate mortalities for the woody plants.
        ! [- Assume only plants in the upper most canopy layer will be directly harvested.] Not so?
        ! - They also both calculate the understory mortality, in the case of logging as a function
        ! of harvest.  The understory here is everything below the top canopy layer, not by PFT.
        !
        ! For management generally this is not true.  We allow direct killing of any PFT and the
        ! mortality will be stored in the cohort.  (So don't skip any)
        ! 
        ! The traditional logging module breaks down mortality into direct and indirect classes.
        ! While the initial set of vegetation management activities do not currently result in
        ! indirect mortality they may in the future so we maintain compatibility  with this feature.
        ! ------------------------------------------------------------------------------------------
        
        if (flux_profile == null_profile) then
          ! This only means this cohort didn't experience mortality, others in the patch did.
          ! All the flux calculations below will cancel out if there is no mortality.
          direct_dead   = 0.0_r8
          indirect_dead = 0.0_r8

        else if (flux_profile == logging_traditional) then
          ! Traditional logging module functionality:-----------------------------------------------
          if (current_cohort%canopy_layer == 1) then
            direct_dead = current_cohort%n * current_cohort%lmort_direct
            indirect_dead = current_cohort%n * &
                            (current_cohort%lmort_collateral + current_cohort%lmort_infra)
          else
            ! This routine is only called during disturbance.  The litter fluxes from
            ! non-disturbance generating mortality are handled in EDPhysiology.  Disturbance
            ! generating mortality are those cohorts in the top canopy layer, or those plants that
            ! were impacted. Thus, no direct dead can occur here, and indirect = collateral impacts.
            if (int(prt_params%woody(current_cohort%pft)) == itrue) then
              direct_dead  = 0.0_r8
              indirect_dead = logging_coll_under_frac * &
                              (1.0_r8 - current_patch%fract_ldist_not_harvested) * &
                              current_cohort%n * (patch_site_areadis/current_patch%area) !kgC/site/day
            else
              ! If the cohort of interest is grass, it will not experience any mortality associated
              ! with the logging disturbance
              direct_dead   = 0.0_r8
              indirect_dead = 0.0_r8
            end if ! Woody
          end if ! Canopy layer
          
        else
          ! All extended Vegetation management profiles:--------------------------------------------
          
          ! Currently only one management activity should be applied at a time, i.e. only one should
          ! be non-zero, but in the future more than one may be.  Fractional mortalties are additive
          ! so the following is theoretically future safe:
          direct_dead = (current_cohort%vm_mort_in_place + current_cohort%vm_mort_bole_harvest) * &
                        current_cohort%n ! Check this, effective_n()?????, cohort area adjustment?
          indirect_dead = 0.0_r8 ! May be added in the future for some modes.
          
        endif ! flux_profile
        
        all_dead = direct_dead + indirect_dead
        
        ! ------------------------------------------------------------------------------------------
        ! Update water balance by removing dead plant water:
        ! We only need to do this once so we use the carbon element id, which should always be
        ! present.
        ! ------------------------------------------------------------------------------------------
        if ((element_id == carbon12_element) .and. (hlm_use_planthydro == itrue)) then
          call AccumulateMortalityWaterStorage(current_site, current_cohort, all_dead)
        end if
        
        ! ------------------------------------------------------------------------------------------
        ! Woody fluxes:
        ! In all current scenarios the non-bole woody organs, i.e. branches, course roots, etc. go
        ! to litter pools.  In planned future modes some may be harvested or possibly burned.
        ! In all foreseeable harvest modes the bole will be exported.
        !
        ! Non-woody plants have these pools, though they may be empty, so this is safe for all PFTs.
        ! [Confirm this!!!!!]
        ! ------------------------------------------------------------------------------------------
        
        ! If this is a harvesting flux profile reserve the stem for now:
        if (flux_profile == logging_traditional .or. flux_profile == bole_harvest) then
          harvest_mode = .true.
          cwd_pool_stop = cwd_trunk -1
        else
          harvest_mode = .false.
          cwd_pool_stop = ncwd
        endif
        
        ! Update the root fractions:
        call set_root_fraction(current_site%rootfrac_scr, pft, current_site%zi_soil)
        
        ! Calculate total, above and below ground, structural and sapwood, stem biomass:
        ag_wood = all_dead * (struct_m + sapw_m) * prt_params%allom_agb_frac(pft)
        bg_wood = all_dead * (struct_m + sapw_m) * (1.0_r8 - prt_params%allom_agb_frac(pft))
        
        ! Transfer (most to all) woody necromass to debris pools:
        ! This is exactly what is done in natural mortality as well.
        do cwd_pool = 1, cwd_pool_stop
          ! Transfer deadwood to aboveground CWD pools:
          new_litt%ag_cwd(cwd_pool) = new_litt%ag_cwd(cwd_pool) +&
                                      ag_wood * SF_val_CWD_frac(cwd_pool) * donate_m2
          
          cur_litt%ag_cwd(cwd_pool) = cur_litt%ag_cwd(cwd_pool) + &
                                      ag_wood * SF_val_CWD_frac(cwd_pool) * retain_m2
          
          ! Transfer deadwood to below-ground CWD soil layer pools:
          do soil_layer = 1, nlevsoil
            new_litt%bg_cwd(cwd_pool, soil_layer) = new_litt%bg_cwd(cwd_pool, soil_layer) + &
                                      bg_wood * current_site%rootfrac_scr(soil_layer) * &
                                      SF_val_CWD_frac(cwd_pool) * donate_m2
            
            cur_litt%bg_cwd(cwd_pool, soil_layer) = cur_litt%bg_cwd(cwd_pool, soil_layer) + &
                                      bg_wood * current_site%rootfrac_scr(soil_layer) * &
                                      SF_val_CWD_frac(cwd_pool) * retain_m2
          end do ! Soil layer loop.
          
          ! Diagnostics on fluxes into the aboveground (AG) and below-ground (BG) CWD pools
          flux_diags%cwd_ag_input(cwd_pool) = flux_diags%cwd_ag_input(cwd_pool) + &
                                              SF_val_CWD_frac(cwd_pool) * ag_wood
          flux_diags%cwd_bg_input(cwd_pool) = flux_diags%cwd_bg_input(cwd_pool) + &
                                              SF_val_CWD_frac(cwd_pool) * bg_wood
          
          ! Diagnostic specific to resource management code
          if( element_id .eq. carbon12_element) then
            delta_litter_stock  = delta_litter_stock + &
                                  (ag_wood + bg_wood) * SF_val_CWD_frac(cwd_pool)
          end if
          
        enddo ! CWD loop
        
        ! ------------------------------------------------------------------------------------------
        ! If this is a harvest activity take care of the stem wood which was not handled above:
        ! ------------------------------------------------------------------------------------------
        if (harvest_mode) then
          
          ! This recalculation is done in the logging module code but is currently not necessary for
          ! other management activities (but it is not harmful either):
          bg_wood = direct_dead * (struct_m + sapw_m ) * SF_val_CWD_frac(ncwd) * &
                    (1.0_r8 - prt_params%allom_agb_frac(current_cohort%pft))
          
          ! ----------------------------------------------------------------------------------------
          ! Harvest: The below-ground portion of the stem goes to soil:
          ! ----------------------------------------------------------------------------------------
          do soil_layer = 1, nlevsoil
            new_litt%bg_cwd(ncwd, soil_layer) = new_litt%bg_cwd(ncwd, soil_layer) + &
                                         bg_wood * current_site%rootfrac_scr(soil_layer) * donate_m2
            
            cur_litt%bg_cwd(ncwd, soil_layer) = cur_litt%bg_cwd(ncwd, soil_layer) + &
                                         bg_wood * current_site%rootfrac_scr(soil_layer) * retain_m2
          end do ! Soil layer loop.
          
          ! Diagnostics:
          flux_diags%cwd_bg_input(ncwd) = flux_diags%cwd_bg_input(ncwd) +  bg_wood
          
          ! ----------------------------------------------------------------------------------------
          ! Harvest: The aboveground portion of the stem, the bole, is exported (fluxed out of the
          ! site) as harvested wood product.
          ! 
          ! We retain the use of the traditional logging modules export fraction parameter.  This
          ! parameter is a simple fractional multiplier that can be used to reduce the wood exported
          ! from the site.  The balence (1 - export_frac) is sent to litter.  This can be used to
          ! approximate logging inefficiencies like damaged trees abandoned on site.
          ! It could also be used to adjust for the stump, although a more direct treatment of that
          ! might be worthwhile.
          ! 
          ! These original notes are probably referring to the fact that fluxes and wood product are
          ! available as output variables?
            ! Losses to the system as a whole, for C-balancing (kGC/site/day)
            ! Site level product, (kgC/site, accumulated over simulation)
          ! ----------------------------------------------------------------------------------------
          
          ! This recalculation is done in the logging module code but is currently not necessary for
          ! other management activities (but it is not harmful either):
          ag_wood = direct_dead * (struct_m + sapw_m) * &
                    prt_params%allom_agb_frac(current_cohort%pft) * SF_val_CWD_frac(ncwd)
          
          trunk_product_site = trunk_product_site + ag_wood * logging_export_frac
          
          ! This is for checking the total mass balance [kg/site/day]
          site_mass%wood_product = site_mass%wood_product + ag_wood * logging_export_frac
          
          ! Send any export losses to the litter
          new_litt%ag_cwd(ncwd) = new_litt%ag_cwd(ncwd) +&
                                  ag_wood * (1.0_r8 - logging_export_frac) * donate_m2
          
          cur_litt%ag_cwd(ncwd) = cur_litt%ag_cwd(ncwd) + &
                                  ag_wood *(1.0_r8 - logging_export_frac) * retain_m2
          
        endif ! if (harvest_mode)
        
        ! ------------------------------------------------------------------------------------------
        ! Indirect mortality wood fluxes:
        ! Handle litter fluxes of trunk wood resulting from infrastructure and collateral damage.
        !
        ! As noted above indirect mortality is currently only used with traditional logging module
        ! events but may be added to other activities in the future.
        ! When indirect_dead = 0 it all cancels out.  Consider adding a check to skip this for
        ! efficiency.
        ! ------------------------------------------------------------------------------------------
        
        ag_wood = indirect_dead * (struct_m + sapw_m) * &
                  prt_params%allom_agb_frac(current_cohort%pft)
        bg_wood = indirect_dead * (struct_m + sapw_m) * &
                  (1.0_r8 - prt_params%allom_agb_frac(current_cohort%pft))
        
        ! Aboveground wood:
        new_litt%ag_cwd(ncwd) = new_litt%ag_cwd(ncwd) + &
                                ag_wood * SF_val_CWD_frac(ncwd) * donate_m2

        cur_litt%ag_cwd(ncwd) = cur_litt%ag_cwd(ncwd) + &
                                ag_wood * SF_val_CWD_frac(ncwd) * retain_m2
        
        ! Below-ground wood:
        do soil_layer = 1, nlevsoil
          new_litt%bg_cwd(ncwd, soil_layer) = new_litt%bg_cwd(ncwd, soil_layer) + &
                                              bg_wood * current_site%rootfrac_scr(soil_layer) * &
                                              SF_val_CWD_frac(ncwd) * donate_m2
          
          cur_litt%bg_cwd(ncwd, soil_layer) = cur_litt%bg_cwd(ncwd, soil_layer) + &
                                              bg_wood * current_site%rootfrac_scr(soil_layer) * &
                                              SF_val_CWD_frac(ncwd) * retain_m2
        end do ! Soil layer loop.
        
        ! Diagnostics:
        flux_diags%cwd_ag_input(ncwd) = flux_diags%cwd_ag_input(ncwd) + &
                                        SF_val_CWD_frac(ncwd) * ag_wood
        
        flux_diags%cwd_bg_input(ncwd) = flux_diags%cwd_bg_input(ncwd) + &
                                        SF_val_CWD_frac(ncwd) * bg_wood
        
        if( element_id .eq. carbon12_element) then
          delta_litter_stock  = delta_litter_stock + (ag_wood + bg_wood) * SF_val_CWD_frac(ncwd)
        end if
        
        ! ------------------------------------------------------------------------------------------
        ! Move leaf, fine root, and storage carbon into debris pools:
        !
        ! Currently this is the same for the traditional logging behavior and vegetation management.
        ! Natural tree fall mortality differs in that it moves some storage carbon to the seed pool,
        ! which we do not currently do here.
        !
        ! In the future this may change for biomass or burn profiles.
        ! ------------------------------------------------------------------------------------------
        
        leaf_litter = all_dead * (leaf_m + repro_m)
        root_litter = all_dead * (fnrt_m + store_m)
        
        do dcmpy=1, ndcmpy
          
          ! Foliage and seeds into litter:
          dcmpy_frac = GetDecompyFrac(pft, leaf_organ, dcmpy)
          
          new_litt%leaf_fines(dcmpy) = new_litt%leaf_fines(dcmpy) + &
                                       leaf_litter * dcmpy_frac * donate_m2
          
          cur_litt%leaf_fines(dcmpy) = cur_litt%leaf_fines(dcmpy) + &
                                       leaf_litter * dcmpy_frac * retain_m2
          
          ! Fine roots into below-ground decomp pools:
          dcmpy_frac = GetDecompyFrac(pft, fnrt_organ, dcmpy)
          do soil_layer = 1, nlevsoil
            new_litt%root_fines(dcmpy, soil_layer) = new_litt%root_fines(dcmpy, soil_layer) + &
            root_litter * current_site%rootfrac_scr(soil_layer) * dcmpy_frac * donate_m2
            
            cur_litt%root_fines(dcmpy, soil_layer) = cur_litt%root_fines(dcmpy, soil_layer) + &
            root_litter * current_site%rootfrac_scr(soil_layer) * dcmpy_frac * retain_m2
          end do
        end do !dcmpy (decomp) loop
        
        ! Track as diagnostic fluxes:
        flux_diags%leaf_litter_input(pft) = flux_diags%leaf_litter_input(pft) + leaf_litter
        flux_diags%root_litter_input(pft) = flux_diags%root_litter_input(pft) + root_litter
        
        ! ------------------------------------------------------------------------------------------
        ! Accumulate additional diagnostics:
        ! These all come from logging_litter_fluxes().  Natural mortality does it differently.
        ! ------------------------------------------------------------------------------------------
        
        ! Note that litter stock also has terms above in the CWD loop
        if(element_id == carbon12_element) then
          delta_litter_stock  = delta_litter_stock + leaf_litter + root_litter
          delta_biomass_stock = delta_biomass_stock + leaf_litter + root_litter + &
                                all_dead * (struct_m + sapw_m)
          delta_individual = delta_individual + all_dead
        end if
        
        current_cohort => current_cohort%taller
      end do ! Cohort loop.
    
      ! --------------------------------------------------------------------------------------------
      ! Update the amount of carbon exported from the site through harvest:
      !
      ! Note: These diagnostics are straight from the traditional logging module.  That defines the
      ! harvest flux trunk_product_site to just be bole harvest.  With the additional management
      ! activities we have added the types of PFTs and canopy layer that get harvested has changed
      ! but that doe not change the interpretation.  However, soon we may harvest non-bole organs in
      ! some schemes so the there may need to some changes here and possibly new export pools.
      !
      ! Need to make sure this code stays consistent with the x_harvestable_biomass() functions!
      ! --------------------------------------------------------------------------------------------
      
      if (element_id == carbon12_element) then
        current_site%resources_management%trunk_product_site = &
                                            current_site%resources_management%trunk_product_site + &
                                            trunk_product_site
        
        current_site%resources_management%delta_litter_stock = &
                                            current_site%resources_management%delta_litter_stock + &
                                            delta_litter_stock
        
        current_site%resources_management%delta_biomass_stock = &
                                          current_site%resources_management%delta_biomass_stock + &
                                          delta_biomass_stock
        
        current_site%resources_management%delta_individual = &
                                              current_site%resources_management%delta_individual + &
                                              delta_individual
      end if ! (element_id == carbon12_element)
    end do ! num_elements loop
    
    ! ----------------------------------------------------------------------------------------------
    ! As noted in the code of logging_litter_fluxes() by RGK the following code may not be needed:
    ! ----------------------------------------------------------------------------------------------
    current_cohort => new_patch%shortest
    do while(associated(current_cohort))
      call carea_allom(current_cohort%dbh, current_cohort%n, current_site%spread, &
                       current_cohort%pft,current_cohort%c_area)
      current_cohort => current_cohort%taller
    enddo
    
    if (debug) write(fates_log(), *) 'management_fluxes() exiting.'
  end subroutine management_fluxes

  !=================================================================================================

  ! This name is rather long!!!!!
  subroutine spawn_anthro_disturbed_cohorts(parent_site, donor_cohort, new_cohort)
    ! ----------------------------------------------------------------------------------------------
    ! For cohorts that have experienced anthropogenic disturbance initialize the new disturbed
    ! cohort's number density and mortality rates and apply mortalities to adjust the numbers in
    ! the donor cohort.
    !
    ! This routine is called from EDPatchDynamicsMod: spawn_patches() and replaces code from that
    ! routine that handled logging disturbance.
    ! The former code in spawn_patches() contained a lot of logging assumptions.  The addition of
    ! more vegetation management activities with differing assumptions would make an already
    ! complicated piece of code even more so.  To reduce the need for other modules to know the
    ! internal assumptions of this module and to make the code more maintainable this routine was
    ! created.
    !
    ! This routine pairs with management_fluxes() performing complementary modifications to plant
    ! numbers.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use EDParamsMod, only : fates_mortality_disturbance_fraction
    use EDParamsMod, only : ED_val_understorey_death, logging_coll_under_frac
    use FatesInterfaceTypesMod, only : hlm_freq_day
    use FatesConstantsMod, only : days_per_sec, g_per_kg, ha_per_m2, years_per_day
    use PRTGenericMod, only : struct_organ, leaf_organ, fnrt_organ, sapw_organ, store_organ
    use PRTGenericMod, only : all_carbon_elements
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: parent_site ! The parent site of the cohorts.
    type(ed_cohort_type), intent(inout), target :: donor_cohort ! The donor cohort to copy & update.
    type(ed_cohort_type), intent(inout), target :: new_cohort ! The new cohort to initialize.
    
    ! Locals:
    integer :: flux_profile ! Vegetation management flux profile.
    real(r8) :: patch_site_areadis ! Total area disturbed in m2 per patch per day
    type(ed_patch_type), pointer :: parent_patch
    real(r8) :: plant_c ! Total plant carbon [kg]
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'spawn_anthro_disturbed_cohorts() entering.'
    
    parent_patch => donor_cohort%patchptr
    
    ! From EDPatchDynamicsMod: spawn_patches():
    ! This is the amount of patch area that is disturbed, and donated by the donor:
    patch_site_areadis = parent_patch%area * parent_patch%disturbance_rate
    ! Note: patch_site_areadis / parent_patch%area is used a lot below. Just need %patchptr%disturbance_rate?
    
    plant_c = donor_cohort%prt%GetState(sapw_organ, all_carbon_elements) + &
              donor_cohort%prt%GetState(struct_organ, all_carbon_elements) + &
              donor_cohort%prt%GetState(leaf_organ, all_carbon_elements) + &
              donor_cohort%prt%GetState(fnrt_organ, all_carbon_elements) + &
              donor_cohort%prt%GetState(store_organ, all_carbon_elements)
    ! plant_c = donor_cohort%prt%GetState(all_organs, all_carbon_elements) ! Would this work?
    
    ! Get the management activity / flux profile that is occurring in this cohort:
    flux_profile = get_flux_profile(donor_cohort)
    
    ! Temporarily short circuit the above to make sure that the historic logging code is called:
    !flux_profile = logging_traditional
    
    select case (flux_profile)
      case (logging_traditional) !------------------------------------------------------------------
        ! This logging code was exported from EDPatchDynamicsMod: spawn_patches().
        ! [Original logic intact with formatting changes.]
        
        if (debug) write(fates_log(), *) 'spawn_anthro_disturbed_cohorts() traditional logging event.'
        
        ! If this cohort is in the upper canopy. It generated 
        if (donor_cohort%canopy_layer == 1) then
           
           ! Calculate the survivorship of disturbed trees because non-harvested
           new_cohort%n = donor_cohort%n * donor_cohort%l_degrad
           ! nc%n            = (currentCohort%l_degrad / (currentCohort%l_degrad + &
           !      currentCohort%lmort_direct + currentCohort%lmort_collateral +
           !   currentCohort%lmort_infra) ) * &
           !      currentCohort%n * patch_site_areadis/currentPatch%area
           
           ! Reduce counts in the existing/donor patch according to the logging rate
           donor_cohort%n = donor_cohort%n * & (1.0_r8 - min(1.0_r8, (donor_cohort%lmort_direct + &
                            donor_cohort%lmort_collateral + & donor_cohort%lmort_infra + &
                            donor_cohort%l_degrad)))

           new_cohort%cmort            = donor_cohort%cmort
           new_cohort%hmort            = donor_cohort%hmort
           new_cohort%bmort            = donor_cohort%bmort
           new_cohort%frmort           = donor_cohort%frmort
           new_cohort%smort            = donor_cohort%smort
           new_cohort%asmort           = donor_cohort%asmort
           new_cohort%dmort            = donor_cohort%dmort

           ! Since these are the ones that weren't logged, set the logging mortality rates as zero:
           new_cohort%lmort_direct     = 0.0_r8
           new_cohort%lmort_collateral = 0.0_r8
           new_cohort%lmort_infra      = 0.0_r8
           ! JMR_MOD_START:
           ! There could be trouble above!!!!!
           new_cohort%vm_mort_in_place = 0.0_r8
           new_cohort%vm_mort_bole_harvest = 0.0_r8
           new_cohort%vm_pfrac_in_place = 0.0_r8
           new_cohort%vm_pfrac_bole_harvest = 0.0_r8
           ! JMR_MOD_END.
           
        else
           
           ! What to do with cohorts in the understory of a logging generated disturbance patch?
           
           if (int(prt_params%woody(donor_cohort%pft)) == itrue) then
              
              ! Survivorship of undestory woody plants.  Two step process.
              ! Step 1:  Reduce current number of plants to reflect the 
              !          change in area.
              !          The number density per square are doesn't change,
              !          but since the patch is smaller
              !          and cohort counts are absolute, reduce this number.
              new_cohort%n = donor_cohort%n * patch_site_areadis / parent_patch%area
              
              ! because the mortality rate due to impact for the cohorts which had 
              ! been in the understory and are now in the newly-
              ! disturbed patch is very high, passing the imort directly to 
              ! history results in large numerical errors, on account
              ! of the sharply reduced number densities.  so instead pass this info 
              ! via a site-level diagnostic variable before reducing 
              ! the number density.
              parent_site%imort_rate(donor_cohort%size_class, donor_cohort%pft) = &
                   parent_site%imort_rate(donor_cohort%size_class, donor_cohort%pft) + &
                   new_cohort%n * parent_patch%fract_ldist_not_harvested * &
                   logging_coll_under_frac / hlm_freq_day

              parent_site%imort_carbonflux = parent_site%imort_carbonflux + &
                   (new_cohort%n * parent_patch%fract_ldist_not_harvested * &
                   logging_coll_under_frac / hlm_freq_day ) * &
                   plant_c * g_per_kg * days_per_sec * years_per_day * ha_per_m2
              
              ! Step 2:  Apply survivor ship function based on the understory death fraction
              
              ! remaining of understory plants of those that are knocked 
              ! over by the overstorey trees dying...  
              ! LOGGING SURVIVORSHIP OF UNDERSTORY PLANTS IS SET AS A NEW PARAMETER 
              ! in the fatesparameter files 
              new_cohort%n = new_cohort%n * (1.0_r8 - &
                             (1.0_r8 - parent_patch%fract_ldist_not_harvested) * &
                             logging_coll_under_frac)
              
              ! Step 3: Reduce the number count of cohorts in the 
              !         original/donor/non-disturbed patch to reflect the area change
              donor_cohort%n = donor_cohort%n * (1.0_r8 - patch_site_areadis / parent_patch%area)
              
              new_cohort%cmort            = donor_cohort%cmort
              new_cohort%hmort            = donor_cohort%hmort
              new_cohort%bmort            = donor_cohort%bmort
              new_cohort%frmort           = donor_cohort%frmort
              new_cohort%smort            = donor_cohort%smort
              new_cohort%asmort           = donor_cohort%asmort
              new_cohort%dmort            = donor_cohort%dmort
              new_cohort%lmort_direct     = donor_cohort%lmort_direct
              new_cohort%lmort_collateral = donor_cohort%lmort_collateral
              new_cohort%lmort_infra      = donor_cohort%lmort_infra
              ! JMR_MOD_START:
              new_cohort%vm_mort_in_place = donor_cohort%vm_mort_in_place
              new_cohort%vm_mort_bole_harvest = donor_cohort%vm_mort_bole_harvest
              new_cohort%vm_pfrac_in_place = donor_cohort%vm_pfrac_in_place
              new_cohort%vm_pfrac_bole_harvest = donor_cohort%vm_pfrac_bole_harvest
              ! JMR_MOD_END.
              
              ! JMR_NOTE: The disturbance for the old cohort should be reset here.  See below!!!!!
              
           else
              
              ! Grass is not killed by mortality disturbance events. 
              ! Just move it into the new patch area. 
              ! Just split the grass into the existing and new patch structures.
              new_cohort%n = donor_cohort%n * patch_site_areadis / parent_patch%area
              
              ! Those remaining in the existing
              donor_cohort%n = donor_cohort%n * (1.0_r8 - patch_site_areadis / parent_patch%area)
              
              ! No grass impact mortality imposed on the newly created patch
              new_cohort%cmort            = donor_cohort%cmort
              new_cohort%hmort            = donor_cohort%hmort
              new_cohort%bmort            = donor_cohort%bmort
              new_cohort%frmort           = donor_cohort%frmort
              new_cohort%smort            = donor_cohort%smort
              new_cohort%asmort           = donor_cohort%asmort
              new_cohort%dmort            = donor_cohort%dmort
              new_cohort%lmort_direct     = donor_cohort%lmort_direct
              new_cohort%lmort_collateral = donor_cohort%lmort_collateral
              new_cohort%lmort_infra      = donor_cohort%lmort_infra
              ! JMR_MOD_START:
              new_cohort%vm_mort_in_place = donor_cohort%vm_mort_in_place
              new_cohort%vm_mort_bole_harvest = donor_cohort%vm_mort_bole_harvest
              new_cohort%vm_pfrac_in_place = donor_cohort%vm_pfrac_in_place
              new_cohort%vm_pfrac_bole_harvest = donor_cohort%vm_pfrac_bole_harvest
              ! JMR_MOD_END
              
           endif ! Is / is-not woody
        endif ! Select canopy layer

      case (in_place) !-----------------------------------------------------------------------------
        if (debug) write(fates_log(), *) 'spawn_anthro_disturbed_cohorts() in_place VM event.'
        
        ! This is almost identical to bole_harvest and should be combined!
        
        if (abs(parent_patch%disturbance_rate - donor_cohort%vm_pfrac_in_place) > 1.0e-10_r8) then
          write(fates_log(),*) 'parent_patch%disturbance_rate /= donor_cohort%vm_pfrac_in_place'
          ! call endrun(msg = errMsg(__FILE__, __LINE__))
        endif
        
        ! Update the plant numbers:
        ! Any floating point magic needed here?????
        
        ! Give the new patch the proportional number of trees for its area fraction:
        new_cohort%n = donor_cohort%n * (patch_site_areadis / parent_patch%area)
        ! Then apply all the mortality to it.
        new_cohort%n = new_cohort%n - (donor_cohort%n * donor_cohort%vm_mort_in_place)
        
        ! The old patch has no management mortality, its just smaller:
        donor_cohort%n = donor_cohort%n * (1.0_r8 - patch_site_areadis / parent_patch%area)
        
        ! Copy the mortality data members to the new cohort:
        new_cohort%cmort            = donor_cohort%cmort
        new_cohort%hmort            = donor_cohort%hmort
        new_cohort%bmort            = donor_cohort%bmort
        new_cohort%frmort           = donor_cohort%frmort
        new_cohort%smort            = donor_cohort%smort
        new_cohort%asmort           = donor_cohort%asmort
        new_cohort%dmort            = donor_cohort%dmort

        ! This follows the example of the traditional logging module events:
        ! I'm not exactly why we set the new patch to 0.  It may be that being new it has no history.
        new_cohort%lmort_direct     = 0.0_r8
        new_cohort%lmort_collateral = 0.0_r8
        new_cohort%lmort_infra      = 0.0_r8
        
        new_cohort%vm_mort_in_place      = 0.0_r8
        new_cohort%vm_mort_bole_harvest  = 0.0_r8
        new_cohort%vm_pfrac_in_place     = 0.0_r8
        new_cohort%vm_pfrac_bole_harvest = 0.0_r8
        
        ! Following the creation of the new cohort reset the disturbance values in the old cohort:
        ! Note: Repeated below!!!!!
        donor_cohort%lmort_direct     = 0.0_r8
        donor_cohort%lmort_collateral = 0.0_r8
        donor_cohort%lmort_infra      = 0.0_r8
        
        donor_cohort%vm_mort_in_place      = 0.0_r8
        donor_cohort%vm_mort_bole_harvest  = 0.0_r8
        donor_cohort%vm_pfrac_in_place     = 0.0_r8
        donor_cohort%vm_pfrac_bole_harvest = 0.0_r8
        
      case (bole_harvest) !-------------------------------------------------------------------------
        if (debug) write(fates_log(), *) 'spawn_anthro_disturbed_cohorts() bole_harvest VM event.'
        
        if (debug) then
          !write(fates_log(), *) 'Starting conditions:--------------------'
          !write(fates_log(), *) 'parent_patch%area =',  parent_patch%area
          !write(fates_log(), *) 'Patch_site_areadis =', Patch_site_areadis
          !write(fates_log(), *) 'parent_patch%disturbance_rate =', parent_patch%disturbance_rate
          !call dump_cohort(donor_cohort)
          !write(fates_log(), *) 'donor_cohort%n = ', donor_cohort%n
          !write(fates_log(), *) 'donor_cohort%lmort_direct = ', donor_cohort%lmort_direct
          !write(fates_log(), *) 'donor_cohort%vm_mort_bole_harvest = ', donor_cohort%vm_mort_bole_harvest
          !write(fates_log(), *) 'donor_cohort%vm_pfrac_bole_harvest = ', donor_cohort%vm_pfrac_bole_harvest
          ! Revised:
          write(fates_log(), *) 'Initial donor_cohort%n = ', donor_cohort%n
        end if
        
        ! Check the area:
        ! I don't know what a reasonable tolerance is.
        if (abs(parent_patch%disturbance_rate - donor_cohort%vm_pfrac_bole_harvest) > 1.0e-10_r8) then
          write(fates_log(),*) 'parent_patch%disturbance_rate /= donor_cohort%vm_pfrac_bole_harvest'
          ! call endrun(msg = errMsg(__FILE__, __LINE__))
        endif
        
        ! Update the plant numbers:
        ! Any floating point magic needed here?????
        
        ! Give the new patch the proportional number of trees for its area fraction:
        new_cohort%n = donor_cohort%n * (patch_site_areadis / parent_patch%area)
        ! Then apply all the mortality to it.
        new_cohort%n = new_cohort%n - (donor_cohort%n * donor_cohort%vm_mort_bole_harvest)
        
        ! The old patch has no management mortality, it's just smaller:
        donor_cohort%n = donor_cohort%n * (1.0_r8 - patch_site_areadis / parent_patch%area)
        
        ! Copy the mortality data members to the new cohort:
        new_cohort%cmort            = donor_cohort%cmort
        new_cohort%hmort            = donor_cohort%hmort
        new_cohort%bmort            = donor_cohort%bmort
        new_cohort%frmort           = donor_cohort%frmort
        new_cohort%smort            = donor_cohort%smort
        new_cohort%asmort           = donor_cohort%asmort
        new_cohort%dmort            = donor_cohort%dmort

        ! This follows the example of the traditional logging module events:
        ! I'm not exactly sure why we set the new patch to 0.  It may be that being new it has no history.
        new_cohort%lmort_direct     = 0.0_r8
        new_cohort%lmort_collateral = 0.0_r8
        new_cohort%lmort_infra      = 0.0_r8
        
        new_cohort%vm_mort_in_place      = 0.0_r8
        new_cohort%vm_mort_bole_harvest  = 0.0_r8
        new_cohort%vm_pfrac_in_place     = 0.0_r8
        new_cohort%vm_pfrac_bole_harvest = 0.0_r8
        
        ! Following the creation of the new cohort reset the disturbance values in the old cohort:
        ! Note: I'm not positive this is the right place to do this but failing to do so results in
        ! the cohort retaining the values into the next time step where it will be harvested again.
        donor_cohort%lmort_direct     = 0.0_r8
        donor_cohort%lmort_collateral = 0.0_r8
        donor_cohort%lmort_infra      = 0.0_r8
        
        donor_cohort%vm_mort_in_place      = 0.0_r8
        donor_cohort%vm_mort_bole_harvest  = 0.0_r8
        donor_cohort%vm_pfrac_in_place     = 0.0_r8
        donor_cohort%vm_pfrac_bole_harvest = 0.0_r8
        
        if (debug) then
          !write(fates_log(), *) 'Ending conditions:'
          !write(fates_log(), *) 'donor_cohort:'
          !call dump_cohort(donor_cohort)
          !write(fates_log(), *) 'new_cohort:'
          !call dump_cohort(new_cohort)
          !write(fates_log(), *) 'donor_cohort%n = ', donor_cohort%n, 'new_cohort%n = ', new_cohort%n
          ! Revised:
          write(fates_log(), *) 'donor_cohort%n         = ', donor_cohort%n
          write(fates_log(), *) 'new_cohort%n           = ', new_cohort%n
        end if
        
      ! case (burn)
        ! Placeholder.
      
      case (null_profile) !-------------------------------------------------------------------------
        ! This cohort has experienced no managed mortality (but others in the patch may have).
        if (debug) write(fates_log(), *) 'spawn_anthro_disturbed_cohorts() cohort has no managed mortality.'
        ! Note: The following behavior is the same as logging module grass above.
        
        ! Split the trees proportionally among the new and old patches:
        new_cohort%n = donor_cohort%n * patch_site_areadis / parent_patch%area
        donor_cohort%n = donor_cohort%n * (1.0_r8 - patch_site_areadis / parent_patch%area)
        
        ! Copy the mortality data members to the new cohort:
        new_cohort%cmort            = donor_cohort%cmort
        new_cohort%hmort            = donor_cohort%hmort
        new_cohort%bmort            = donor_cohort%bmort
        new_cohort%frmort           = donor_cohort%frmort
        new_cohort%smort            = donor_cohort%smort
        new_cohort%asmort           = donor_cohort%asmort
        new_cohort%dmort            = donor_cohort%dmort
        new_cohort%lmort_direct     = donor_cohort%lmort_direct
        new_cohort%lmort_collateral = donor_cohort%lmort_collateral
        new_cohort%lmort_infra      = donor_cohort%lmort_infra
        new_cohort%vm_mort_in_place      = donor_cohort%vm_mort_in_place
        new_cohort%vm_mort_bole_harvest  = donor_cohort%vm_mort_bole_harvest
        new_cohort%vm_pfrac_in_place     = donor_cohort%vm_pfrac_in_place
        new_cohort%vm_pfrac_bole_harvest = donor_cohort%vm_pfrac_bole_harvest
        
      case default
        write(fates_log(),*) 'Unrecognized flux profile.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
    end select
    
    if (debug) write(fates_log(), *) 'spawn_anthro_disturbed_cohorts() exiting.'
  end subroutine spawn_anthro_disturbed_cohorts

!   !=================================================================================================
!
!   subroutine vegetation_management_close()
!     ! ----------------------------------------------------------------------------------------------
!     ! 
!     ! ----------------------------------------------------------------------------------------------
!     
!     ! Uses:
!     
!     ! Arguments:
!     
!     ! Locals:
!     
!     ! ----------------------------------------------------------------------------------------------
!     
!   end subroutine kill

  !=================================================================================================
  ! Managed Fecundity (Planting) Primitives:
  !   Fecundity primitives are low level subroutines used to introduce plants and seeds to patches.
  ! They provide an interface for higher level functions to build abstractions of human vegetation
  ! management activities without depending as intimately on the the underlying implementation.
  !
  !   Currently only planting is implemented but the addition of seeds will be added in the future.
  !=================================================================================================
  
  subroutine plant(site, patch, bc_in, pft_index, density, dbh, height) !plant_sapling()?
    ! ----------------------------------------------------------------------------------------------
    ! Plant a new cohort with the PFT, density, and size specified into an existing patch.
    ! Planting differs from recruitment in that the new plants come from somewhere else, currently
    ! not specified, rather than the seed bank.  There are currently no constrains of the size of
    ! the plantings, however the expected use will be primarily to add small plants, i.e. seedlings.
    ! ----------------------------------------------------------------------------------------------
    
    use EDCohortDynamicsMod, only : create_cohort, zero_cohort, InitPRTObject
    use EDTypesMod, only : site_massbal_type
    use PRTGenericMod, only : num_elements, num_organ_types, element_list
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: site ! Should be able to get this from patch?????
    ! I'm not sure why this must be pointer rather than target?
    type(ed_patch_type), intent(inout), pointer :: patch
    type(bc_in_type), intent(in) :: bc_in
    integer(i4), intent(in) :: pft_index ! PFT index number to plant.
    
    ! Providing the planting density, DBH, and height for the seedlings to plant is optional.
    ! Default behavior is provided but is likely that in most cases a density and dbh or height
    ! will need to be provided to get the desired behavior.
    
    ! If the density is not provided fates_recruit_initd will be used.
    ! Only DBH or height or should be provided.  If both are provided then height will be ignored.
    ! Consider adding a warning?????
    ! If neither DBH nor height is provided then fates_recruit_hgt_min will be used.
    
    ! Optional parameters:
    ! Note: It appears that the compliler on Cheyenne does not support the value keyword / attribute.
    ! If called due to a VM driver file event 'empty' values may be passed for these parameters.
    ! We check that below.
    real(r8), intent(in), optional :: density ! The planting density (plants / m^2)
    real(r8), intent(in), optional :: dbh ! Sapling diameter at breast height (cm)
    real(r8), intent(in), optional :: height ! Sapling height (m)
    
    ! Locals:
    real(r8) :: the_density ! The actual planting density (plants / m^2)
    real(r8) :: the_dbh ! The actual DBH (cm)
    real(r8) :: the_height ! The actual height (m)
    real(r8) :: area ! Patch / cohort area (m^2)
    logical :: use_dbh ! Specify cohort size by DBH, otherwise use height
    type(ed_cohort_type), pointer :: planting_cohort ! Temporary cohort for calculations
    class(prt_vartypes), pointer :: planting_parteh ! PARTEH object to hold elemental states
    type(site_massbal_type), pointer :: site_mass ! Used to access masses for flux accounting
    real(r8) :: perplant_mass ! The total elemental mass for the cohort
    integer :: el ! Element number counter
    integer :: element_id ! Element ID number
    integer :: organ_id ! Organ ID counter
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Initialize variables:
    the_density = 0
    the_dbh = 0
    the_height = 0
    area = patch%area
    
    ! If not provided use the default initial plant density:
    if (present(density) .and. (density /= vm_empty_real)) then
      the_density = density
    else
      the_density = EDPftvarcon_inst%initd(pft_index)
    end if
    
    ! Size:
    ! If DBH or height if provided (in that order), otherwise use the default height.
    if (present(dbh) .and. (dbh /= vm_empty_real)) then
      the_dbh = dbh
      use_dbh = .true.
    else
      if (present(height) .and. (height /= vm_empty_real)) then
        the_height = height
      else
        ! Get the default from the parameter file:
        the_height = EDPftvarcon_inst%hgt_min(pft_index)
      end if
      
      use_dbh = .false.
      
      ! Calculate the plant diameter from height:
      call h2d_allom(the_height, pft_index, the_dbh)
    end if
    
    if (debug) then
      write(fates_log(),*) 'FatesVegetationManagementMod: plant() running with:'
      write(fates_log(), '(A,F5.3)') 'Density =', the_density
      write(fates_log(), '(A,F5.3)') 'DBH =', the_dbh
      write(fates_log(), '(A,F5.3)') 'Height =', the_height
      write(fates_log(), '(A,F5.3)') 'Patch area =', area
      write(fates_log(), '(A,L)') 'Use DBH =', use_dbh
    end if
    
    ! Instantiate a cohort object with the desired properties:
    allocate(planting_cohort)
    
    ! Zeroing out the cohort may not be necessary and only EDPhysiologyMod: recruitment() does it
    ! but it seems like a good idea.
    call zero_cohort(planting_cohort)
    
    ! Create a PARTEH object to hold the rest of the cohorts states:
    planting_parteh => null()
    call InitPRTObject(planting_parteh)
    
    if (use_dbh) then
      call init_temporary_cohort(site, planting_cohort, planting_parteh, &
                                 pft_index, the_density, area, dbh = the_dbh)
    else
      call init_temporary_cohort(site, planting_cohort, planting_parteh, &
                                 pft_index, the_density, area, height = the_height)
    end if
    
    ! Use the cohort's mass quantities to determine the amounts to bring in from the generic fluxes:
    ! Note: Follows code in EDPhysiologyMod: recruitment().
    do el = 1,num_elements
      element_id = element_list(el)
      
      if (debug) then
        write(fates_log(), '(A,I)') 'Element number = ', el
        write(fates_log(), '(A,I)') 'Element ID = ', element_id
      end if
      
      ! For each element we get the total mass across all organs (& subpools / positions) within it:
      ! Note: In EDPhysiologyMod: recruitment() this is done manually:
      ! perplant_mass = (m_struct + m_leaf + m_fnrt + m_sapw + m_store + m_repro)
      perplant_mass = 0
      do organ_id = 1,num_organ_types
        perplant_mass = perplant_mass + planting_parteh%GetState(organ_id, element_id)
      end do
      
      if (debug) then
        write(fates_log(), '(A,F8.3)') 'Per-plant mass = ', perplant_mass
      end if
      
      site_mass => site%mass_balance(el)
      site_mass%flux_generic_in = site_mass%flux_generic_in + (planting_cohort%n * perplant_mass)
    end do
    
    ! Add the cohort to the patch:
    ! clayer:
    ! Assume we are planting on bareground and are therefore in the top layer. This is fine for now
    ! but will not be a safe assumption in all cases.  If a canopy is established a value of clayer
    ! = 2 may be better.
    ! recruitstatus:
    ! I'm not sure whether to consider these recruits.
    ! EDPhysiologyMod: recruitment(): 1.  This is straightforward.
    ! FatesInventoryInitMod: set_inventory_edcohort_type1() uses a value of 0, which makes sense for
    ! inventory (note: it has constants rstatus and recruitstatus, the later of which is unused).
    ! However, the flag is not used much so I will start by treating them as recruits.
    call create_cohort(site, patch, planting_cohort%pft, planting_cohort%n, planting_cohort%hite, &
                       planting_cohort%coage, planting_cohort%dbh, planting_parteh, &
                       planting_cohort%laimemory, planting_cohort%sapwmemory, &
                       planting_cohort%structmemory, planting_cohort%status_coh, &
                       1, planting_cohort%canopy_trim, 1, site%spread, bc_in)
                       ! Can't mix labeled and unlabeled arguments.
                       ! recruitstatus = 1, planting_cohort%canopy_trim, &
                       ! clayer = 1, site%spread, bc_in)
    
    deallocate(planting_cohort) ! Free the temporary cohort.
    
  end subroutine plant

  !=================================================================================================
  
  ! Revise name????? InitUnattached cohort?
  subroutine init_temporary_cohort(site, cohort_obj, cohort_parteh, pft, density, cohort_area, dbh, height) 
    ! ----------------------------------------------------------------------------------------------
    ! Populate a cohort based on the the cohort statistics passed in.
    !
    ! Code to initialize a cohort occurs in several places including:
    ! - EDInitMod: init_cohorts()
    ! - EDPhysiologyMod: recruitment()
    ! - FatesInventoryInitMod: set_inventory_edcohort_type1()
    ! I extracted the common code from these subroutines.  The code is modified primarily from 
    ! FatesInventoryInitMod: set_inventory_edcohort_type1().  I have mainly altered indenting,
    ! whitespace, some variable names, and added some comments.
    ! This code is called from a loop over PFTs and in some cases patches in the existing
    ! subroutines, but for planting it may be only called once.
    !
    ! If this subroutine persists it could be modified to replace the duplicated code with a few
    ! changes. To replace the current code this subroutine should be split into two halves.  I have
    ! made notes about the differences to facilitate this.
    ! ----------------------------------------------------------------------------------------------
    
    use EDTypesMod, only : leaves_on, leaves_off
    use EDTypesMod, only : phen_cstat_nevercold, phen_cstat_iscold
    use EDTypesMod, only : phen_dstat_timeoff, phen_dstat_moistoff
    use FatesAllometryMod, only : bagw_allom
    use FatesAllometryMod, only : bbgw_allom
    use FatesAllometryMod, only : bleaf
    use FatesAllometryMod, only : bfineroot
    use FatesAllometryMod, only : bsap_allom
    use FatesAllometryMod, only : bdead_allom
    use FatesAllometryMod, only : bstore_allom
    use FatesGlobals, only : endrun => fates_endrun
    use FatesInterfaceTypesMod, only : hlm_parteh_mode, nleafage
    use PRTGenericMod, only : SetState
    use PRTGenericMod, only : carbon12_element, nitrogen_element,  phosphorus_element
    use PRTGenericMod, only : struct_organ, leaf_organ, fnrt_organ, sapw_organ, store_organ
    use PRTGenericMod, only : repro_organ
    use PRTGenericMod, only : prt_carbon_allom_hyp, prt_cnp_flex_allom_hyp
    use PRTGenericMod, only : num_elements, element_list
    
    ! Arguments:
    ! The site is only used to check the leaf status tags cstatus and dstatus.  It might be better
    ! to pass these in or find some other way to avoid dependance on the site information.
    type(ed_site_type), intent(inout), target :: site
    !Name: temp_cohort, cohort, the_cohort?????
    type(ed_cohort_type), intent(inout), target :: cohort_obj ! Pointer to the cohort to initialize.
    class(prt_vartypes), intent(inout), target :: cohort_parteh ! Pointer to PARTEH object to initialize.
    integer(i4), intent(in) :: pft ! PFT index number to for cohort
    real(r8), intent(in) :: density ! The planting density in plants / m^2
    real(r8), intent(in) :: cohort_area ! The area occupied by the cohort
    real(r8), intent(in), optional :: dbh ! Plant diameter at breast height in cm
    real(r8), intent(in), optional :: height ! Plant height in meters
    
    ! Locals:
    ! integer :: cstatus ! Cohort leaf status  [Use object member]
    integer :: iage ! Leaf age counter
    integer :: el ! Element number counter
    integer :: element_id ! Element ID number
    real(r8) :: b_agw    ! Biomass above ground non-leaf [kgC]
    real(r8) :: b_bgw    ! Biomass below ground non-leaf [kgC]
    real(r8) :: b_leaf   ! Biomass in leaves [kgC]
    real(r8) :: b_fnrt   ! Biomass in fine roots [kgC]
    real(r8) :: b_sapw   ! Biomass in sapwood [kgC]
    real(r8) :: b_struct ! Structural biomass [kgC]
    real(r8) :: b_store  ! Storage pool carbon [kgC]
    real(r8) :: a_sapw   ! Area of sapwood at reference height [m2]
    real(r8) :: m_struct ! Generic (any element) mass for structure [kg]
    real(r8) :: m_leaf   ! Generic mass for leaf  [kg]
    real(r8) :: m_fnrt   ! Generic mass for fine-root  [kg]
    real(r8) :: m_sapw   ! Generic mass for sapwood [kg]
    real(r8) :: m_store  ! Generic mass for storage [kg]
    real(r8) :: m_repro  ! Generic mass for reproductive tissues [kg]
    real(r8) :: stem_drop_fraction
    
    ! ----------------------------------------------------------------------------------------------
    ! Part 1:
    ! Set the central properties.
    ! ----------------------------------------------------------------------------------------------
    
    ! Plant functional type:
    cohort_obj%pft = pft
    
    ! Age: (as cohort age bin, not years)
    ! EDInitMod: init_cohorts(): 0 (but set later)
    ! EDPhysiologyMod: recruitment(): 0
    ! FatesInventoryInitMod: set_inventory_edcohort_type1(): Not set
    ! There may be cases where a non-zero age is appropriate.  We could for example use the actual
    ! age of saplings used (~1yr).  If so it likely needs to be passed in because there is
    ! probably no great way to calculate it.
    cohort_obj%coage = 0.0_r8
    
    ! Number of trees / density:
    ! This should be straight forward but it is a bit confusing.
    ! EDTypesMod defines n as the "number of individuals in cohort per 'area' (10000m2 default)"
    ! If the area is invariant (as onw hectare), this is really a density.  However, both
    ! EDInitMod: init_cohorts() & FatesInventoryInitMod: set_inventory_edcohort_type1() initialize
    ! it using the patch area.  Based on previous work I think it is the number of trees per the
    ! the whole nominal hectare, even though the cohort may not occupy the whole hectare, which is
    ! not very intuitive.
    cohort_obj%n = density * cohort_area
    
    ! DBH & Height:
    ! The different functions start with either DBH or height and we allow either:
    ! EDInitMod: init_cohorts(): height
    ! EDPhysiologyMod: recruitment(): height
    ! FatesInventoryInitMod: set_inventory_edcohort_type1(): DBH

    ! If both are present allow height to override DBH:
    ! Note: It may be better to throw and error here?????
    if (present(height)) then
      cohort_obj%hite = height
      
      ! Calculate the plant diameter from height:
      call h2d_allom(height, pft, cohort_obj%dbh)
    else
      if (.not. present(dbh)) then
        ! One should be provided.
        write(fates_log(),*) 'DBH or height must be provided.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      end if
      
      cohort_obj%dbh = dbh
      
      ! Calculate the plant height from diameter:
      call h_allom(dbh, pft, cohort_obj%hite)
    end if
    
    ! The different subroutines make different assumptions for the canopy status.
    ! EDInitMod: init_cohorts() & FatesInventoryInitMod: set_inventory_edcohort_type1() = 1.0
    ! EDPhysiologyMod: recruitment() = 0.8.
    ! For planting an open canopy would be expected so I will use the lower value for now.
    ! Why not use 0?  I need to confirm that trim is the opposite of spread
    cohort_obj%canopy_trim = 1.0_r8 ! ?????
    
    ! ----------------------------------------------------------------------------------------------
    ! Part 2:
    ! Initialize the live pools.
    !
    ! The following code is functionally identical in:
    ! EDInitMod: init_cohorts() & FatesInventoryInitMod: set_inventory_edcohort_type1().
    ! Only the some variable names are different.
    ! EDPhysiologyMod: recruitment() also does the same thing in a different order and calls
    ! zero_cohort() first and only initializes sapwood and structural biomass if the PFT is woody.
    ! 
    ! I have only altered indenting & some comments.
    ! ----------------------------------------------------------------------------------------------

    ! Calculate total above-ground biomass from allometry
    call bagw_allom(cohort_obj%dbh, cohort_obj%pft, b_agw)

    ! Calculate coarse root biomass from allometry
    call bbgw_allom(cohort_obj%dbh, cohort_obj%pft, b_bgw)

    ! Calculate the leaf biomass (calculates a maximum first, then applies canopy trim
    ! and sla scaling factors)
    call bleaf(cohort_obj%dbh, cohort_obj%pft, cohort_obj%canopy_trim,b_leaf)

    ! Calculate fine root biomass
    call bfineroot(cohort_obj%dbh, cohort_obj%pft, cohort_obj%canopy_trim,b_fnrt)!...

    ! Calculate sapwood biomass
    call bsap_allom(cohort_obj%dbh, cohort_obj%pft, cohort_obj%canopy_trim, a_sapw, b_sapw)

    call bdead_allom( b_agw, b_bgw, b_sapw, cohort_obj%pft, b_struct )

    call bstore_allom(cohort_obj%dbh, cohort_obj%pft, cohort_obj%canopy_trim, b_store)

    ! The default assumption is that leaves are on:
    cohort_obj%laimemory = 0.0_r8
    cohort_obj%sapwmemory = 0.0_r8
    cohort_obj%structmemory = 0.0_r8
    !cstatus = leaves_on
    cohort_obj%status_coh = leaves_on

    stem_drop_fraction = EDPftvarcon_inst%phen_stem_drop_fraction(cohort_obj%pft)

    if (prt_params%season_decid(cohort_obj%pft) == itrue .and. &
        any(site%cstatus == [phen_cstat_nevercold, phen_cstat_iscold])) then
      cohort_obj%laimemory = b_leaf
      cohort_obj%sapwmemory = b_sapw * stem_drop_fraction
      cohort_obj%structmemory = b_struct * stem_drop_fraction	    
      b_leaf  = 0.0_r8
      b_sapw = (1.0_r8 - stem_drop_fraction) * b_sapw
      b_struct  = (1.0_r8 - stem_drop_fraction) * b_struct
      !cstatus = leaves_off
      cohort_obj%status_coh = leaves_off
    endif

    if (prt_params%stress_decid(cohort_obj%pft) == itrue .and. &
        any(site%dstatus == [phen_dstat_timeoff, phen_dstat_moistoff])) then
      cohort_obj%laimemory = b_leaf
      cohort_obj%sapwmemory = b_sapw * stem_drop_fraction
      cohort_obj%structmemory = b_struct * stem_drop_fraction	    
      b_leaf  = 0.0_r8
      b_sapw = (1.0_r8 - stem_drop_fraction) * b_sapw
      b_struct  = (1.0_r8 - stem_drop_fraction) * b_struct	    
      !cstatus = leaves_off
      cohort_obj%status_coh = leaves_off
    endif

    ! End Part 2.
    
    ! ----------------------------------------------------------------------------------------------
    ! Part 3:
    ! Some of the subroutines do something here:
    ! EDInitMod: init_cohorts() sets the age and a debug message here.
    ! EDPhysiologyMod: recruitment(): The limiting element is determined.
    ! FatesInventoryInitMod: set_inventory_edcohort_type1(): Nothing.
    ! ----------------------------------------------------------------------------------------------
    
    ! ----------------------------------------------------------------------------------------------
    ! Part 4:
    ! Initialize the mass of every element in every organ of the plant / cohort and store_organ
    ! in a PARTEH object.
    !
    ! The following code is functionally identical in:
    ! EDInitMod: init_cohorts(), EDPhysiologyMod: recruitment(),  & FatesInventoryInitMod:
    ! set_inventory_edcohort_type1(), except init_cohorts() & recruitment() put all the leaf mass
    ! into the first leaf age bin while set_inventory_edcohort_type1() spread it across all age
    ! bins.
    ! Since we are only really trying to do planting currently we will use former approach.
    ! ----------------------------------------------------------------------------------------------
    
    do el = 1, num_elements
      element_id = element_list(el)

      ! If this is carbon12, then the initialization is straight forward
      ! otherwise, we use stoichiometric ratios
      select case(element_id)
      case(carbon12_element)

        m_struct = b_struct
        m_leaf   = b_leaf
        m_fnrt   = b_fnrt
        m_sapw   = b_sapw
        m_store  = b_store
        m_repro  = 0.0_r8

      case(nitrogen_element)

        m_struct = b_struct * prt_params%nitr_stoich_p1(cohort_obj%pft, struct_organ)
        m_leaf   = b_leaf * prt_params%nitr_stoich_p1(cohort_obj%pft, leaf_organ)
        m_fnrt   = b_fnrt * prt_params%nitr_stoich_p1(cohort_obj%pft, fnrt_organ)
        m_sapw   = b_sapw * prt_params%nitr_stoich_p1(cohort_obj%pft, sapw_organ)
        m_store  = b_store * prt_params%nitr_stoich_p1(cohort_obj%pft, store_organ)
        m_repro  = 0.0_r8

      case(phosphorus_element)

        m_struct = b_struct * prt_params%phos_stoich_p1(cohort_obj%pft, struct_organ)
        m_leaf   = b_leaf * prt_params%phos_stoich_p1(cohort_obj%pft, leaf_organ)
        m_fnrt   = b_fnrt * prt_params%phos_stoich_p1(cohort_obj%pft, fnrt_organ)
        m_sapw   = b_sapw * prt_params%phos_stoich_p1(cohort_obj%pft, sapw_organ)
        m_store  = b_store * prt_params%phos_stoich_p1(cohort_obj%pft, store_organ)
        m_repro  = 0.0_r8
      end select

      select case(hlm_parteh_mode)
      case (prt_carbon_allom_hyp,prt_cnp_flex_allom_hyp )

        ! Put all of the leaf mass into the first bin
        call SetState(cohort_parteh, leaf_organ, element_id, m_leaf,1)
        do iage = 2, nleafage
          call SetState(cohort_parteh, leaf_organ, element_id, 0.0_r8, iage)
        end do

        call SetState(cohort_parteh, fnrt_organ, element_id, m_fnrt)
        call SetState(cohort_parteh, sapw_organ, element_id, m_sapw)
        call SetState(cohort_parteh, store_organ, element_id, m_store)
        call SetState(cohort_parteh, struct_organ, element_id, m_struct)
        call SetState(cohort_parteh, repro_organ, element_id, m_repro)

      case default
        write(fates_log(),*) 'Unspecified PARTEH module during inventory intitialization'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      end select

    end do

    call cohort_parteh%CheckInitialConditions()
    
    ! End Part 4.

  end subroutine init_temporary_cohort

  !=================================================================================================
  ! Managed Mortality Primitives:
  !   Mortality primitives are low level subroutines used to initiate human caused mortality events
  ! at the cohort and patch levels.  They provide an interface for higher level functions to build
  ! abstractions of human vegetation management activities without depending as intimately on the
  ! the underlying implementation of managed mortality representation.
  !
  !   The kill() functions provided share arguments and concepts...
  !
  ! Arguments:
  !   The kill() functions provide shared arguments that together control how much vegetation is
  ! killed and what happens to the dead and surviving trees.  These are the fraction arguments and
  ! flux profile argument.
  !
      ! Fraction Arguments:
  ! Fractions: kill_fraction vs. patch_fraction
  !   Two arguments control how much of a cohort or patch to kill and what effect this has on patch
  ! disturbance behavior.
      ! How much mortality is there?  How is it spread out?  Disturbing mortality all goes to a new
      ! patch.  Non-disturbing goes everywhere...
      ! Just need to figure out the math and make sure it works and allows clear fluxes...
  !
  ! 1. kill_fraction specifies the amount of the cohort to kill as a fraction of the number of
  ! plants in the cohort.
  !
  ! 2. patch_fraction specifies the area of the patch (as a fraction) to remove the mortality from.
    ! It does not modify the number of plants that will be killed, it just gives them context.  It
    ! can also be thought of as an evenness factor.  Are the plants being removed evenly from the
    ! whole patch or more from one side.
    !   Since patched are inherently spatially homogenous this value can only influence the patch
    ! structure through splitting.  Thus the patch_fraction is also a request to break off a new
    ! patch of this (fractional) size.
    ! 
    ! Note: This is straight forward to think about at the patch level but less so at the cohort level...
  
  ! Fraction examples:
    !   For simplicity of visualization assume all the cohorts have the same number of plants.
    !
    ! Common starting point for examples:
    ! Cohort
    !      1  TTTTTTTTTTTTTTTTTTTT
    !      2  TTTTTTTTTTTTTTTTTTTT
    !      3  TTTTTTTTTTTTTTTTTTTT
    
            ! Example 1: kill()
            !      1  TTTTTTTTTTTTTTTTTTTT
            !      2  TTTTTTTTTT_T_T_T_T_T  kill_fraction = 0.5, patch_fraction = 0.5
            !      3  TTTTTTTTTT___T___T__  kill_fraction = 0.8, patch_fraction = 0.5
            !
            ! Patch splitting:
            !      1  TTTTTTTTTT    TTTTTTTTTT
            !      2  TTTTTTTTTT    _T_T_T_T_T
            !      3  TTTTTTTTTT    ___T___T__
    
! Condensed form:
    ! Example 1:                    -Fractions-
    ! Cohort  -----Patch Area-----  kill  patch         -----Split Patches-----
    !      1  TTTTTTTTTTTTTTTTTTTT                      TTTTTTTTTT   TTTTTTTTTT
    !      2  TTTTTTTTTT_T_T_T_T_T   0.25  0.5     =>   TTTTTTTTTT   _T_T_T_T_T
    !      3  TTTTTTTTTT___T___T__   0.4   0.5          TTTTTTTTTT   ___T___T__
    !
    !
    ! Example 2:                    -Fractions-
    ! Cohort  -----Patch Area-----  kill  patch         -----Split Patches-----
    !      1  TTT_TTT_TTT_TTT_TTT_   0.25  1.0          TTT_TTT_TT   T_TTT_TTT_
    !      2  T_T_T_T_T_T_T_T_T_T_   0.5   1.0     =>   T_T_T_T_T_   T_T_T_T_T_
    !      3  TTTTTTTTTT___T___T__   0.4   0.5          TTTTTTTTTT   ___T___T__  Although both patches should be "new" here!
    !
    !
    ! Example 3:                    -Fractions-
    ! Cohort  -----Patch Area-----  kill  patch         -----Split Patches------
    !      1  TTTTTTTTTTTTTTTTTTTT                      TTTTTTTT  TTTTTT  TTTTTT
    !      2  TTTTTTTT_T_T_T_T_T_T   0.3   0.6     =>   TTTTTTTT  _T_T_T  _T_T_T
    !      3  TTTTTTTTTTTTTT______   1.0   0.3          TTTTTTTT  TTTTTT  ______
    !
    !
  
  !
  ! (Death) Flux Profiles:--------------------------------------------------------------------------
  ! [Note: consider changing the name to reflect that is about practices disturbances as much as fluxes.]
  !   Flux profiles control what happens to plant materials after death.  Flux profiles succinctly
  ! and clearly state what happens to necromass and are used to specify and track mortalities,
  ! execute resulting fluxes, and inform patch disturbance.  More that one management activity may
  ! map to the same profile.  Identifying an existing or making a new profile is also a useful part
  ! of developing a management activity or regime.
  !
  !   At this stage I'm not trying to create an exhaustive list of flux profiles or get everything
  ! right.  I just need to get the structural details figure our and draft behaviors for my test
  ! scenarios.
  !
  ! Initial profiles:
  ! - Logging module classic / legacy:
  !   Depending on the logging module settings...
  !
  ! - Harvest (new) [class]
  !     For the new vegetation harvest routines we make some changes from the logging module.
  ! Initial we won't worry about calculating collateral and infrastructure mortality.
  ! We don't assume that wood will only be harvested from trees in the top canopy layer.
  ! Disturbance calculation differences?????
  !
  !   Use for thinning and new harvest subroutines.
  !
  !   There will ultimately be multiple flavors but start with bole harvest:
  ! Bole_harvest:
  ! The bole is removed to the wood product pool.
  ! The remaining aboveground biomass = slash ->  course woody debris.
  ! Below-ground -> soil (consider stump sprouting in future).
      ! Is disturbing:
      ! Disturbance is proportional to area logged = fractional footprint of canopy trees...
      ! Unlike in the logging module trees can be harvested from any layer...
  !
  ! - In place:
  !   Kill and leave the dead.the plant where it is.
  ! Fluxes:
  ! Aboveground -> CWD & litter
  ! Below-ground -> soil
  !   Use for understory control.  Currently there is no real way to distinguish between mowing and
  ! herbicide.  Maybe in time there we will be a standing dead pool and we can flow herbicide
  ! through that pool.
  !
  ! - Burn: [ToDO]
  !  Some part of the aboveground (a fraction TBD) to atmosphere, the rest to litter and soil.
  !
  !
  !
  ! Sequential mortality and effective values:------------------------------------------------------
  !
      !Subsequent mortalities operate on a reduced (effective) number of individuals but on the same
      ! area, potentially....
  !
  ! Natural mortalities are considered non-disturbing and will have a 'virtual patch fraction' of 1?????
  
  !
  ! Sequential mortality scenarios:
    !   For simplicity we show a single cohort patch (or a patch were all cohorts experience the
    ! same mortalities).
    !
    ! Example 1:                     -Fractions-
    ! -Stage-  -----Patch Area-----  kill  patch
    ! Initial  TTTTTTTTTTTTTTTTTTTT               Prior to mortality
    ! 1st      TTT_TTT_TTT_TTT_TTT_   0.25  1.0   25% thinning
    ! 2nd      T_T_T_T_T_T_T_T_T_T_   0.333 1.0   1/3 mortality => (1 - 0.25) * 1/3 = 0.5 effective
    !
    ! This resolves to one patch so it is not disturbing...
    
    ! Example 2:                     -Fractions-
    ! -Stage-  -----Patch Area-----  kill  patch
    ! Initial  TTTTTTTTTTTTTTTTTTTT               Prior to mortality
    ! 1st      TTT_TTT_TTT_TTT_TTT_   0.25  1.0   25% thinning       [Not divisible by 2! t = 1/2T]
    ! 2nd      TTT_TTT_Tt__________   0.333 0.5   1/3 mortality => (1 - 0.25) * 1/3 = 0.5 effective
    !
    ! This resolves to two patches unambiguously because the first mortality occurred evenly so that
    ! it didn't matter 'where' in the patch subsequent mortality occurs.
    
    ! Example 3:                     -Fractions-
    ! -Stage-  -----Patch Area-----  kill  patch
    ! Initial  TTTTTTTTTTTTTTTTTTTT               Prior to mortality
    ! 1st      TTTTTTTTTT_T_T_T_T_T   0.25  0.5   50% thinning on 1/2 of the patch
    ! 2nd      ????????????????????   0.333 0.5
    !  This scenario is problematic for several reasons.  First
    
        ! -----Stage-----    -----Patch Area-----  kill  patch
        ! Initial cohort     TTTTTTTTTTTTTTTTTTTT
        ! Mortality event 1  T_T_T_T_T_T_T_T_T_T_   0.5   1.0
        ! Mortality event 2  __T___T___T___T___T_   0.5   1.0
        
				! Example:
				! TTTTTTTTTTTTTTTTTTTT : Prior to mortality
				! TTT_TTT_TTT_TTT_TTT_ : 25% thinning, non-disturbing
				! T_T_T_T_T_T_T_T_T_T_ : 1/3 mortality on top of that... = 50% total, non-disturbing
				! Or:
				! TTT_TTT_TTT_________ : 100% mortality on 40% of patch, disturbing
  
        ! Example:
        ! TTTTTTTTTTTTTTTTTTTT : Prior to mortality
        ! 111111TTTTTTTTTTTTTT : 1st mortality = 0.3 (= 30%) => 0.3
        ! 1111112222222TTTTTTT : 2nd mortality = 0.5 => (1 - 0.3) * 0.5 = 0.35 effective
  !
  !=================================================================================================

  subroutine kill(cohort, flux_profile, kill_fraction, area_fraction) ! Argument order?
    ! ----------------------------------------------------------------------------------------------
    ! 'Kill' a fraction of the specified cohort.
    ! This is accomplished by staging a mortality fraction...
    ! This mortality fraction is calculated on top of any mortality which has already been staged
    ! during this time step. (More below...)
    ! The mortality is classified in a manner that conveys both the fate of the dead plant material
    ! and its disturbing effect.
    !
    ! Mechanically, kill at the cohort level is pretty trivial!
    ! At the cohort level the PFT and size is implied.  Higher level calls should do any necessary checking.
    !
    ! kill() may be called more than once in a single timestep...
    ! Need to make aware of the effective state of the cohort, or pass in a value of n to cull
    ! rather than a fraction, or change the way cohorts are updated.
    !
    ! Collateral damage should be handled at a higer level?
    !
    !
    ! Return harvest?
    ! What about fract_ldist_not_harvested?
    !
    ! Move:
    ! Note: Calculating effective states and manipulating them could be done with matrix operations
    ! but that would require copying of all cohort states and copying them back when done.  I'm
    ! not sure that would be any more efficient and it seems error prone. A second alternative is to
    ! make set of duplicate cohorts that are modified.  That definitely seems costly.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_cohort_type), intent(inout), target :: cohort ! The cohort to apply mortality to.
    integer, intent(in) :: flux_profile ! A code specifying where the dead material goes.
    real(r8), intent(in) :: kill_fraction ! Fraction of plants to remove.
    ! Make optional and default to 1.0 / everything?????
    
    ! The fraction of the patch that the mortality is removed from / the size of the disturbance:
    real(r8), intent(in), optional :: area_fraction ! Not currently handled as optional!!!!!
    
    ! Locals:
    real(r8) the_area_fraction
    real(r8), dimension(3) :: prev_mortalities
    integer :: num_mortalities ! How many different mortality factions have already been staged?
    real(r8) :: prev_vm_mortalities
    
    ! ----------------------------------------------------------------------------------------------
    
    ! If not provided assume the mortality is spread over the entire patch:
    if (present(area_fraction)) then
      the_area_fraction = area_fraction
    else
      the_area_fraction = 1.0_r8
    endif
    
    num_mortalities = 0 ! Optional...
    prev_mortalities = [cohort%lmort_direct, cohort%vm_pfrac_in_place, &
                        cohort%vm_pfrac_bole_harvest]
    
    ! Check that the fractions are valid values:
    if (kill_fraction <= 0.0_r8 .or. kill_fraction > 1.0_r8) then
      write(fates_log(),*) 'Invalid value for kill_fraction argument.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    if (the_area_fraction <= 0.0_r8 .or. the_area_fraction > 1.0_r8) then
      write(fates_log(),*) 'Invalid value for area_fraction argument.', the_area_fraction
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    ! Since kill_fraction and area_fraction are both relative to the cohort/patch kill_fraction
    ! limits the amount that can be practically removed.
    ! Note: Currently the flux code will try to execute requests that don't make sense leading to 
    ! carbon balence errors.
    if (kill_fraction > the_area_fraction) then
      write(fates_log(),*) "The fraction of cohort vegetation being killed can't exceed the area fraction."
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
        
    ! Only allow for one mortality type that matches that being requested here:
    num_mortalities = count(prev_mortalities /= 0)
    if (num_mortalities > 1) then
      write(fates_log(),*) 'There is more that one mortality staged.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
      
    else if (num_mortalities == 1) then
      ! Make sure the flux_profile matches the existing one:
      ! In the future differing profiles may be allowed.
      if (flux_profile /= get_flux_profile(cohort)) then
        write(fates_log(),*) 'Previous flux profile does not match the requested one.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
    endif ! (num_mortalities > 1)
    
    ! If there was a previous mortality applied this will overwrite the existing data so make sure
    ! you know what you are doing.  Consider using kill_disturbed().
    
    ! Initialize the kill process by recording the event attributes in the cohort object:
    select case (flux_profile)
      case (logging_traditional) ! Legacy...
        ! Traditional logging module events use LoggingMortality_frac().
        write(fates_log(),*) 'kill() is not for use with traditional logging event. '
        call endrun(msg = errMsg(__FILE__, __LINE__))
      
      case (in_place)
        
        ! Draft multi-mortality code:
!         if (cohort%vm_mort_in_place == 0)
!           cohort%vm_mort_in_place = kill_fraction
!         else
!           cohort%vm_mort_in_place = (1 - cohort%vm_mort_in_place) * kill_fraction
!           ! if cohort%vm_mort_in_place = 1 then it will be impossible to kill more, which should
!           ! probably generate and error.
!         endif
!         
!         if (cohort%vm_pfrac_in_place == 0.0_r8 .or. cohort%vm_pfrac_in_place == 1.0_r8)
!           cohort%vm_pfrac_in_place = the_area_fraction
!         else if (the_area_fraction /= 1.0_r8)
!           ! The previous mortality and this mortality both effect an area less than one, which
!           ! presents and ambiguous result. 
!         endif
        ! if (the_area_fraction == 1.0_r8) leave cohort%vm_pfrac_in_place unchanged.
        
        cohort%vm_mort_in_place = kill_fraction
        cohort%vm_pfrac_in_place = the_area_fraction
        
      case (bole_harvest)
        ! We may be able to use cohort%lmort_direct here but that will require care.
        ! For now use a new profile.
        
        cohort%vm_mort_bole_harvest = kill_fraction
        cohort%vm_pfrac_bole_harvest = the_area_fraction
        
      ! case (burn)
        ! Placeholder.
      
      case default
        write(fates_log(),*) 'Unrecognized flux profile.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
    end select
    
  end subroutine kill

!=================================================================================================

  subroutine kill_disturbed(cohort, flux_profile, kill_fraction) ! , area_fraction)
    ! ----------------------------------------------------------------------------------------------
    ! Apply a mortality event to the disturbed subset of the current cohort.
    ! This routine treats any staged mortalities like they have already happened and kills only the
    ! surviving plants on the nascent disturbed patch.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_cohort_type), intent(inout), target :: cohort ! The cohort to apply mortality to.
    integer, intent(in) :: flux_profile ! A code specifying where the dead material goes.
    real(r8), intent(in) :: kill_fraction ! Fraction of plants to remove from disturbed sub-patch.
    ! Make optional and default to 1.0 / everything?????
    
    ! If we include the area_fraction argument it needs to be interpreted differently here.
    
    ! Locals:
    integer :: num_mortalities ! How many different mortality factions have already been staged?
    real(r8), dimension(2) :: prev_vm_mortalities ! num_profiles?????
    real(r8), dimension(2) :: prev_area_fractions
    
    real(r8) :: prev_mortalilty
    real(r8) :: prev_pfrac
    real(r8) :: total_mortality
    
    ! ----------------------------------------------------------------------------------------------
    
    num_mortalities = 0 ! Optional...
    prev_vm_mortalities = [cohort%vm_mort_in_place, cohort%vm_mort_bole_harvest]
    prev_area_fractions = [cohort%vm_pfrac_in_place, cohort%vm_pfrac_bole_harvest]
    
    ! Validity checking:
    if (kill_fraction <= 0.0_r8 .or. kill_fraction > 1.0_r8) then
      write(fates_log(),*) 'Invalid value for kill_fraction argument.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    ! Determine if any mortalities have been staged:
    
    ! We only allow vegetation management mortalities:
    if (cohort%lmort_direct > 0.0_r8) then
      write(fates_log(),*) "Can't currently mix traditional logging module events with other vegetation management."
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    ! Only allow for one mortality type that matches that being requested here:
    num_mortalities = count(prev_area_fractions /= 0)
    if (num_mortalities > 1) then
      write(fates_log(),*) 'There is more that one mortality staged.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    if (num_mortalities == 0) then
      ! If no mortality has previously occured just pass the values to kill():
      call kill(cohort, flux_profile, kill_fraction) ! , area_fraction) ! Need names?????
      ! Or make changes directly?
      
    else if (num_mortalities == 1) then
      ! Make sure the flux_profile matches the existing one:
      ! In the future differing profiles may be allowed.
      if (flux_profile /= get_flux_profile(cohort)) then
        write(fates_log(),*) 'Previous flux profile does not match the requested one.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
      
      ! Should also check the area fraction matches pfrac!
      
      
      ! Review!!!!!:
      ! Convert the kill_fraction from numbers relative to the nasent disturbed patch to the numbers
      ! relative to full cohort:
      
      ! disturbed_n * kill fraction / n
      ! rescaled_fraction = disturbed_n(cohort) * kill_fraction / cohort%n
      
      ! disturbed_n = n in area that is disturbed * mortality rate
      ! disturbed_n = n * pfrac * mort
      
      
      
      ! Convert kill_fraction to a cohort-wide mortality rate and add it the existing staged rate:
      prev_mortalilty = maxval(prev_vm_mortalities)
      prev_pfrac = maxval(prev_area_fractions)
      total_mortality = prev_mortalilty + kill_fraction * (prev_pfrac - prev_mortalilty)
      
      ! I suspect there is a problem with this calculation: [Probably no longer necessary]
      if (debug) then
        write(fates_log(), *) 'prev_mortalilty = ', prev_mortalilty
        write(fates_log(), *) 'prev_pfrac = ', prev_pfrac
        write(fates_log(), *) 'total_mortality = ', total_mortality
        
        write(fates_log(), *) 'cohort%vm_mort_in_place = ', cohort%vm_mort_in_place
        write(fates_log(), *) 'cohort%vm_mort_bole_harvest = ', cohort%vm_mort_bole_harvest
        write(fates_log(), *) 'prev_vm_mortalities = ', prev_vm_mortalities
        write(fates_log(), *) 'cohort%vm_pfrac_in_place = ', cohort%vm_pfrac_in_place
        write(fates_log(), *) 'cohort%vm_pfrac_bole_harvest = ', cohort%vm_pfrac_bole_harvest
        write(fates_log(), *) 'prev_area_fractions = ', prev_area_fractions
      end if
      
      ! Call kill() or make the changes directly?
      ! Passing it to kill() reduces the code overhead but it may require us to break the validity checks.
      !call kill(cohort = cohort, flux_profile = flux_profile, kill_fraction = total_mortality)
      call kill(cohort = cohort, flux_profile = flux_profile, kill_fraction = total_mortality, &
                area_fraction = 1.0_r8) ! Temporarily hardwire!!!!!
      
    endif ! (num_mortalities == X)
  end subroutine kill_disturbed

  !=================================================================================================
  
  ! Argument order?????
  subroutine kill_patch(patch, flux_profile, pfts, dbh_min, dbh_max, ht_min, ht_max, &
                        kill_fraction, area_fraction) ! kill_patch_kill_disturbed? cull? harvest_patch
    ! ----------------------------------------------------------------------------------------------
    ! Kill plants that match the PFT and size specifications passed across cohorts of a patch.
    !
    ! I need to add a disturbed flag or a disturbed version!
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: None
    
    ! Arguments:
    type(ed_patch_type), intent(inout), target :: patch
    integer, intent(in) :: flux_profile ! A code specifying where the dead material goes.
    integer(i4), dimension(:), intent(in) :: pfts ! Optional?
    
    ! Size specification:
    ! Size range of to plants to kill.  Defaults to everything, otherwise some range of sizes, e.g. > 10cm DBH, < 15 m in height, etc.
    real(r8), intent(in), optional :: dbh_min
    real(r8), intent(in), optional :: dbh_max
    real(r8), intent(in), optional :: ht_min
    real(r8), intent(in), optional :: ht_max
    
    real(r8), intent(in) :: kill_fraction ! Fraction of plants to remove from disturbed sub-patch.
    ! Make optional and default to 1.0 / everything?????
    
    ! The fraction of the patch that the mortality is removed from / the size of the disturbance:
    real(r8), intent(in) :: area_fraction ! Make optional?
    
    ! Locals:
    real(r8) :: the_dbh_min
    real(r8) :: the_dbh_max
    real(r8) :: the_ht_min
    real(r8) :: the_ht_max
    
    type(ed_cohort_type), pointer :: current_cohort
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Check and set the size specifications:
    call validate_size_specifications(the_dbh_min, the_dbh_max, the_ht_min, the_ht_max, &
                                      dbh_min, dbh_max, ht_min, ht_max)
    
    ! Similar to understory_control?????
    current_cohort => patch%shortest
    do while(associated(current_cohort))
    
      if (any(pfts == current_cohort%pft) .and. &
          current_cohort%dbh >= dbh_min .and. current_cohort%dbh <= dbh_max .and. &
          current_cohort%hite >= ht_min .and. current_cohort%hite <= ht_max) then
        
        call kill(cohort = current_cohort, flux_profile = flux_profile, &
                  kill_fraction = kill_fraction, & area_fraction = area_fraction)
        endif
      current_cohort => current_cohort%taller
    end do ! Cohort loop.
    
  end subroutine kill_patch

  !=================================================================================================
  ! Management Operations:
  !   Operations are specific (atomci) real world management interventions that include activities
  ! such as planting, fertilization, and harvest
  !   These routines build on primitives to implement operations ...
  !
  ! Establishment:
  ! - Seeding and planting
  !
  ! Intermediate Operations:
  ! - Competition control
  ! - Fertilization
  ! - Thinning
  !
  ! Harvest
  !
  ! Site Level Routines:
  !   These routines perform management activities at the level of the site / grid cell.
  ! They may target only part of the site but they do so without any a priori knowledge of the patch
  ! structure.
  !=================================================================================================

  !=================================================================================================
  ! Establishment Subroutines:
  !=================================================================================================

  subroutine plant_site(site, bc_in, pfts, density, dbh, height) ! where = "everywhere"
    ! ----------------------------------------------------------------------------------------------
    ! Plant the specified pft across all patches of a site.
    !
    ! In the future we may allow targeting of thinning behavior to specific patches via the 'where'
    ! parameter.
    !
    ! VM Event Interface plant(pfts, [density], [dbh], [height])
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: site
    type(bc_in_type), intent(in) :: bc_in
    ! The PFT index number to plant.  Allows an array but only on will be used (see pft_index):
    integer(i4), intent(in), target :: pfts(:)
    
    ! Optional parameters:
    ! If called due to a VM driver file event 'empty' values may be passed for these parameters.
    ! The question is whether to catch that in managed_fecundity(), here, or in plant(). Handling
    ! it in plant() seemed simplest.
    real(r8), intent(in), optional :: density ! The planting density (plants / m^2)
    real(r8), intent(in), optional :: dbh ! Sapling diameter at breast height (cm)
    real(r8), intent(in), optional :: height ! Sapling height (m)
    
    ! Add 'where' parameter.
    
    ! Locals:
    type(ed_patch_type), pointer :: current_patch
    type(ed_cohort_type), pointer :: current_cohort
    integer(i4) :: pft_index
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'plant_site() beginning.'
    
    ! Many other driver events take more than one PFTs.  We could loop over multiple PFTs if that
    ! functionality were desired.  Now we just pass on one:
    if (size(pfts) > 1) then
      write(fates_log(),*) 'plant_site() currently only accepts a single PFT.  Only the first will be used.'
    endif
    
    pft_index = pfts(1)
    
    current_patch => site%oldest_patch
    do while (associated(current_patch))
      
      call plant(site = site, patch = current_patch, bc_in = bc_in, pft_index = pft_index, &
                 density = density, dbh = dbh, height = height)
      
      current_patch => current_patch%younger
    end do ! Patch loop.
    
  end subroutine plant_site

  !=================================================================================================
  ! Competition Control Subroutines:
  !=================================================================================================

  subroutine understory_control(patch, method) ! REVIEW! Add area_fraction!!!!!
    ! ----------------------------------------------------------------------------------------------
    ! Remove all understory plants, grasses and shrubs, from the specified patch using the specified
    ! method.
    ! We treat understory control as a composition altering but not disturbance generating.
    ! We would need to split a patch if we wanted to apply understory_control to part of a patch.
    !
    ! Could add option to remove small trees as well, possibly by PFT.
    !
    ! Note: This was written before kill_patch() and could be rewritten to use it.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use FatesInterfaceTypesMod, only : numpft
    
    ! Arguments:
    type(ed_patch_type), intent(inout), target :: patch
    integer, intent(in) :: method ! A code value that controls the resulting fluxes.
    
    ! Locals:
    type(ed_cohort_type), pointer :: current_cohort
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Check that the patch is secondary forest!!!!!
    
    ! Kill the understory component:
    ! Consider adding and efficiency factor => fraction argument in kill().
    select case(method)
      case(method_mow, method_herbicide)
        ! For now these should be the same but in the future they could differ if there were a
        ! standing dead stock.
        
        current_cohort => patch%shortest
        do while(associated(current_cohort))
          if (any(current_cohort%pft == understory_pfts)) then
            call kill(cohort = current_cohort, flux_profile = in_place, kill_fraction = 1.0_r8, &
                      area_fraction = 1.0_r8) ! Leaving out the area_fraction right now won't work.  Fix that.
          endif
          current_cohort => current_cohort%taller
        end do ! Cohort loop.
        
      case(method_burn)
        write(fates_log(),*) 'understory_control(): method_burn not yet implemented.'
        
      case default
        write(fates_log(),*) 'understory_control(): Invalid method requested.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
    end select
    
    ! Log the event:  Master process only?  Debug only??????
    write(fates_log(),*) 'Vegetation management: Understory control performed...'
    
  end subroutine understory_control

  !=================================================================================================

  ! Are shorthand routines understory_mow(), understory_treat/herbicide(), and understory_burn() needed?
!   subroutine understory_mow(patch)
!     ! ----------------------------------------------------------------------------------------------
!     ! 
!     ! ----------------------------------------------------------------------------------------------
!     
!     ! Uses:
!     
!     ! Arguments:
!     type(ed_patch_type), intent(inout), target :: patch
!     
!     ! Locals:
!     
!     ! ----------------------------------------------------------------------------------------------
!     
!     call understory_control(patch, method_mow)
!     
!   end subroutine understory_mow

  !=================================================================================================
  ! Thinning Subroutines:
  !   Thinning is an intermediate operation in forest management designed to reduce uncontrolled
  ! mortality from natural self-thinning and allow for removal of less vigorous and poorly formed
  ! individuals.
  !
  !  In reality some trees removed in a thinning may too small or of insufficient quality to yeild
  ! lumber.  Therefore some or all harvest from a thinning be used for pulp, or other uses (biomass
  ! is an emerging possible use).  It is also possible that some small trees may be left on site.
  !
  ! There is more than one way to perform a thinning.  The simplest example would be to remove a
  ! fixed percentage.  In practice the methods for thinning are a function of management goals, the
  ! the stand structure and available machinery.
  ! I'm starting with a thinning methodology that is common for southern pines both because I need
  ! it and because it is complex enough that it will force me to figure out most of the components
  ! I will need to build other thinning types.
  !  
  !=================================================================================================

  subroutine thin_row_low(patch, pfts, row_fraction, patch_fraction, &
                          final_basal_area, final_stem_density, harvest_estimate) ! Return the harvest amount!
    ! ----------------------------------------------------------------------------------------------
    ! Perform a row thinning followed by a low thinning to achieve a desired final basal area or
    ! stem density.
    !
    ! Define row and low thinning!!!!!
    ! Wood is harvested, no matter the size.
    !
    ! We currently add no infrastructure mortality (the row thinning is in part to allow equipment
    ! access) or collateral damage from thinning.
    !
    ! This is a fairly specific thinning routine.  It was chosen to directly mimic an actual
    ! management practice and is parameterized to be compatible with an operation perception.  The
    ! two phase thinning process yields a more complex demography and the final result will be
    ! somewhat variable, not always hitting the goal value perfectly, which is also realistic.  It
    ! alos requires an iterative algorithm, which was useful to help develop a robust set of low
    ! level utilities that will allow other practices to be developed.
    !
    ! Note: In general a low thinning can not be performed perfectly from the "bottom up" due to
    ! mechanical limitations, clustering of small and large trees, and the desire to remove ill
    ! formed or damaged trees, regardless of size.  This imperfection could be approximated using 
    ! an inverse proportional probably of taking larger trees, or something like that.
    !
    ! Ideas:
    ! This could be performed at a site level as well.  Would that just be a loop?
    ! Alternatively, supply as a row frequency?
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use FatesInterfaceTypesMod, only : numpft
    !use EDPftvarcon, only : EDPftvarcon_inst ! Will change to prt_params soon when merged to master.
    
    ! Arguments:
    type(ed_patch_type), intent(inout), target :: patch ! Or pointer?????
    ! Which tree PFTs are harvestable trees?  If omitted all woody trees will be used.  Otherwise
    ! any omitted will be ignored for the purpose of calculating BAI and stem density.
    ! Be careful when excluding PFTs that may compete in the canopy.
    integer(i4), intent(in), optional, target :: pfts(:) ! An array of PFT IDs to thin.
    real(r8), intent(in), optional :: row_fraction       ! The fraction of rows to initially thin,
                                                         ! e.g. every 5th row = 0.2, often every 4th or 5th row...
                                                         ! Alternatively, supply as a row frequency?
    real(r8), intent(in), optional :: patch_fraction     ! The fraction of the patch to thin.
                                                         ! Defaults to 1.0 / whole patch.
                                                         ! Values > 1 result in the patch being 
                                                         ! split into thinned and unthinned patches.
    real(r8), intent(in), optional :: final_basal_area   ! Goal final basal area index (m^2 / ha ?????)
    real(r8), intent(in), optional :: final_stem_density ! Goal final stem density (trees / ha)
    real(r8), intent(out), optional :: harvest_estimate  ! The wood harvested by this operation.
    
    ! Locals:
    integer(i4), pointer, dimension(:) :: thin_pfts ! Holds the PFTs to thin, computed from arguments.
    real(r8) :: the_row_fraction ! Or change row_fraction -> row_fraction_in or row_fraction_opt
    real(r8) :: the_patch_fraction
    
    logical :: use_bai ! If true use basal area index as the criteria, otherwise use stem density.
    real(r8) :: patch_bai ! Basal area of the patch (including only PFTs in pfts)
    real(r8) :: patch_sd ! Stem density of the patch (including only PFTs in pfts)
    real(r8) :: cohort_ba ! Basal area of a cohort
    real(r8) :: cohort_stems
    real(r8) :: thin_ba_remaining
    real(r8) :: thin_sd_remaining
    real(r8) :: cohort_fraction
    real(r8) :: harvest ! Accumulator
    
    type(ed_cohort_type), pointer :: current_cohort
    
    integer :: i ! Iterator
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'thin_row_low() entering.'
    
    ! Initialize variables:
    harvest = 0.0_r8
    thin_ba_remaining = 0.0_r8
    
    ! Determine which arguments were passed in, validity check them, and supply default values:
    
    ! If present check that the PFTs are valid:
    if (present(pfts)) then
      ! Check if PFTs are valid:
!       if (.not. any(tree_pfts == pfts)) then
!         write(fates_log(),*) 'thin_row_low(): Cannot thin non-tree PFTs.'
!         write(fates_log(),*) 'PFTs =', pfts
!         write(fates_log(),*) 'tree_pfts =', tree_pfts
!         write(fates_log(),*) '(tree_pfts == pfts) = ', (tree_pfts == pfts)
!         write(fates_log(),*) '(pfts == tree_pfts) = ', (pfts == tree_pfts)
!         write(fates_log(),*) 'any(tree_pfts == pfts) = ', any(tree_pfts == pfts)
!         write(fates_log(),*) 'any(pfts == tree_pfts) = ', any(pfts == tree_pfts)
!         
!         !write(fates_log(),*) 'any(woody_pfts == tree_pfts) = ', any(woody_pfts == tree_pfts)
!         write(fates_log(),*) '(2 == tree_pfts) = ', (2 == tree_pfts)
!         write(fates_log(),*) 'any(2 == tree_pfts) = ', any(2 == tree_pfts)
!         write(fates_log(),*) '([2,2,2,2,2,2] == tree_pfts) = ', ([2,2,2,2,2,2] == tree_pfts)
!         
!         !write(fates_log(),*) '(woody_pfts == tree_pfts) = ', (woody_pfts == tree_pfts)
!         
!         call endrun(msg = errMsg(__FILE__, __LINE__))
!       end if
      
      ! Revised and corrected:
      ! Check if PFTs to thin are valid:
      do i = 1, size(pfts)
        if (.not. any(pfts(i) == tree_pfts)) then
          write(fates_log(),*) 'thin_row_low(): Cannot thin non-tree PFTs.'
          write(fates_log(),*) 'Tree PFTs =    ', tree_pfts
          write(fates_log(),*) 'Selected PFTs =', pfts
          call endrun(msg = errMsg(__FILE__, __LINE__))
        endif
      end do
      
      thin_pfts => pfts
    else ! Otherwise thin all tree PFTs:
      thin_pfts => tree_pfts
    endif
    
    ! Determine the thinning amount:
    if (present(row_fraction)) then
      if (row_fraction > 1.0_r8 .or. row_fraction <= 0.0_r8) then ! Better high and low values?
        write(fates_log(),*) 'thin_row_low(): Invalid row fraction provided.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      end if
    
      the_row_fraction = row_fraction
    else
      the_row_fraction = 0.2_r8 ! Default value.
    end if
    
    ! Note: Since all the calculations below are per area the patch fraction does not effect them.
    ! It is only passed on to kill().
    if (present(patch_fraction)) then
      if (patch_fraction > 1.0_r8 .or. patch_fraction <= 0.0_r8) then
        write(fates_log(),*) 'thin_row_low(): Invalid patch fraction provided.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      end if
    
      the_patch_fraction = patch_fraction
    else
      the_patch_fraction = 1.0_r8
    end if
    
    ! Need only BAI or stem density:
    if (present(final_basal_area) .and. present(final_stem_density)) then
      write(fates_log(),*) 'thin_row_low(): Provide basal area index or stem density, not both.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    else if (present(final_basal_area)) then
      use_bai = .true.
    else if (present(final_stem_density)) then
      use_bai = .false.
    else
      write(fates_log(),*) 'thin_row_low(): Must provide basal area index or stem density.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    ! Determine the basal area and stem densities for the relavant PFTs:
    !patch_bai = effective_basal_area(patch, thin_pfts)
    !patch_sd = effective_stem_density(patch, thin_pfts)
    patch_bai = disturbed_basal_area(patch, thin_pfts)
    patch_sd = patch_disturbed_n(patch, thin_pfts) !disturbed_stem_density(patch, thin_pfts)
    
    ! If the stand is above the goal value (BAI or stem density) start thinning:
    if ((use_bai .and. (patch_bai > final_basal_area)) .or. &
        ((.not. use_bai) .and. (patch_sd > final_stem_density))) then
      
      ! Thin every X rows = remove 1/X of each cohort:----------------------------------------------
      if (debug) write(fates_log(), *) 'thin_row_low() starting row thinning.'
      
      current_cohort => patch%shortest
      do while(associated(current_cohort))
        
        if (any(pfts == current_cohort%pft)) then
          ! Call kill() here to set the area fraction.  Use kill_disturbed() for subsequent calls:
          call kill(cohort = current_cohort, flux_profile = bole_harvest, &
                    kill_fraction = the_row_fraction, area_fraction = the_patch_fraction)
          
          ! Accumulate harvest estimate:
          harvest = harvest + cohort_harvestable_biomass(current_cohort) ! staged = true!!!!
        endif
        
        current_cohort => current_cohort%taller
      end do ! Cohort loop.
      
      ! If still over the goal value thin further:--------------------------------------------------
      ! The kill() call above does not result in a change to the cohort numbers yet, this will
      ! happen later.  Therefore we need to keep track of the changes to BAI, density, and plant
      ! number's manually from here on.
      patch_bai = disturbed_basal_area(patch, thin_pfts)
      patch_sd = patch_disturbed_n(patch, thin_pfts) !disturbed_stem_density(patch, thin_pfts)
      
      ! If still above our goal BAI / density thin out the smallest trees (low thinning) recursively
      ! until the goal has been reached:
      
      if (use_bai) then ! Thin to a goal basal area:
        if (debug) write(fates_log(), *) 'thin_row_low() starting low thinning by BAI.'
        
        ! Given the fixed allometry this will also give us cohorts from the lowest DBH:
        current_cohort => patch%shortest
        do while(associated(current_cohort) .and. patch_bai > final_basal_area)
          
          if (any(pfts == current_cohort%pft)) then
            ! Get the effective (after row thinning) basal area of the cohort and determine if it
            ! can be removed in part or in whole:
            cohort_ba = disturbed_basal_area(current_cohort)
            
            ! Remaining basal area:
            thin_ba_remaining = patch_bai - final_basal_area
            
            ! If the cohort basal area is less that what still needs to be removed kill all of it:
            if (cohort_ba <= thin_ba_remaining) then
              if (debug) write(fates_log(), *) 'thin_row_low() cut whole cohort.'
              
              call kill_disturbed(cohort = current_cohort, flux_profile = bole_harvest, &
                                  kill_fraction = 1.0_r8)
              
            else ! Otherwise only take part of the cohort:
              if (debug) write(fates_log(), *) 'thin_row_low() cut part of cohort.'
              
              cohort_fraction = thin_ba_remaining / cohort_ba
              
              if (debug) then
                write(fates_log(), *) 'thin_ba_remaining = ', thin_ba_remaining
                write(fates_log(), *) 'cohort_ba = ', cohort_ba
                write(fates_log(), *) 'cohort_fraction = ', cohort_fraction
              end if
              
              call kill_disturbed(cohort = current_cohort, flux_profile = bole_harvest, &
                                  kill_fraction = cohort_fraction)
              
            end if ! (cohort_ba <= thin_ba_remaining)
            
            ! The harvest estimate and BAI only need to updated if we harvested something:
            ! Accumulate harvest estimate:
            harvest = harvest + cohort_harvestable_biomass(current_cohort) ! staged = true!!!!
            patch_bai = disturbed_basal_area(patch, thin_pfts) ! Update the BAI calculation.
          end if ! (any(pfts == current_cohort%pft))
          
          current_cohort => current_cohort%taller
        end do ! Cohort loop.
        
      !else if ((.not. use_bai) .and. (patch_sd > final_stem_density)) then
      else ! Thin to a goal stem density:
        if (debug) write(fates_log(), *) 'thin_row_low() starting row by stem density.'
        
        ! This loop is identically structured to the above so the two could be combined by
        ! generalizing the comparator variables.
        
        current_cohort => patch%shortest
        do while(patch_sd > final_stem_density)
          
          if (any(pfts == current_cohort%pft)) then
            ! Get the effective number of stems in cohort and determine if they can be removed in
            ! part or in whole:
            
            ! Because n is per nominal hectare n = stem density (n/ha)
            cohort_stems = cohort_disturbed_n(current_cohort)
            
            ! Remaining stems to remove:
            thin_sd_remaining = patch_sd - final_stem_density
            
            ! If the cohort stem count is less that what still needs to be removed kill all of it:
            if (cohort_stems <= thin_sd_remaining) then
              if (debug) write(fates_log(), *) 'thin_row_low() cut whole cohort.'
              
              call kill_disturbed(cohort = current_cohort, flux_profile = bole_harvest, &
                                kill_fraction = 1.0_r8)
              
            else ! Otherwise only take part of the cohort:
              if (debug) write(fates_log(), *) 'thin_row_low() cut part of cohort.'
              
              cohort_fraction = thin_sd_remaining / cohort_stems
              call kill_disturbed(cohort = current_cohort, flux_profile = bole_harvest, &
                        kill_fraction = cohort_fraction)
              
            end if
            
            ! The harvest estimate and stem density only need to updated if we harvested something:
            ! Accumulate harvest estimate:
            harvest = harvest + cohort_harvestable_biomass(current_cohort) ! staged = true!!!!
            patch_sd = patch_disturbed_n(patch, thin_pfts) ! disturbed_stem_density(patch, thin_pfts) ! Update the stem density.
          end if ! (any(pfts == current_cohort%pft))
          
          current_cohort => current_cohort%taller
        end do ! Cohort loop.
      endif ! (use_bai)
    endif ! Stand is above goal loop.
    ! Consider adding warning here if patch is below goal?
    
    if (present(harvest_estimate)) then
      harvest_estimate = harvest
    endif
    ! Should anything be done or reported if no thinning was needed?
    
    if (debug) write(fates_log(), *) 'thin_row_low() exiting.'
  end subroutine thin_row_low

  !=================================================================================================

  function thinnable_patch(patch, pfts, goal_basal_area) ! result(thinnable_patch) ! REVIEW!
    ! ----------------------------------------------------------------------------------------------
    ! Rough Draft!!!!!
    ! Determine if a patch is is ready to be thinned based on some set of criteria.
    !
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_patch_type), intent(in), target :: patch
    ! An array of PFT IDs to include in the basal area calculation:
    integer(i4), dimension(:), intent(in) :: pfts
    real(r8), intent(in) :: goal_basal_area ! The goal basal area that we should thin to (m^2/ha).
    
    ! Add criteria arguments.  Some reasonable criteria would be:
    !   A basal area.  There is no point in thinning if the patch is at the target already.
    ! In fact it probably needs to exceed the goal by at least 25% to be commercial.
    !   A minimum age.
    ! The PFTs to consider in the basal area.  If omitted consider all trees.
    
    ! Locals:
    logical :: thinnable_patch ! Return value.
    real(r8) :: patch_basal_area ! The basal area of the patch.
    
    ! ----------------------------------------------------------------------------------------------
    
    thinnable_patch = .false.
    
    patch_basal_area = effective_basal_area(patch, pfts) ! Add PFTs
    
    if (patch_basal_area >= goal_basal_area * 1.25_r8) then
      thinnable_patch = .true.
    endif
    
  end function thinnable_patch

  !=================================================================================================

  subroutine thin_proportional(site, pfts, thin_fraction) ! Return the harvest amount in harvest_estimate! Rename?
    ! ----------------------------------------------------------------------------------------------
    ! Thin all size classes of the PFT(s) passed proportionally across a stie, i.e by the same
    ! fraction for all sizes.  This approximates a random thinning.
    !
    ! Unlike thin_row_low() this thinning doesn't use specific thinning goals, the goal is relative.
    ! In the future we may allow targeting of thinning behavior to specific patches via the 'where'
    ! parameter.
    !
    ! VM Event Interface thin_proportional([pfts], thin_fraction)
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_site_type), intent(in), target :: site ! The current site object.
    integer(i4), dimension(:), intent(in) :: pfts ! Make optional?
    real(r8), intent(in) :: thin_fraction
    ! Add where
    
    ! Locals:
    type(ed_patch_type), pointer :: current_patch
    type(ed_cohort_type), pointer :: current_cohort
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'thin_proportional() beginning.'
    
    current_patch => site%oldest_patch
    do while (associated(current_patch))
      current_cohort => current_patch%shortest
      do while(associated(current_cohort))
        
        if (any(pfts == current_cohort%pft)) then
          ! Harvest a fraction of all trees:
          call kill(cohort = current_cohort, flux_profile = bole_harvest, &
                    kill_fraction = thin_fraction, area_fraction = 1.0_r8)
          
          ! Accumulate harvest estimate:
          !harvest = harvest + cohort_harvestable_biomass(current_cohort) ! staged = true!!!!
        endif
        
        current_cohort => current_cohort%taller
      end do ! Cohort loop.
      current_patch => current_patch%younger
    end do ! Patch loop.
    
  end subroutine thin_proportional

  !=================================================================================================

  subroutine thin_patch_low_perfect(patch, pfts, patch_fraction, thin_fraction, final_basal_area, &
                                    final_stem_density, harvest_estimate)
    ! ----------------------------------------------------------------------------------------------
    ! Thin a patch harvesting the trees from smallest to largest until the thinning goal is met.
    ! For maximum flexibility the amount to thin may be specified as a fraction of stems (relative)
    ! or a goal (stem density or basal area).
    ! The low thinning is 'perfect' because only exact number of the smallest trees are removed,
    ! which is not terrible realistic in most cases.
    !
    ! This routine duplicates and expands much of the code from thin_row_low(), except the row
    ! thinning.  thin_row_low() should probably be rewritten to use this routine.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use FatesInterfaceTypesMod, only : numpft
    
    ! Arguments:
    type(ed_patch_type), intent(inout), target :: patch
    ! Which tree PFTs are harvestable trees.  If omitted all woody trees will be used.  Otherwise
    ! any omitted will be ignored for the purpose of calculating BAI and stem density.
    ! Be careful when excluding PFTs that may compete in the canopy.
    integer(i4), intent(in), optional, target :: pfts(:) ! An array of PFT IDs to thin.
    real(r8), intent(in), optional :: patch_fraction     ! The fraction of the patch to thin.
                                                         ! Defaults to 1.0 / whole patch.
                                                         ! Values > 1 result in the patch being 
                                                         ! split into thinned and unthinned patches.
    ! Three ways to specify how much to thin:
    real(r8), intent(in), optional :: thin_fraction      ! Fraction of trees (in pfts) to thin.
    real(r8), intent(in), optional :: final_basal_area   ! Goal final basal area index (m^2 / ha ?????)
    real(r8), intent(in), optional :: final_stem_density ! Goal final stem density (trees / ha)
    
    real(r8), intent(out), optional :: harvest_estimate  ! The wood harvested by this operation.
    
    ! Locals:
    integer(i4), pointer, dimension(:) :: thin_pfts ! Holds the PFTs to thin, computed from arguments.
    real(r8) :: the_patch_fraction ! Patch fraction after application of defaults.
    real(r8) :: the_final_stem_density ! Stem density goal possibly converted from thin_fraction.
    
    logical :: use_bai ! If true use basal area index as the criteria, otherwise use stem density.
    real(r8) :: patch_bai ! Basal area of the patch (including only PFTs in pfts)
    real(r8) :: patch_sd ! Stem density of the patch (including only PFTs in pfts)
    real(r8) :: cohort_ba ! Basal area of a cohort
    real(r8) :: cohort_stems
    real(r8) :: thin_ba_remaining
    real(r8) :: thin_sd_remaining
    real(r8) :: cohort_fraction
    real(r8) :: harvest ! Accumulator
    
    type(ed_cohort_type), pointer :: current_cohort
    
    integer :: i ! Iterator
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'thin_patch_low_perfect() entering.'
    
    ! Initialize variables:
    harvest = 0.0_r8
    thin_ba_remaining = 0.0_r8
    
    ! Determine which arguments were passed in, validity check them, and supply default values:
    
    ! If present check that the PFTs are valid:
    if (present(pfts)) then
      ! Check if PFTs to thin are valid:
      do i = 1, size(pfts)
        if (.not. any(pfts(i) == tree_pfts)) then
          write(fates_log(),*) 'thin_patch_low_perfect(): Cannot thin non-tree PFTs.'
          write(fates_log(),*) 'Tree PFTs =    ', tree_pfts
          write(fates_log(),*) 'Selected PFTs =', pfts
          call endrun(msg = errMsg(__FILE__, __LINE__))
        endif
      end do
      
      thin_pfts => pfts
    else ! Otherwise thin all tree PFTs:
      thin_pfts => tree_pfts
    endif
    
    ! Note: Since all the calculations below are per area the patch fraction does not effect them.
    ! It is only passed on to kill().
    if (present(patch_fraction)) then
      if (patch_fraction > 1.0_r8 .or. patch_fraction <= 0.0_r8) then
        write(fates_log(),*) 'thin_patch_low_perfect(): Invalid patch fraction provided.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      end if
    
      the_patch_fraction = patch_fraction
    else
      the_patch_fraction = 1.0_r8
    end if
    
    ! Determine the basal area and stem densities for the relavant PFTs:
    patch_bai = disturbed_basal_area(patch, thin_pfts)
    patch_sd = patch_disturbed_n(patch, thin_pfts) !disturbed_stem_density(patch, thin_pfts)
    
    ! Only thin_fraction, final_basal_area, or final_stem_density should be provided:
    if ((present(thin_fraction) .and. present(final_basal_area)) .or. &
        (present(thin_fraction) .and. present(final_stem_density)) .or. &
        (present(final_basal_area) .and. present(final_stem_density))) then
      write(fates_log(),*) 'thin_patch_low_perfect(): Provide thinning fraction, basal area index, or stem density, not more than one.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    else if (present(final_basal_area)) then
      use_bai = .true.
    else if (present(final_stem_density)) then
      use_bai = .false.
      the_final_stem_density = final_stem_density
    else if (present(thin_fraction)) then
      ! Validity checking:
      if (thin_fraction <= 0.0_r8 .or. thin_fraction > 1.0_r8) then
        write(fates_log(),*) 'Invalid value for thin_fraction argument.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
      
      ! Convert the thinning fraction to a stem density goal:
      use_bai = .false.
      the_final_stem_density = patch_sd * (1.0_r8 - thin_fraction)
    else
      write(fates_log(),*) 'thin_patch_low_perfect(): Must provide thinning fraction, basal area index, or stem density.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    ! If the stand is above the goal value (BAI or stem density) start thinning:
    if ((use_bai .and. (patch_bai > final_basal_area)) .or. &
        ((.not. use_bai) .and. (patch_sd > the_final_stem_density))) then
      
      ! There may be previously staged mortality. Therefore we need to keep track of the changes to
      ! BAI, density, and plant number's:
      patch_bai = disturbed_basal_area(patch, thin_pfts)
      patch_sd = patch_disturbed_n(patch, thin_pfts) !disturbed_stem_density(patch, thin_pfts)
      
      ! Cull the smallest trees (low thinning) recursively until the goal has been reached:
      
      if (use_bai) then ! Thin to a goal basal area:
        if (debug) write(fates_log(), *) 'thin_patch_low_perfect() starting low thinning by BAI.'
        
        ! Given the fixed allometry this will also give us cohorts from the lowest DBH:
        current_cohort => patch%shortest
        do while(associated(current_cohort) .and. patch_bai > final_basal_area)
          
          if (any(pfts == current_cohort%pft)) then
            ! Get the effective (after row thinning) basal area of the cohort and determine if it
            ! can be removed in part or in whole:
            cohort_ba = disturbed_basal_area(current_cohort)
            
            ! Remaining basal area:
            thin_ba_remaining = patch_bai - final_basal_area
            
            ! If the cohort basal area is less that what still needs to be removed kill all of it:
            if (cohort_ba <= thin_ba_remaining) then
              if (debug) write(fates_log(), *) 'thin_patch_low_perfect() cut whole cohort.'
              
              call kill_disturbed(cohort = current_cohort, flux_profile = bole_harvest, &
                                  kill_fraction = 1.0_r8)
              
            else ! Otherwise only take part of the cohort:
              if (debug) write(fates_log(), *) 'thin_patch_low_perfect() cut part of cohort.'
              
              cohort_fraction = thin_ba_remaining / cohort_ba
              
              if (debug) then
                write(fates_log(), *) 'thin_ba_remaining = ', thin_ba_remaining
                write(fates_log(), *) 'cohort_ba = ', cohort_ba
                write(fates_log(), *) 'cohort_fraction = ', cohort_fraction
              end if
              
              call kill_disturbed(cohort = current_cohort, flux_profile = bole_harvest, &
                                  kill_fraction = cohort_fraction)
              
            end if ! (cohort_ba <= thin_ba_remaining)
            
            ! The harvest estimate and BAI only needs to be updated if we harvested something:
            ! Accumulate harvest estimate:
            harvest = harvest + cohort_harvestable_biomass(current_cohort) ! staged = true!!!!
            patch_bai = disturbed_basal_area(patch, thin_pfts) ! Update the BAI calculation.
          end if ! (any(pfts == current_cohort%pft))
          
          current_cohort => current_cohort%taller
        end do ! Cohort loop.
        
      else ! Thin to a goal stem density:
        if (debug) write(fates_log(), *) 'thin_patch_low_perfect() starting row by stem density.'
        
        ! This loop is identically structured to the above so the two could be combined by
        ! generalizing the comparator variables.
        
        current_cohort => patch%shortest
        do while(patch_sd > the_final_stem_density)
          
          if (any(pfts == current_cohort%pft)) then
            ! Get the effective number of stems in cohort and determine if they can be removed in
            ! part or in whole.  Because n is per nominal hectare n = stem density (n/ha):
            cohort_stems = cohort_disturbed_n(current_cohort)
            
            ! Remaining stems to remove:
            thin_sd_remaining = patch_sd - the_final_stem_density
            
            ! If the cohort stem count is less that what still needs to be removed kill all of it:
            if (cohort_stems <= thin_sd_remaining) then
              if (debug) write(fates_log(), *) 'thin_patch_low_perfect() cut whole cohort.'
              
              call kill_disturbed(cohort = current_cohort, flux_profile = bole_harvest, &
                                  kill_fraction = 1.0_r8)
              
            else ! Otherwise only take part of the cohort:
              if (debug) write(fates_log(), *) 'thin_patch_low_perfect() cut part of cohort.'
              
              cohort_fraction = thin_sd_remaining / cohort_stems
              call kill_disturbed(cohort = current_cohort, flux_profile = bole_harvest, &
                                  kill_fraction = cohort_fraction)
              
            end if
            
            ! The harvest estimate and stem density only need to updated if we harvested something:
            ! Accumulate harvest estimate:
            harvest = harvest + cohort_harvestable_biomass(current_cohort) ! staged = true!!!!
            patch_sd = patch_disturbed_n(patch, thin_pfts) ! disturbed_stem_density(patch, thin_pfts) ! Update the stem density.
          end if ! (any(pfts == current_cohort%pft))
          
          current_cohort => current_cohort%taller
        end do ! Cohort loop.
      endif ! (use_bai)
    endif ! Stand is above goal loop.
    ! Consider adding warning here if patch is below goal?
    
    if (present(harvest_estimate)) then
      harvest_estimate = harvest
    endif
    ! Should anything be done or reported if no thinning was needed?
    
    if (debug) write(fates_log(), *) 'thin_patch_low_perfect() exiting.'
  end subroutine thin_patch_low_perfect

  !=================================================================================================

  subroutine thin_low_perfect(site, pfts, thin_fraction) ! where = everywhere
    ! ----------------------------------------------------------------------------------------------
    ! Perform a 'perfect' low thinning for a site.
    !
    ! This routine applies thin_patch_low_perfect() to all the patches in the site using a limited
    ! set of its options.
    ! In the future we will (may) allow targeting of thinning behavior to specific patches via the
    ! 'where' parameter.
    !
    ! This routine does not yet return a harvest estimate.
    !
    ! VM Event Interface thin_low_perfect([pfts], thin_fraction)
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_site_type), intent(in), target :: site ! The current site object.
    integer(i4), dimension(:), intent(in) :: pfts ! An array of PFT IDs to thin.
    real(r8), intent(in) :: thin_fraction ! Fraction of trees to thin.
    
    ! Locals:
    type(ed_patch_type), pointer :: current_patch
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'thin_low_perfect() beginning.'
    
    ! Add handling for 'empty' pfts value.
    
    current_patch => site%oldest_patch
    do while (associated(current_patch))
      
      ! Thinning is applied to all of each patch by omitting the patch_fraction argument:
      call thin_patch_low_perfect(patch = current_patch, pfts = pfts, thin_fraction = thin_fraction)
      
      current_patch => current_patch%younger
    end do ! Patch loop.
    
  end subroutine thin_low_perfect

  !=================================================================================================

  subroutine thin_patch_low_probabilistic(patch, pfts, thin_fraction)!Name? thin_low_realistic  thin_low_X_patch
    ! ----------------------------------------------------------------------------------------------
    ! Perform a low thinning on a patch such that most thinned trees are those with smaller
    ! diameters but some trees are also removed from the larger size classes.  This is more
    ! realistic than a 'perfect' low thinning since foresters generally have to balence spacing with
    ! size during thinning.
    ! 
    ! The Thinning Model:
    !   A model of thinning was constructed based on observations and some basic principles and is
    ! able to reproduce a realistic demographic pattern.  The probability of thinning at a given DBH
    ! is modeled as a function of the fraction of trees to be thinned.  The probability density
    ! function is a logistic equation of the form:
    !
    ! 1 / (1 + e^-k(DBH_trans - DBH_0))
    !
    ! Where k is the steepness parameter, DBH_0 is the midpoint parameter, and DBH_trans is
    ! transformed diameter at breast height.
    !
    ! Steepness:
    !   The steepness parameter controls how 'imperfect' the low thinning is. We use a fixed average
    ! value derived from data. However, the steepness can be adjusted to capture a wider range of
    ! thinning behavior.  At a steepness of zero the thinning is proportional.  At the other extreem
    ! at a very high steepness the behavior approaches a 'perfect' low thinning.  this routine
    ! should not be used for these extremes as such thinnings can be simulated with more efficiently
    ! with other calls.  To allow the behavior to be tuned to different systems that may differ from
    ! that to derive the default value the parameter could be made an optional parameter.
    !
    ! Midpoint:
    !   The midpoint or inflection point of the logistic curve occurs at the DBH where the chance of
    ! being thinned is 50%.  Moving the curve to the right increases the number of trees that will
    ! be thinned.  We derived a relationship between the midpoint value and the intensity of
    ! thinning from a observations.
    !
    ! DBH Transformation:
    !   To apply a simple model to forests of different ages and demographics the thinning model
    ! uses DBHs that are transformed so that mean of the distribution is centered at 0.
    !
    ! Parameterization:
    !   The model was parameterized using a set of thinning trial experiments in loblolly pine
    ! plantations across the Southerestern U.S. (manuscript in progress).
    !
    !   This model predicts the number of trees to be thinning quite accurately for observations.
    ! Using the results of the calculation from the initial calculated midpoint value may do a
    ! reasonable job for some purposes.  However, it is important for comparably that the thinning
    ! occur very close to the rate specified.  To ensure this we iterate the calculation modifying
    ! the midpoint value until the requested thinning fraction is achieved.
    !
    ! Realism and Weaknesses:
    !   When performing a low thinning (thinning from below) foresters preferentially remove the
    ! smaller trees but this in not perfect for several reasons.  First, trees are rarely evenly
    ! spaced with regard to their size.  If small trees are clustered removing them will create
    ! large canopy gaps and other trees will be closely spaced.  In such cases some larger trees may
    ! be removed to maintain optimal spacing.  Also larger trees may be removed due to damage or
    ! poor form.  And of couse some stochastic errors are to be expected.
    !   This routine captures this demographic effect but is unrealistic in its smoothness.  The
    ! smooth PDF used to calculate thinning results in some thinning at all sizes, where in reality
    ! some small and many large trees escape thinning completely.
    !   This algorithm was designed for single aged stands.  For stands that have a lot smaller
    ! trees result may not be as intended.  Since canopy trees are the focus of thinning it is
    ! probably appropriate in such cases to have a size cut-off under which all trees are either
    ! removed or ignored.
    !   In reality some very small trees may be left on site.
    !
    ! Notes:
    ! - This routine does use effective values in its calculations!
    ! - Instead of looping through all the cohorts multiple times we could make a list of the valid cohorts once and loop through those...
    !
    ! Use fraction or specific number of trees?????
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_patch_type), intent(in), target :: patch
    integer(i4), dimension(:), intent(in) :: pfts ! An array of PFT IDs to thin.
    real(r8), intent(in) :: thin_fraction ! Could compute from a number of tree to thin?
    
    ! Locals:
    real(r8), parameter :: midpoint_slope = 16.48_r8 ! Exp. 7/58
    real(r8), parameter :: midpoint_intercept = -8.34_r8 ! Exp. 7/58
    real(r8), parameter :: model_steepness = -0.4428_r8 ! -0.442796647865622 Exp. 7/58
    real(r8), parameter :: thin_tolerance = 0.1_r8 ! Thin trees to with 0.1 trees (n) of the goal.
    
    real(r8) :: sum_dbh ! Accumulator
    real(r8) :: num_trees ! Accumulator
    real(r8) :: mean_dbh ! Mean (weighted) DBH of the patch.
    
    real(r8) :: model_midpoint ! In centered DBH space.
    real(r8) :: midpoint_step ! Step for iteratively solving the thinning weights.
    real(r8) :: thin_goal ! The number of trees that should be thinned.
    real(r8) :: thin_total ! Accumulator for thinning calculation
    real(r8) :: thin_last ! The number of trees to thin calculated by the last iteration. 
    real(r8) :: thin_prob ! The probability of thinning for a cohort.
    real(r8) :: dbh_tranformed ! Cohort DBH adjusted so the mean DBH is at 0.
    
    integer :: cycles ! Loop counter for debugging purposes.
    
    type(ed_cohort_type), pointer :: current_cohort
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'thin_patch_low_probabilistic() entering.'
    
    sum_dbh = 0.0_r8
    num_trees = 0.0_r8
    thin_total = 0.0_r8
    thin_last = 0.0_r8
    !midpoint_step = 0.05_r8 ! The initial step size is arbitrary. We assume we are pretty close when we start.
    midpoint_step = 0.5_r8 ! Start by adjusting the midpoint by 0.5 cm if we are not within the tolerance.
    cycles = 0
    
    ! Validity checking for thin_fraction!!!!!
    
    ! Calculate the mean DBH of the trees (cohort DBH weighted by number) to be thinned (exclude some?)
    current_cohort => patch%shortest
    do while(associated(current_cohort))
      if (any(pfts == current_cohort%pft)) then
        ! Should use patch_effective_n() here!!!!!
        sum_dbh = sum_dbh + (current_cohort%dbh * current_cohort%n)
        num_trees = num_trees + current_cohort%n
      endif
      current_cohort => current_cohort%taller
    end do ! Cohort loop.
    mean_dbh = sum_dbh / num_trees
    
    ! Calculate the number of trees to be removed based on the thinning fraction:
    thin_goal = num_trees * thin_fraction
    
    ! Calculate the midpoint parameter based on the fraction of trees to be thinned:
    model_midpoint = midpoint_slope * thin_fraction + midpoint_intercept
    
    if (debug) write(fates_log(), *) 'Initial model_midpoint:', model_midpoint ! Temporary reporting.
    
    ! Starting with the initial midpoint value calculate thinning weights and repeat the process
    ! with adjusted values until we are within the tolerance:
    do while(abs(thin_total - thin_goal) > thin_tolerance)
      cycles = cycles + 1 ! Counter ...
      thin_total = 0.0_r8
      
      ! For each valid cohort determine the probability of thinning and number of trees to thin:
      current_cohort => patch%shortest
      do while(associated(current_cohort))
        
        if (any(pfts == current_cohort%pft)) then
          dbh_tranformed = current_cohort%dbh - mean_dbh
          !thin_prob = 1 / (1 + exp(-1.0r8 * model_steepness * (dbh_tranformed - model_midpoint))) ! -model_steepness?
          thin_prob = 1 / (1 + exp(-model_steepness * (dbh_tranformed - model_midpoint)))
          thin_total = thin_total + (thin_prob * current_cohort%n)
        endif
        
        current_cohort => current_cohort%taller
      end do ! Cohort loop.
      
      ! Adjust the midpoint for the next cycle (may not be used if we are close enough already):
      if (cycles > 1) then
        ! If we passed over the goal since the last step decrease the step size:
        if ((thin_last < thin_goal .and. thin_total > thin_goal) .or. &
            (thin_last > thin_goal .and. thin_total < thin_goal)) then
          midpoint_step = midpoint_step / 2.0_r8
        endif
      endif
      
      ! Step up or down:
      if (thin_total < thin_goal) then
        model_midpoint = model_midpoint + midpoint_step
      elseif (thin_total > thin_goal) then
        model_midpoint = model_midpoint - midpoint_step
      endif
      
      ! Are there invalid midpoint values that will let us know the search has gone wrong?
      
      thin_last = thin_total ! Record the trees thinned for comparison during the next cycle.
      ! cycles = cycles + 1 ! Counter for debugging purposes only.
      
      if (debug) then ! Temporary reporting:
        write(fates_log(), *) 'cycles:', cycles
        write(fates_log(), *) 'model_midpoint:', model_midpoint
        write(fates_log(), *) 'midpoint_step:', midpoint_step
        write(fates_log(), *) 'thin_total:', thin_total
      endif
      
    end do ! Thinning calculation
    
    if (debug) write(fates_log(), *) 'Thinning calculation completed.' ! Temporary reporting.
    
    ! Once the proper thinning weights have been solved apply the thinning mortality:
    ! This is a bit repetetive but there is little avoiding that this must be done in a loop.
    current_cohort => patch%shortest
    do while(associated(current_cohort))
      if (debug) write(fates_log(), *) 'Starting cohort.' ! Temporary reporting.
      
      if (any(pfts == current_cohort%pft)) then
        dbh_tranformed = current_cohort%dbh - mean_dbh
        thin_prob = 1 / (1 + exp(-model_steepness * (dbh_tranformed - model_midpoint)))
        
        if (debug) write(fates_log(), *) 'thin_prob:', thin_prob! Temporary reporting.
        
        ! Some or all trees could be left in place but we assume a commercial harvest:
        call kill(cohort = current_cohort, flux_profile = bole_harvest, kill_fraction = thin_prob)
      endif
      
      current_cohort => current_cohort%taller
    end do ! Cohort loop.
    
    if (debug) then
      ! Report algorithm stats:
      write(fates_log(), *) 'Thinning fraction specified: ', thin_fraction
      write(fates_log(), *) 'As trees:                    ', thin_goal
      write(fates_log(), *) 'Number of trees thinned:     ', thin_total
      ! Starting midpoint?
      write(fates_log(), *) 'Solution midpoint:          ' , model_midpoint
      write(fates_log(), *) 'Cycles to solve:             ', cycles
      write(fates_log(), *) 'thin_patch_low_probabilistic() exiting.'
    endif
  end subroutine thin_patch_low_probabilistic

  !=================================================================================================

  subroutine thin_low_probabilistic(site, pfts, thin_fraction) ! where = everywhere
    ! ----------------------------------------------------------------------------------------------
    ! Perform a low thinning for a site such that most thinned trees are those with smaller
    ! diameters but some trees are also removed from the larger size classes.
    !
    ! In it's current form this routine just applies thin_patch_low_probabilistic() to all patches.
    ! In the future we will (may) allow targeting of thinning behavior to specific patches via the
    ! 'where' parameter.
    !
    ! VM Event Interface thin_low_probabilistic([pfts], thin_fraction)
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_site_type), intent(in), target :: site ! The current site object.
    integer(i4), dimension(:), intent(in) :: pfts ! An array of PFT IDs to thin.
    real(r8), intent(in) :: thin_fraction ! Fraction of trees to thin.
    
    ! Locals:
    type(ed_patch_type), pointer :: current_patch
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'thin_low_probabilistic() beginning.'
    
    ! Add handling for 'empty' pfts value.
    
    current_patch => site%oldest_patch
    do while (associated(current_patch))
      
      call thin_patch_low_probabilistic(current_patch, pfts, thin_fraction)
      
      current_patch => current_patch%younger
    end do ! Patch loop.
    
  end subroutine thin_low_probabilistic

  !=================================================================================================
  ! Harvest Subroutines:
  ! ...
  !
  ! We may want to provide a non-disturbing harvest routines for continuous cover forestry and
  ! agriculture.
  !
  !=================================================================================================

  ! Patch level: Clear-cut, non-marketable trees are cut and left on site...
  ! Site level: Preferred species, preferred land type...
  ! What is the correct most generic harvest function?????

  !subroutine harvest(patch, pfts, dbh_min, dbh_max, ht_min, ht_max, fraction) ! log()??????
    ! ----------------------------------------------------------------------------------------------
    ! 
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    
    ! Arguments:
    
    ! Locals:
    
    ! ----------------------------------------------------------------------------------------------
    
    ! This may end up being a convenience wrapper of a patch level kill function.
    
  !end subroutine harvest

  !=================================================================================================

  ! Better name? harvest_by_mass()?  This should be moved to the site level functions!
  subroutine harvest_mass_min_area(site_in, harvest_c_primary, harvest_c_secondary, & ! REVIEW!
                                   pfts, dbh_min, dbh_max, ht_min, ht_max)
    ! ----------------------------------------------------------------------------------------------
    ! Harvest the requested harvest amounts from the site while attempting to use the minimum area.
    ! The harvest amounts are specified as product biomass in carbon (g or kg??????).
    !   [I need to confirm that LUH data uses wood product...]
    ! Primary and secondary lands are harvested separately.  A shortfall in one will not be made up
    ! from the other.
    ! Optional criteria specify the types and sizes of plants that qualify to be harvested.  If not
    ! specified all PFTs and sizes will be harvested.
    ! 
    ! An estimate of the harvest amount (the amount staged prior to maximum disturbance calculation)
    ! is also calculated and returned in the harvest amount arguments.
    !
    ! Heuristic: (using a logging metaphor)
    !   Let's assume that foresters will seek out stands that have the highest yield of marketable
    ! timber per unit area and exploit them first, before moving on the next best locations.
    !   Let's assume that a patch = a stand, or a set of similar stands that will be exploited as a
    ! group (and to the same extent?????).
    !   We need a measure of site / stand quality so we can figure out where to start.
    ! (There must certainly be methods from natural resource economics that could be used here but
    ! I'm not familiar enough with them, so I'm using my own approximation.)
    ! We will use the number of stems / ha that meet the harvest criteria as the metric of stand
    ! quality and rank patches according to this.
    ! 
    !   While we allow the mix of PFTs that will be harvested to be specified we do not priorize the
    ! harvest of any listed PFT over another.  We also don't harvest the biggest trees first.  We
    ! assume that any tree meeting the harvest criteria is as good as any other.  Therefore, all
    ! PFTs and sizes will be harvested proportionally according to abundance.  This is consist with
    ! the idea of a patch representing a heterogeneous and well mixed population.
    !
    !   Also, in a departure from the historic logging module, we don't assume that harvest can only
    ! come from the top layer of the canopy.  We change this assumption because the top canopy layer
    ! is not magic in FATES.  If a patch is mostly trees there may be many cohorts of appropriate
    ! sizer for harvest in the second layer.  While FATES terms these plants 'understory' they might
    ! they may include codominant or suppressed trees.
    !
    !   While we use a logging metaphor above, this subroutine could also be used for other purposes.
    ! Harvest of fodder would be a good application. To fully realize this a harvest mode parameter
    ! needs to be added here and in harvestable_biomass().
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use FatesConstantsMod, only : primaryforest, secondaryforest
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: site_in
    ! Amounts to harvest in product carbon by land class, the actual amount harvested is returned:
    ! Note: The internal units of the model are kgC/site and these are currently in kg.  These
    ! inputs however should probably be in g/m^2.  Need to convert.
    real(r8), intent(inout) :: harvest_c_primary
    real(r8), intent(inout) :: harvest_c_secondary
    ! An array of PFT IDs to include in the basal area calculation: Make optional?
    integer(i4), dimension(:), intent(in) :: pfts
    
    ! Size range of to trees to harvest.  Defaults to everything, otherwise some range of sizes,
    ! e.g. > 10cm DBH, < 15 m in height, etc.
    real(r8), intent(in), optional :: dbh_min
    real(r8), intent(in), optional :: dbh_max
    real(r8), intent(in), optional :: ht_min
    real(r8), intent(in), optional :: ht_max
    
    ! Locals:
    type(ed_patch_type), pointer :: current_patch
    type(ed_cohort_type), pointer :: current_cohort
    type(ed_patch_type), pointer :: best_patch
    
    real(r8) :: the_dbh_min
    real(r8) :: the_dbh_max
    real(r8) :: the_ht_min
    real(r8) :: the_ht_max
    
    real(r8) :: patch_harvestable_stems ! Number of plants available for harvest.
    real(r8) :: cohort_harvest ! cohort_harvestable_biomass
    real(r8) :: patch_harvestable_biomass
    real(r8) :: best_patch_harvestable_stems
    real(r8) :: best_patch_harvestable_biomass
    real(r8) :: patch_harvest_fraction
    real(r8) :: harvest_remaining
    real(r8) :: harvest_total ! Running total of wood product harvested.
    
    integer :: forest_class ! Primary or secondary
    integer :: num_patches ! The number of patches of a certain forest class.
    integer :: search_cycles ! Counter
    integer :: i ! Counter
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'harvest_mass_min_area() entering.'
    
    ! An estimate of the harvest flux will be made and stored: ?????????????????
    !site_in%harvest_carbon_flux = 0.0_r8
    
    ! Check and set the size specifications:
    call validate_size_specifications(the_dbh_min, the_dbh_max, the_ht_min, the_ht_max, &
                                      dbh_min, dbh_max, ht_min, ht_max)
    
    ! Temporary reporting:
    if (debug) then
      write(fates_log(), *) 'harvest_mass_min_area(): the_dbh_min = ', the_dbh_min
      write(fates_log(), *) 'harvest_mass_min_area(): the_dbh_max = ', the_dbh_max
      write(fates_log(), *) 'harvest_mass_min_area(): the_ht_min = ', the_ht_min
      write(fates_log(), *) 'harvest_mass_min_area(): the_ht_max = ', the_ht_max
    end if
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Repeat for primary and secondary land:
    do i = 1, 2
      if (i == 1) then
        if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Starting primary forest.'
        
        harvest_remaining = harvest_c_primary
        forest_class = primaryforest
        !harvest_c_primary = 0 ! Repurpose as return value.
      else
        if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Starting secondary forest.'
        
        harvest_remaining = harvest_c_secondary
        forest_class = secondaryforest
        !harvest_c_secondary = 0 ! Repurpose as return value.
      endif
      
      harvest_total = 0.0_r8
      
      ! Count patches in this land class:
      ! This is probably wrong as we increment the counter below for all patches!
      num_patches = 0
      current_patch => site_in%oldest_patch
      do while (associated(current_patch) .and. &
                current_patch%anthro_disturbance_label == forest_class)
        num_patches = num_patches + 1
        
        current_patch => current_patch%younger
      end do ! Patch loop
      
      if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Through counting.' ! Temp
      if (debug) write(fates_log(), *) 'harvest_mass_min_area(): num_patches = ', num_patches ! Temp
      
      ! Find the patch of the highest quality and harvest from it, then find the next best patch,
      ! and so on iteratively until either we harvest enough or we run out of patches:
      
      ! It would be most efficient to discard / skip the previous best patch after each search and
      ! harvest cycle.  However the fact that the patches are organized as a linked list makes this
      ! difficult.  Since the appropriate trees are removed from patches as they are harvested they
      ! won't match again as the best (assuming we use the 'effective' accessors).
      ! A better way to do this would be to make our own list of patches.
      !
      ! For safety we make sure that we stop after we have done as many cycle as there are patches.
      ! However, this logic is probably superflous since the loop will stop if we harvest enough or
      ! run out of patches with appropriate PFTs.
      
      search_cycles = 0 ! Change to 1?
      !do while (harvest_remaining > 0 .and. search_cycles < num_patches) !Tiny????
      do while (search_cycles < num_patches)
        best_patch_harvestable_stems = 0.0_r8
        best_patch => null()
        
        if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Search loop.' ! Temp
        if (debug) write(fates_log(), *) 'harvest_mass_min_area(): search_cycles = ', search_cycles ! Temp
        
        ! Walk the patches from oldest to youngest, on the assumption that older patches will likely
        ! have larger and more valuable trees:
        current_patch => site_in%oldest_patch
        do while (associated(current_patch))
          patch_harvestable_stems = 0.0_r8
          
          if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Patch loop.' ! Temp
          
          if (current_patch%anthro_disturbance_label == forest_class) then
            ! Calculate the number of trees that meet the harvest criteria for the patch:
            ! Since n is per nominal hectare the count = stem density, which is what we want to compare.
            current_cohort => current_patch%shortest
            do while(associated(current_cohort))
              
              if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Cohort loop.' ! Temp
              
              ! Temporary reporting:
              if (debug) then
                write(fates_log(), *) 'harvest_mass_min_area(): cohort DBH = ', current_cohort%dbh
                write(fates_log(), *) 'harvest_mass_min_area(): cohort height = ', current_cohort%hite
              end if
              
              ! Note: The DBH and height ranges have be initialized such that those conditions not
              ! actually used will pass the following conditionals:
              if (any(pfts == current_cohort%pft) .and. &
                  current_cohort%dbh >= the_dbh_min .and. current_cohort%dbh <= the_dbh_max .and. &
                  current_cohort%hite >= the_ht_min .and. current_cohort%hite <= the_ht_max) then
                
                if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Valid cohort.' ! Temp
                
                patch_harvestable_stems = patch_harvestable_stems + effective_n(current_cohort)
                
                ! Accumulated the potentially harvestable biomass:
                ! Currently this the bole X logging module modifiers
                
                ! Get the harvestable mass for the entire cohort:
                cohort_harvest = cohort_harvestable_biomass(cohort = current_cohort, &
                                                            harvest_profile = bole_harvest, &
                                                            staged = .true.)
                patch_harvestable_biomass = patch_harvestable_biomass + cohort_harvest
              endif
              
              current_cohort => current_cohort%taller
            end do ! Cohort loop.
          endif ! (current_patch%anthro_disturbance_label == forest_class)
          
          if (debug) write(fates_log(), *) 'harvest_mass_min_area(): patch_harvestable_biomass: ', patch_harvestable_biomass
          
          ! If this is the best patch yet, store it:
          ! It is possible that no patch will have harvestable trees.
          if (patch_harvestable_stems > best_patch_harvestable_stems) then
            best_patch => current_patch
            best_patch_harvestable_stems = patch_harvestable_stems
            best_patch_harvestable_biomass = patch_harvestable_biomass
          endif
          
          current_patch => current_patch%younger
        end do ! Patch loop
        
        ! Harvest from the best available patch:
        if (best_patch_harvestable_stems > 0.0_r8) then ! Tiny?????? Add tolerance....
          
          if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Harvesting from best patch.'
          
          ! Note: Similar to logic in thinning:
          if (harvest_remaining >= best_patch_harvestable_biomass) then
            ! If less than the remaining demand harvest all the trees in the patch:
            call kill_patch(patch = best_patch, flux_profile = bole_harvest, pfts = pfts, &
                            dbh_min = the_dbh_min, dbh_max = the_dbh_max, &
                            ht_min = the_ht_min, ht_max = the_ht_max, &
                            kill_fraction = 1.0_r8, area_fraction = 1.0_r8) ! kill_fraction is not currently optional but should be.
            
            harvest_total = harvest_total + best_patch_harvestable_biomass ! Record the harvest.
            
            harvest_remaining = harvest_remaining - best_patch_harvestable_biomass ! Is this needed with harvest_total?  Move down?
            ! It might be better to store the requested amount and subtract harvest_total each time?
            
          else
            ! Otherwise harvest a fraction across all sizes.  This is equivalent to harvesting the
            ! smallest area of a heterogeneous area represented by the patch:
            
            if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Harvesting from best patch (partial).'
            
            patch_harvest_fraction = harvest_remaining / best_patch_harvestable_biomass
            
            ! Temporary reporting:
            if (debug) then
              write(fates_log(), *) 'harvest_mass_min_area(): harvest_remaining = ', harvest_remaining
              write(fates_log(), *) 'harvest_mass_min_area(): best_patch_harvestable_biomass = ', best_patch_harvestable_biomass
              write(fates_log(), *) 'harvest_mass_min_area(): patch_harvest_fraction = ', patch_harvest_fraction
            end if
            
            ! Remove the wood in a clear-cut like fashion. This results in the smallest distubance:
            ! It could be removed in another manner.  Think about this.
            call kill_patch(patch = best_patch, flux_profile = bole_harvest, pfts = pfts, &
                            dbh_min = the_dbh_min, dbh_max = the_dbh_max, &
                            ht_min = the_ht_min, ht_max = the_ht_max, &
                            kill_fraction = patch_harvest_fraction, &
                            area_fraction = patch_harvest_fraction)
            
            harvest_total = harvest_total + (best_patch_harvestable_biomass * patch_harvest_fraction)
            
            ! We have matched the harvest amount so stop:
            exit
            ! We could include the conditional (harvest_remaining > 0) in the bounding do loop but
            ! using exit avoids any potential floating point comparison issues.
          endif ! (harvest_remaining > 0)
        else
          ! If we have run out of patches that have harvestable trees terminate early:
          if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Out of patches.'
          exit
        endif ! (best_patch_harvestable_stems > 0)
        
        search_cycles = search_cycles + 1
      end do ! (search_cycles < num_patches)
      
      ! Return the actual amounts harvested:
      if (i == 1) then
        harvest_c_primary = harvest_total ! Repurpose as return value.
      else
        harvest_c_secondary = harvest_total ! Repurpose as return value.
      endif
      
      if (debug) write(fates_log(), *) 'harvest_mass_min_area(): Harvested: ', harvest_total
      
    end do ! Primary and secondary loop.
    
    if (debug) write(fates_log(), *) 'harvest_mass_min_area() exiting.'
  end subroutine harvest_mass_min_area

  !=================================================================================================

  function plant_harvestable_biomass(cohort, harvest_profile) result(harvest)  ! REVIEW!
    ! ----------------------------------------------------------------------------------------------
    ! Return the harvestable biomass per plant for a given cohort.
    ! 
    ! This function compartmentalizes this calculation for a level of abstraction, to reduce
    ! dependancies and code repetion.  The calculation used historically for the logging module
    ! is complex and depends on some assumptions that might change in the future.
    !
    ! This function is also a step towards having more that one type of harvest modality.
    ! The logging module use a bole harvest.  Whole tree harvest might be more appropriate for
    ! biomass extraction etc.  We may also add biomass calculations for non-trees in the future and
    ! this provided a centralized interface.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use EDParamsMod, only : logging_export_frac
    use FatesLitterMod, only : ncwd ! The number of coarse woody debris pools.
    use SFParamsMod, only : SF_VAL_CWD_FRAC
    use PRTGenericMod, only : all_carbon_elements
    use PRTGenericMod, only : struct_organ, sapw_organ
    
    ! Arguments:
    type(ed_cohort_type), intent(in), target :: cohort
    integer(i4), intent(in), optional :: harvest_profile ! The type of harvest to perform.
                                                         ! Valid values are any harvest type flux
                                                         ! profile).  Defaults to bole_harvest if
                                                         ! not provided.
    
    ! Locals:
    real(r8) :: harvest ! Return value.
    integer(i4) :: the_profile ! Name?????
    integer(i4), pointer, dimension(:) :: valid_pfts
    
    ! ----------------------------------------------------------------------------------------------
    
    harvest = 0.0_r8
    
    ! Determine the type of harvest:
    if (present(harvest_profile)) then
      the_profile = harvest_profile
    else
      the_profile = bole_harvest
    endif
    
    ! For the implemented modes the only difference is what is considered a harvestable 'tree':
    select case(the_profile)
      case(logging_traditional)
        valid_pfts => woody_pfts
        
      case(bole_harvest)
        valid_pfts => tree_pfts
        
      case default
        write(fates_log(),*) 'Invalid harvest profile passed in.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
    end select
    
    if (any(valid_pfts == cohort%pft)) then
       ! For bole harvest we maintain the historic logging module calculation:
       !   (Originally found in EDPatchDynamicsMod: disturbance_rates().)
       ! bole harvest = trunk wood * export fraction * Spitfire Coarse Woody Debris fraction
       ! Where:
       ! harvest = (sapwood + structual wood) * aboveground biomass fraction
       ! export fraction =  The fraction of harvested wood that actually makes it off the site.
       !                    Is this completely open to interpretation?
       !                    Not sure if this includes the transportation carbon losses?????
       ! SF_val_CWD_frac = Talk to Jakie about this spitfire parameter!
       harvest = (cohort%prt%GetState(sapw_organ, all_carbon_elements) + &
                  cohort%prt%GetState(struct_organ, all_carbon_elements)) * &
                 prt_params%allom_agb_frac(cohort%pft) * &
                 SF_val_CWD_frac(ncwd) * &
                 logging_export_frac
      
      ! Add other options here:
    endif
    ! Otherwise the cohort has 0 harvestable biomass, which is a valid result.
    
  end function plant_harvestable_biomass

  !=================================================================================================

  function cohort_harvestable_biomass(cohort, harvest_profile, staged) result(harvest) ! Review!
    ! ----------------------------------------------------------------------------------------------
    ! Return the harvestable biomass for the whole cohort or, if staged = true, the yield based on
    ! the harvest rates stored (staged) in the cohort object passed in.
    !
    ! The following is not future safe because it does not adjust n for any staged mortalities that
    ! would occur before harvest.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_cohort_type), intent(in), target :: cohort
    ! The type of harvest to perform (harvest activity / flux profile).  Defaults to bole_harvest if not provided.  Name?
    integer(i4), intent(in), optional :: harvest_profile ! The type of harvest to perform.
                                                         ! Valid values are any harvest type flux
                                                         ! profile).  Defaults to bole_harvest if
                                                         ! not provided.
    logical, intent(in), optional :: staged              ! If true return the harvest yield given
                                                         ! the staged harvest rates.  Otherwise
                                                         ! return the yield if the whole cohort was
                                                         ! harvested. Defaults to false.
    
    ! Locals:
    real(r8) :: harvest ! Return value.
    type(ed_cohort_type), pointer :: current_cohort
    real(r8) :: harvest_fraction
    integer :: the_profile
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Determine the type of harvest:
    if (present(harvest_profile)) then
      the_profile = harvest_profile
      ! Validity checking is performed in plant_harvestable_biomass().
    else
      the_profile = bole_harvest
    endif
    
    ! effective_n() won't work here because the staged harvest mortalities will be included in the
    ! the calculation that reduces n.
    harvest = plant_harvestable_biomass(cohort, harvest_profile) * cohort%n
    
    if (present(staged) .and. staged == .true.) then ! If staged is not passed in default to false.
      ! Get the harvest rate (assumes only one, which should be safe):
      ! Also assumes only two potential harvest types which is probably not safe!!!!!
      harvest_fraction = max(cohort%lmort_direct, cohort%vm_mort_bole_harvest)
      
!       if (harvest_fraction == 0.0_r8) then
!         write(fates_log(),*) 'No harvest is currently staged.'
!         call endrun(msg = errMsg(__FILE__, __LINE__))
!       endif
      
      if (harvest_fraction /= 0.0_r8) then
        ! Reduce the harvest amount to reflect the staged mortally:
        harvest = harvest * (1.0_r8 - harvest_fraction)
      endif
    endif
    
  end function cohort_harvestable_biomass

  !=================================================================================================

  subroutine clearcut_patch(patch, pfts, dbh_min, ht_min, patch_fraction)
    ! ----------------------------------------------------------------------------------------------
    ! Perform a clearcut harvest on a patch.
    !   All the trees matching the PFT and size specification are bole harvested and moved to the
    ! wood product pool.  Other trees and understory are killed and left on site. For simplicity we
    ! assume that there will be no maximum tree size limit for harvest, though there may be minimum
    ! size, which can be specified by either diameter or height.
    !
    ! It could be more realistic to leave some harvestable biomass on site as a harvest
    ! inefficiency (the logging module behavior does this).  Also leaving some understory alive
    ! would be appropriate if a reasonable fraction could be determined.
    !
    ! This routine uses kill() calls that assume this is the first disturbing event that has
    ! occurred.  It should be revised used 'disturbed' calls.
    !
    !-----------------------------------------------------------------------------------------------
    ! Harvest Operation: Patch Level
    ! Operation Type: Harvest
    ! Operation Level: Patch
    ! Regime: (single age) clearcut harvest rotation
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_patch_type), intent(in), target :: patch
    ! Which tree PFTs are harvestable trees?  If omitted all woody trees will be used:
    !integer(i4), dimension(:), intent(in) :: pfts ! Defaults to trees?
    integer(i4), intent(in), optional, target :: pfts(:)
    
    ! Size specification:
    ! Minimum size of to trees to harvest.  Only provide a minimum DBH or height, not both:
    real(r8), intent(in), optional :: dbh_min ! Defaults to 0.
    real(r8), intent(in), optional :: ht_min ! Defaults to 0.
    
    ! The fraction of the patch that the mortality is removed from / the size of the disturbance:
    real(r8), intent(in), optional :: patch_fraction ! Defaults to 1 (whole patch).
    
    ! Locals:
    integer(i4), pointer, dimension(:) :: the_pfts ! Holds the PFTs to harvest, computed from arguments.
    real(r8) :: the_dbh_min
    real(r8) :: the_ht_min
    real(r8) :: the_patch_fraction
    
    integer :: i, j ! Counters
    integer:: all_pfts(12) = (/ (integer :: k, k = 1, 12) /) ! Generalize and make global!!!!!
    integer:: other_pfts(12)
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'clearcut_patch() entering.'
    
    ! Validate the PFTs:
    !if (present(pfts) .and. (pfts /= vm_empty_array)) then
    if (present(pfts) .and. all(pfts /= vm_empty_integer)) then
      ! Confirm PFTs to harvest are all trees:
      do i = 1, size(pfts)
        if (.not. any(pfts(i) == tree_pfts)) then
          write(fates_log(),*) 'clearcut_patch(): Only tree PFTs are expected.'
          write(fates_log(),*) 'Tree PFTs =    ', tree_pfts
          write(fates_log(),*) 'Selected PFTs =', pfts
          call endrun(msg = errMsg(__FILE__, __LINE__)) ! We could just warn here?
        endif
      end do
      
      the_pfts => pfts
    else
      the_pfts => tree_pfts
    endif ! present(dbh)...
    
    ! Validate the size specifications:
    call validate_size_specifications(dbh_min_out = the_dbh_min, ht_min_out = the_ht_min, &
                                      dbh_min = dbh_min, ht_min = ht_min)
    
    ! Validate the patch fraction:
    if (present(patch_fraction)) then
      the_patch_fraction = patch_fraction
    else
      the_patch_fraction = 1.0_r8
    endif
    
    if (the_patch_fraction <= 0.0_r8 .or. the_patch_fraction > 1.0_r8) then
      write(fates_log(),*) 'Invalid value for the_patch_fraction argument.', the_patch_fraction
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    ! Harvest the appropriate trees:
    call kill_patch(patch = patch, flux_profile = bole_harvest, pfts = pfts, &
                    !dbh_min = the_dbh_min, dbh_max = the_dbh_max, &
                    !ht_min = the_ht_min, ht_max = the_ht_max, &
                    dbh_min = the_dbh_min, ht_min = the_ht_min, &
                    kill_fraction = 1.0_r8, area_fraction = 1.0_r8)
    
    ! Kill everything else in place:
    
    ! If there were size limits set for the harvested PFTs then there may be some trees remaining.
    ! Assume those were killed in the process of harvest but were left on site:
    if (present(dbh_min)) then
      call kill_patch(patch = patch, flux_profile = in_place, pfts = pfts, dbh_max = the_dbh_min, &
                      kill_fraction = 1.0_r8, area_fraction = 1.0_r8)
      
    else if (present(ht_min)) then
      call kill_patch(patch = patch, flux_profile = in_place, pfts = pfts, ht_max = the_ht_min, &
                      kill_fraction = 1.0_r8, area_fraction = 1.0_r8)
    endif
    
    ! Kill the unharvested PFTs and leave them on site:
    ! It might be more realistic to leave some fraction of plants, especially grasses, alive.
    
    ! Get the reciprocal PFTs:
    j = 0
    do i = 1, size(all_pfts)
      if (.not. any(all_pfts(i) == pfts)) then
        j = j + 1
        other_pfts(j) = all_pfts(i)
      endif
    end do
    
    if (j > 0) then ! If there are any remaining PFTs kill them:
      call kill_patch(patch = patch, flux_profile = in_place, pfts = other_pfts(1:j), &
                      kill_fraction = 1.0_r8, area_fraction = 1.0_r8)
    endif
    
    if (debug) write(fates_log(), *) 'clearcut_patch() exiting.'
  end subroutine clearcut_patch

  !=================================================================================================

  subroutine clearcut(site, pfts, dbh_min, ht_min) ! patch_fraction
    ! ----------------------------------------------------------------------------------------------
    ! Perform a clearcut harvest for all patches of a site.
    !
    ! In the future we will (may) allow targeting of thinning behavior to specific patches via the
    ! 'where' parameter.
    !
    ! VM Event Interface: clearcut(pfts, [dbh_min], [ht_min]) pfts optional?
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    
    type(ed_site_type), intent(in), target :: site ! The current site object.
    integer(i4), dimension(:), intent(in), optional :: pfts ! An array of PFT IDs to harvest.
    
    ! Size specification:
    ! Minimum size of to trees to harvest.  Only provide a minimum DBH or height, not both:
    real(r8), intent(in), optional :: dbh_min ! Defaults to 0.
    real(r8), intent(in), optional :: ht_min ! Defaults to 0.
    
    ! Locals:
    type(ed_patch_type), pointer :: current_patch
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'clearcut() beginning.'
    
    current_patch => site%oldest_patch
    do while (associated(current_patch))
      
      call clearcut_patch(current_patch, pfts, dbh_min, ht_min)
      
      current_patch => current_patch%younger
    end do ! Patch loop.
    
  end subroutine clearcut

  !=================================================================================================
  ! Utilities:
  !
  ! Effective values:
  ! At the time management activities are calculated the patch and cohort data structures are in a
  ! state of flux.  Some mortality rates may have been calculated but have not been applied, which
  ! means the number of plants, and other dependent metrics, do not yet reflect these mortalities.
  ! This is due to the way the event sequence is structured and because only one disturbance type
  ! can currently occur in a patch per time step.  Rates are calculated and stored at one point in
  ! the event loop and the plant numbers are updated later. Some rates that are calculated and
  ! 'staged' many be overridden during disturbance comparison and ultimately not occur.
  !
  ! This complicates the calculation of management activities that depend on stand properties in a
  ! couple ways.
  ! 1. Iterative management decisions:
  !   Some management actives may require an iterative or step by step approach.  This requires that
  ! we be able to make a change to an ecosystem and see the result (within the same time step) before
  ! proceeding.
  ! 
  ! 2. We would like to be able to apply multiple management activities at the same time step in
  ! some cases in order to be able to simulate dynamic and heterogeneous sites / grid cells.
  ! Not all management activities are compatible with each other but some are.  We really need to
  ! identify where conflicts may occur and come up with logic and rules to prevent them.  We also
  ! need to be careful about the order that activities occur in.
  !
  ! 3.  There is no reason fundamental reason that management should be incompatible with
  ! simultaneous natural mortality.  Natural mortality should be happening alongside management.
  ! Even if we only allow for one dominant disturbance type it should be possible to apply low level
  ! natural mortality on top of management mortality.
  !
  !   In order to be able to combine these different sources of mortality we need to be able to get
  ! state numbers that reflect any mortality, natural or managed, that has already 'occurred', i.e.
  ! been staged, previously in this timestep.
  !
  !   To help with this we provide several functions below that return stand properties for patches
  ! and cohorts assuming all staged mortalities will be actually applied.  They should be used
  ! rather than accessing cohort states directly.
  !
  !   These functions are a work in progress.  I haven't figured out quite what I need to do yet in
  ! its entirety, but by using function calls I code the higher level behavior without dependance
  ! on the final solution details.
  !   Currently we are only trying to support item 1 above.  Handling multiple mortality types or
  ! management events will be left for later.
  !
  !   Because n is expressed in terms of the nominal hectare BAI = BA and stem density = n.
  ! An interface alias is provided for the later relationship to reduce dependency on that
  ! assumption and to make for more readable code.
  !
  ! Note: Need to make kill() aware of the effective state as well!
  !
  ! Can disturbed_x() functions replace all calls to effective_x()?
  !
  !=================================================================================================

  function cohort_effective_basal_area(cohort) result(effective_basal_area) ! Currently not called.
    ! ----------------------------------------------------------------------------------------------
    ! Return the effective basal area for a cohort after applying any existing staged potential
    ! mortalities.
    !
    ! While minimal this fucntion does reduce calling code complexity and readability.  It may be
    ! inlined anyway.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_cohort_type), intent(in), target :: cohort
    
    ! Locals:
    real(r8) :: effective_basal_area ! Return value
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Only the adjusted stem count is need, the rest of the equation is unchanged: 
    effective_basal_area = pi_const * (cohort%dbh / 200.0_r8)**2.0_r8 * cohort_effective_n(cohort)! effective_n(cohort)
    
  end function cohort_effective_basal_area

  !=================================================================================================

  function patch_effective_basal_area(patch, pfts) result(effective_basal_area) ! Currently not called.
    ! ----------------------------------------------------------------------------------------------
    ! Return the effective basal area for a patch after applying any existing staged potential
    ! mortalities, for only the PFTs passed in.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_patch_type), intent(in), target :: patch ! The patch to be calculated.
    ! An array of PFT IDs to include in the basal area calculation:
    integer(i4), dimension(:), intent(in) :: pfts
    !integer(i4), dimension(:), intent(in), optional :: pfts ! Should it be optional?
    
    ! Locals:
    real(r8) :: effective_basal_area ! Return value
    type(ed_cohort_type), pointer :: current_cohort
    
    ! ----------------------------------------------------------------------------------------------
    
    effective_basal_area = 0.0_r8 ! Initialize.
    
    current_cohort => patch%shortest
      do while(associated(current_cohort))
        
        if (any(pfts == current_cohort%pft)) then
          effective_basal_area = effective_basal_area + cohort_effective_basal_area(current_cohort) ! effective_basal_area(current_cohort)
        endif
        
        current_cohort => current_cohort%taller
      end do ! Cohort loop.
    
  end function patch_effective_basal_area

  !=================================================================================================

  function cohort_effective_n(cohort) !result(effective_n)
    ! ----------------------------------------------------------------------------------------------
    ! Return the effective stem count for a cohort after applying any existing staged potential
    ! mortalities.
    !
    ! This function is currently limited and expects only to be called in the context of a iterative
    ! process involving a single vegetation management mortality type.  The ability to handle
    ! multiple mortality types may be added in the future.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use EDTypesMod, only : dump_cohort
    
    ! Arguments:
    type(ed_cohort_type), intent(in), target :: cohort
    
    ! Locals:
    real(r8) :: cohort_effective_n ! Return value
    real(r8) :: staged_mortality
    
    ! ----------------------------------------------------------------------------------------------
    
    staged_mortality = 0.0_r8 ! Initialize.
    
    ! The following are pretty conservative checks:
    
    ! This is not currently intended to be used along with logging module harvests:
    if (cohort%lmort_direct /= 0.0_r8) then
      write(fates_log(),*) 'effective_n() is not currently intended to be use in conjunction with traditional logging module events.'
      call dump_cohort(cohort)
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    if (cohort%vm_mort_bole_harvest > 0.0_r8 .and. cohort%vm_mort_in_place > 0.0_r8) then
      write(fates_log(),*) 'effective_n(): more than one management mortality type staged.'
      call dump_cohort(cohort)
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    staged_mortality = max(cohort%vm_mort_bole_harvest, cohort%vm_mort_in_place)
    if (staged_mortality == 0.0_r8) then
      cohort_effective_n = cohort%n
    else
      cohort_effective_n = cohort%n * (1 - staged_mortality)
    endif
    
  end function cohort_effective_n

  !=================================================================================================

  function patch_effective_n(patch, pfts) ! result(effective_n)
    ! ----------------------------------------------------------------------------------------------
    ! Return the effective stem count for a patch after applying any existing staged potential
    ! mortalities, for only the PFTs passed in.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_patch_type), intent(in), target :: patch ! The patch to be calculated.
    integer(i4), dimension(:), intent(in) :: pfts ! Array of PFT IDs to include in the calculation.
    
    ! Locals:
    real(r8) :: patch_effective_n ! Return value
    type(ed_cohort_type), pointer :: current_cohort
    
    ! ----------------------------------------------------------------------------------------------
    
    patch_effective_n = 0.0_r8 ! Initialize.
    
    current_cohort => patch%shortest
    do while(associated(current_cohort))
    
      if (any(pfts == current_cohort%pft)) then
        patch_effective_n = patch_effective_n + cohort_effective_n(current_cohort) ! effective_n(current_cohort)
      endif
      
      current_cohort => current_cohort%taller
    end do ! Cohort loop.
    
  end function patch_effective_n

  !=================================================================================================

  function cohort_disturbed_n(cohort) result(disturbed_n)
    ! ----------------------------------------------------------------------------------------------
    ! Return the effective stem count contributed by a cohort to a nascent disturbed patch based on
    ! the staged mortalities.
    ! Because of n is per nominal hectare n = stem density (plants/ha).
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use EDTypesMod, only : dump_cohort
    
    ! Arguments:
    type(ed_cohort_type), intent(in), target :: cohort
    
    ! Locals:
    real(r8) :: disturbed_n ! Return value
    real(r8) :: staged_mortality
    real(r8) :: staged_pfrac
    
    ! ----------------------------------------------------------------------------------------------
    
    staged_mortality = 0.0_r8 ! Initialize.
    staged_pfrac = 0.0_r8
    
    ! This is not currently intended to be used along with logging module harvests:
    ! Note: To allow the use of patch_disturbed_basal_area() in development reporting I have
    ! temporarily lowered this to a warning.
    if (cohort%lmort_direct /= 0.0_r8) then
      write(fates_log(),*) 'disturbed_n() is not currently intended to be use in conjunction with traditional logging module events.'
      call dump_cohort(cohort)
      !call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    if (cohort%vm_mort_bole_harvest > 0.0_r8 .and. cohort%vm_mort_in_place > 0.0_r8) then
      write(fates_log(),*) 'disturbed_n(): more than one management mortality type staged.'
      call dump_cohort(cohort)
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    staged_mortality = max(cohort%vm_mort_bole_harvest, cohort%vm_mort_in_place)
    staged_pfrac = max(cohort%vm_pfrac_bole_harvest, cohort%vm_pfrac_in_place)
    if (staged_mortality == 0.0_r8) then
      disturbed_n = cohort%n
    else
      disturbed_n = (cohort%n * staged_pfrac) - (cohort%n * staged_mortality)
    endif
    
  end function cohort_disturbed_n

  !=================================================================================================

  function patch_disturbed_n(patch, pfts)! result(disturbed_n)
    ! ----------------------------------------------------------------------------------------------
    ! Return the stem count for the nascent disturbed patch based on the staged mortalities.
    ! Because n is per nominal hectare n = stem density (plants/ha).
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_patch_type), intent(in), target :: patch ! The patch to be calculated.
    integer(i4), dimension(:), intent(in) :: pfts ! Array of PFT IDs to include in the calculation.
    
    ! Locals:
    real(r8) :: patch_disturbed_n ! Return value
    type(ed_cohort_type), pointer :: current_cohort
    
    ! ----------------------------------------------------------------------------------------------
    
    patch_disturbed_n = 0.0_r8 ! Initialize.
    
    current_cohort => patch%shortest
    do while(associated(current_cohort))
    
      if (any(pfts == current_cohort%pft)) then
        patch_disturbed_n = patch_disturbed_n + cohort_disturbed_n(current_cohort) !disturbed_n(cohort)
      endif
      
      current_cohort => current_cohort%taller
    end do ! Cohort loop.
    
  end function patch_disturbed_n

  !=================================================================================================

  function cohort_disturbed_basal_area(cohort) result(disturbed_basal_area)
    ! ----------------------------------------------------------------------------------------------
    ! Return the effective basal area contributed by a cohort in a nascent disturbed patch based on
    ! the staged potential mortalities.
    !
    ! While minimal this fucntion does reduce calling code complexity and readability.  It may be
    ! inlined anyway.
    !
    !
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_cohort_type), intent(in), target :: cohort
    
    ! Locals:
    real(r8) :: disturbed_basal_area ! Return value
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Only the adjusted stem count is need, the rest of the equation is unchanged: 
    ! pi * ((diameter (cm) / 2) / 100 cm/m)^2 ... 
    disturbed_basal_area = pi_const * (cohort%dbh / 200.0_r8)**2.0_r8 * disturbed_n(cohort)
    
  end function cohort_disturbed_basal_area

  !=================================================================================================

  function patch_disturbed_basal_area(patch, pfts) result(disturbed_basal_area)
    ! ----------------------------------------------------------------------------------------------
    ! Return the effective basal area for a patch after applying any existing staged potential
    ! mortalities, for only the PFTs passed in.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    type(ed_patch_type), intent(in), target :: patch ! The patch to be calculated.
    ! An array of PFT IDs to include in the basal area calculation:
    integer(i4), dimension(:), intent(in) :: pfts
    !integer(i4), dimension(:), intent(in), optional :: pfts ! Should it be optional?
    
    ! Locals:
    real(r8) :: disturbed_basal_area ! Return value
    real(r8) :: total_basal_area ! Temporary!
    type(ed_cohort_type), pointer :: current_cohort
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'patch_disturbed_basal_area():'
    
    disturbed_basal_area = 0.0_r8 ! Initialize.
    total_basal_area = 0.0_r8 ! Temporary!
    
    current_cohort => patch%shortest
      do while(associated(current_cohort)) ! Bad indenting!
        
        if (any(pfts == current_cohort%pft)) then
          disturbed_basal_area = disturbed_basal_area + cohort_disturbed_basal_area(current_cohort) !disturbed_basal_area(current_cohort)
          
          total_basal_area = total_basal_area + (pi_const * (current_cohort%dbh / 200.0_r8)**2.0_r8 * current_cohort%n) ! Temporary!
        !endif
          !if (debug) write(fates_log(), *) 'Cohort is in +++++++++++, PFT = ', current_cohort%pft
        !else
          !if (debug) write(fates_log(), *) 'Cohort is out ----------, PFT = ', current_cohort%pft
        endif
        
        !if (debug) then
          !total_basal_area = total_basal_area + (pi_const * (current_cohort%dbh / 200.0_r8)**2.0_r8 * current_cohort%n) ! Move if kept.
          !write(fates_log(), *) 'Next cohort:'
          !write(fates_log(), *) 'current_cohort%dbh = ', current_cohort%dbh
          !write(fates_log(), *) 'current_cohort%n = ', current_cohort%n, &
          !                      'Effective n = ', cohort_effective_n(current_cohort), &
          !                      'Disturbed n = ', cohort_disturbed_n(current_cohort)
          !write(fates_log(), *) 'Basal Area = ', pi_const * (current_cohort%dbh / 200.0_r8)**2.0_r8 * current_cohort%n, &
          !                      'Disturbed BA = ', disturbed_basal_area, &
          !                      'Total EBA = ', total_basal_area
        !endif
        
        current_cohort => current_cohort%taller
      end do ! Cohort loop.
    
     if (debug) then
       !write(fates_log(), *) 'patch_disturbed_basal_area(): basal area = ', disturbed_basal_area
       !write(fates_log(), *) 'patch_disturbed_basal_area():-----------'
       write(fates_log(), *) 'Patch BA = ', total_basal_area, 'Disturbed BA = ', disturbed_basal_area
     end if
  end function patch_disturbed_basal_area
  
  !=================================================================================================

  function get_flux_profile(cohort) result(flux_profile)
    ! ----------------------------------------------------------------------------------------------
    ! Return the flux profile, if any, associated with mortalities staged in the cohort.
    ! This function will only return a sensible result if one mortality type is present.
    !!
    ! This routine may not identify all cohorts that are subject to logging module mortalities
    ! because some are calculated without mortality rates stored in the cohort.  Work is ongoing to
    ! solve this issue. 
    !
    ! Note: If we changed vm_mort_in_place etc. into arrays using indexes matching the flux profile
    ! IDs these values might be easier to manipulate and extend.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use EDLoggingMortalityMod, only : logging_time ! Make global?
    use EDTypesMod, only : dtype_ilog ! Just for debugging.
    
    ! Arguments:
    type(ed_cohort_type), intent(in), target :: cohort
    
    ! Locals:
    integer :: flux_profile ! Vegetation management flux profile.
    
    ! ----------------------------------------------------------------------------------------------
    
    ! We assume this routine is only called when a managed disturbance has occurred.
    ! This is a safety check since logging_traditional could be returned even for a cohort in an
    ! undisturbed patch given the logic that follows.
    ! This condition occurs in all cases when the this is called during mortality assignment because
    ! cohort%patchptr%disturbance_mode has not been set and equals fates_unset_int.
    ! Is this still helpful in any way then?
    if (debug) then
      if (cohort%patchptr%disturbance_mode /= dtype_ilog) then
        write(fates_log(),*) 'get_flux_profile() called for cohort in patch without managed disturbance.'
        !write(fates_log(),*) 'cohort%patchptr%disturbance_mode = ', cohort%patchptr%disturbance_mode
        ! call dump_patch(cohort%patchptr) ! Overkill?
      end if
    endif
    
    ! Notes from management_fluxes(), integrate:
        ! ------------------------------------------------------------------------------------------
        ! Determine the flux profile corresponding to the managed mortality:
        ! For one mortality we can just identify the only non-zero mortality (or corresponding pfrac
        ! value).
        !
        ! When more that one mortality is applied we would need to identify them all and get each
        ! patch fraction (pfrac) value as well.
        ! In that case 0 = no mortality, 1 = non-disturbing mortality, and the rest are disturbing
        ! mortalities that should = patch_site_areadis.
        ! ------------------------------------------------------------------------------------------
    
    if (cohort%lmort_direct > 0.0_r8 .or. cohort%lmort_collateral > 0.0_r8 .or. &
        cohort%lmort_infra > 0.0_r8 .or. cohort%l_degrad > 0.0_r8) then
      flux_profile = logging_traditional
    else if (cohort%vm_mort_in_place > 0.0_r8) then
      flux_profile = in_place
    else if (cohort%vm_mort_bole_harvest > 0.0_r8) then
      flux_profile = bole_harvest
    else
      ! If logging_time is true and we have not detected any other management mortalities assume
      ! this cohort is subject to traditional logging:
      ! This is a temporary hack and not future safe.  It could eat unmodified cohorts from other
      ! VM events that occur during the same time step as a traditional logging module event.
      if (logging_time) then
        flux_profile = logging_traditional
      else
        flux_profile = null_profile
      endif
    endif
    
    ! If we assume this routine is only called when a managed disturbance has occurred then report
    ! on any cohorts during a logging event that don't have immediately obvious mortality.
    !if (logging_time .and. flux_profile == null_profile) then
    !  write(fates_log(),*) 'Logging event cohort without expected mortalities / degradation.'
    !  call dump_cohort(cohort)
    !end if
    
  end function get_flux_profile

  !=================================================================================================

  subroutine validate_size_specifications(dbh_min_out, dbh_max_out, ht_min_out, ht_max_out, &
                                          dbh_min, dbh_max, ht_min, ht_max)
    ! ----------------------------------------------------------------------------------------------
    ! This routine takes a set of size specifications, performs validity checks on the values and 
    ! combinations passed, and returns appropriate defaults for any missing arguments.
    !
    ! The main purpose of this routine is to reduce code repetition...
    !
    ! The arguments and argument order may seem a bit strange due to the optional nature of the
    ! arguments. It would be most elegant to have the four size specifications as inout but it is
    ! not possible to return a value if it is optional and having to check the values for presence
    ! prior to passing them in would defeat one of the main purposes for this routine.
    !
    ! All size specifications are returned.  Rather that trying to determine if DBH or height is
    ! being used the calling code should just compare if a cohort meets all the criteria. Inclusive
    ! values are given for the unused specifications so only the used specification will actually
    ! matter.
    ! [Add example?????]
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: None
    
    ! Arguments:
    real(r8), intent(out), optional :: dbh_min_out
    real(r8), intent(out), optional :: dbh_max_out
    real(r8), intent(out), optional :: ht_min_out
    real(r8), intent(out), optional :: ht_max_out
    
    ! Size range of to trees to harvest.  Defaults to everything, otherwise some range of sizes,
    ! e.g. > 10cm DBH, < 15 m in height, etc.
    real(r8), intent(in), optional :: dbh_min
    real(r8), intent(in), optional :: dbh_max
    real(r8), intent(in), optional :: ht_min
    real(r8), intent(in), optional :: ht_max
    
    ! Locals:
    real(r8), parameter :: imposibly_small = 0.0_r8 ! Impossibly small value
    
    ! Massive upper size limits based on the largest trees known:
    real(r8), parameter :: dbh_massive = 14050.0_r8 ! Arbol del Tule: Taxodium mucronatum, 14.05 m
    real(r8), parameter :: ht_massive = 115.92_r8 ! Hyperion: Sequoia sempervirens
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Only allow specification of harvest size by DBH or height:
    ! Note: For a single PFT it would be acceptable to specify a size range using a single DBH value
    ! and a single height value.  However, for a mix of PFTs differing allometries means a DBH can
    ! not be mapped uniquely to height and it becomes imposible to perform comparisons safely.
    if ((present(dbh_min) .or. present(dbh_max)) .and. (present(ht_min) .or. present(ht_max))) then
      
      ! Only allow a mix of DBH and height if the values for one essentially include all possible
      ! values.  It is pretty safe to assume that this will only be the case when the values have
      ! already been set by this routine.
      if (.not. ((dbh_min == imposibly_small .and. dbh_max == dbh_massive) .or. &
                 (ht_min == imposibly_small .and. ht_max == ht_massive))) then
        write(fates_log(),*) 'Cannot specify harvest range as a mix of DBH and height.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      else if (debug) then
        write(fates_log(),*) 'It appears these size specifications have been validated previously.'
      end if
    endif
    
    ! Check validity of values that are present and provide sensibel defaults for the missing ones:
    if (present(dbh_min)) then
      if (dbh_min < 0) then
        write(fates_log(),*) 'dbh_min cannot be less than 0. Leave blank for no lower limit.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
      dbh_min_out = dbh_min
    else
      dbh_min_out = imposibly_small ! Impossibly small value.
    endif
    
    if (present(dbh_max)) then
      if (dbh_max > dbh_massive) then
        write(fates_log(),*) 'dbh_max is unrealistically large. Leave blank for no upper limit.'
        ! We could just warn here as this is unlikely to have any negative effects.
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
      dbh_max_out = dbh_max
    else
      dbh_max_out = dbh_massive ! Impossibly large value.
    endif
    
    if (present(ht_min)) then
      if (ht_min < 0) then
        write(fates_log(),*) 'ht_min cannot be less than 0. Leave blank for no lower limit.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
      ht_min_out = ht_min
    else
      ht_min_out = imposibly_small ! Impossibly small value.
    endif
    
    if (present(ht_max)) then
      if (ht_max > ht_massive) then
        write(fates_log(),*) 'ht_max is unrealistically large. Leave blank for no upper limit.'
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
      ht_max_out = ht_max
    else
      ht_max_out = ht_massive ! Impossibly large value.
    endif
    
    ! Check that values define valid ranges:
    if (dbh_min_out > dbh_max_out) then
      write(fates_log(),*) 'The maximum DBH specification is smaller that the minimum.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    if (ht_min_out > ht_max_out) then
      write(fates_log(),*) 'The maximum height specification is smaller that the minimum.'
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
  end subroutine validate_size_specifications

  !=================================================================================================
  ! Patch and Cohort Search Subroutines:
  !=================================================================================================

!   function patch_is_bare(patch)
!     ! ----------------------------------------------------------------------------------------------
!     ! 
!     ! ----------------------------------------------------------------------------------------------
!     
!     ! Uses:
!     
!     ! Arguments:
!     
!     ! Locals:
!     
!     ! ----------------------------------------------------------------------------------------------
!     
!   end function path_is_bare

! patch_is_forest?
! patch_is_managed_pine?
! patch_is_grassland?

  !=================================================================================================
  ! VM Prescribed Event Driver File Subroutines:
  !   Vegetation Management Events may be prescribed using an optional event driver file.
  !
  ! Driver File Format Draft 2:---------------------------------------------------------------------
  !   The driver file uses a custom text format designed to be human readable and easily written by
  ! hand.  This driver file is intended for hands on experimentation at the single point to small
  ! regional simulation scale.  Therefore the format attempts to be somewhat flexible and tolerant
  ! of small format irregularities that users might introduce.
  !
  !   Global gridded input will be accomplished via extension of the existing land use driver file
  ! behaviors.  (In development, see heuristics.)
  !
  ! General File Structure:
  !   There should be at least two lines.  The first (non-comment) line should be a header line that
  ! specifies the field names.  The following lines are event lines.
  !
  ! Example:
  ! Date        Lat  Lon  EventSpec
  ! YYYY-MM-DD  X.X  Y.Y  do_something(x = 3.0, y = 2, z = 13.1)
  !
  ! Blank lines:
  !   Blank lines are ignored and can occur anywhere in the file.
  !
  ! Comments:
  !   Any content following a '!' on a line is ignored.
  !
  ! Header Line:
  !   The header line is a soft requirement.  It is only used for readability and we do minimal
  ! checking of it.
  !   Ideally the the header will occur as the first line of content.  Comment lines describing the
  ! file (e.g. creator, date, project, etc) may precede it (a good idea) but all events should
  ! follow it.  This is a style guidline and is not currently enforced.
  !   The header may be missing (bad style) or may occur more than once (might be useful for really
  ! long files).
  !   The column names beyond 'Date' are not checked and their values will not change anything about
  ! how the file will be read in, but it is good style to include them.
  !   The header line must start with 'Date".  It can be be preceded by white space and some case
  ! variants will be accepted but both are bad style and may be deprecated in the future.
  !
  ! Field Delimiters:
  !   The header columns and event fields are delimited by one or more spaces.  For flexibility
  ! column widths are not specified.  While aligning columns for readability is encouraged it is not
  ! required. Any number of spaces between fields will be accepted.
  !   Tabs or other whitespace may cause weirdness.  Do not use them.
  !
  ! Field Order:
  !   The field order is the date of the event, the latitude and longitude (in model grid units)
  ! where the event should occur, followed by an event specification, which is written like a
  ! Fortran function call (see example above).
  !
  ! Event Fields:
  !   The date and lat/lon formats are flexible.  See is_now() and is_here() for more information.
  !   
  !   The event specification is formated like a Fortran function call and maps the event to a site
  ! level subroutine that in most cases will have the same name.  See vm_event%load() for more
  ! information.
  !
  ! Namelist Variables:-----------------------------------------------------------------------------
  !   The use of the VM Event Driver File is turned by setting the use_fates_vm_driver_file = .true.
  ! in the user namelist.  When turned on the path to a valid VM Event Driver File should be
  ! specified using fates_vm_driver_filepath in the user namelist file.
  !
  ! Event Globals: (vm_generative_event & vm_mortality_event)---------------------------------------
  !   Events read from the driver file are stored in (module) globals that are then executed along
  ! with events generated by other means.
  !
  ! Notes:------------------------------------------------------------------------------------------
  ! - Important: This format is a work in progress and may change further.  The file format is
  ! denoted a draft 2 not because it is not working, but because it is prerelease and may change.
  !
  ! - The line length should probably be kept to under 100 characters to be compatible will
  ! different Fortran implementations.  The code does not enforce that though.  In fact it uses
  ! string lengths internally that are questionably long.
  !
  ! ToDo: Scriptable Format Variant.
  !   The format is designed to be easy to read, but is not perfect for writing with a script.  The
  ! leading fields are very compatible with a script or spreadsheet but writing the event
  ! specification would require a bit of work if it contains arguments that must be scripted.
  ! [Actually I have done this in R and it isn't too bad.]
  ! The solution I foresee is a format variant that is the same for the leading fields but specifies
  ! the event and its arguments using named columns.  Some lines would have blank entries depending
  ! on the event type.   Adding handling of comma separated and tab delimited formats would be
  ! pretty easy and would finish up the variant.
  !
  ! Example: (spaces added for readability)
  ! Date,       Lat, Lon, EventCode,    x,   y, z
  ! YYYY-MM-DD, X.X, Y.Y, do_something, 3.0, 2, 13.1
  !
  !   Implementing this will wait on the need however.
  !=================================================================================================

  subroutine load_prescribed_events(site)
    ! ----------------------------------------------------------------------------------------------
    ! Check for the presence of an input file that prescribes a series of VM events.  If present
    ! load any events for this timestep and location into memory.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use FatesInterfaceTypesMod, only : hlm_use_vm_driver_file, hlm_vm_driver_filepath
    use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit
    
    ! Arguments:
    type(ed_site_type), intent(in), target :: site ! The current site object.
    
    ! Locals:
    logical :: driver_file_exists ! Does the VM driver file exist?
    logical :: driver_file_open ! Is the VM driver open?
    integer :: driver_file_unit ! File unit for the VM driver file
    integer :: io_status ! IO error flag
    integer :: comment_index ! Where does a comment start on a line?
    
    character(len = line_strlen) :: line_str ! Holds lines of data from the file
    character(len = line_strlen) :: date_str, lat_str, lon_str ! Field values
    type(vm_event) :: the_event
    
    ! ----------------------------------------------------------------------------------------------
    if (debug) write(fates_log(), *) 'load_prescribed_events() entering.'
    
    ! Initialize event globals in any case:
    call vm_generative_event%zero()
    call vm_mortality_event%zero()
    
    ! Note: Indent the following block!
    if (hlm_use_vm_driver_file) then
      
      ! Initialize an empty event to hold any events found:
      ! Note: There could be more than one event.  We probably need to reinitialize again below.
      call the_event%zero()
      
      ! Check if the driver file exists and is not in use:
      inquire(file = trim(hlm_vm_driver_filepath), exist = driver_file_exists, opened = driver_file_open)
      
      ! If the file exists proceed otherwise don't. We don't use a flag:
      if (driver_file_exists) then
        
        if (driver_file_open) then
          write(fates_log(),*) 'The vegetation management driver file is open already.'
          call endrun(msg = errMsg(__FILE__, __LINE__))
        else
          
          ! Open the file:
          driver_file_unit = shr_file_getUnit()
          open(unit = driver_file_unit, file = trim(hlm_vm_driver_filepath), status = 'OLD', &
               action = 'READ', form = 'FORMATTED')
          rewind(driver_file_unit)
          
          ! Read the first event line:
          read(driver_file_unit, fmt='(A)',iostat = io_status) line_str
          
          ! An empty file will not cause problems but is likely not expected so warn:
          if (io_status /= 0) then
            write(fates_log(),*) 'The vegetation management prescribed event driver file was found but contains no events.'
          else
            
            ! Read through line by line looking for events specified for this timestep and location:
            do while (io_status == 0)
              if (debug) write(fates_log(),*) 'Parsing event line:'
              
              ! Ignore anything following a '!' as a comment:
              comment_index = index(line_str, '!')
              if (comment_index /= 0) then
                line_str = line_str(:comment_index-1)
                if (debug) write(fates_log(), *) 'Comment ignored.'
              endif
             
              !There shouldn't be any leading whitespace but allow it. Trailing whitespace could
              ! cause parsing issues:
              line_str = trim(line_str)
              
              ! Ignore blank lines (whitespace only) and lines that contain only comments:
              ! Will this handle tabs?
              if (len_trim(line_str) == 0) then
                if (debug) write(fates_log(), *) 'Skipping blank or commented line.' ! Temporary!
              else
                
                ! Extract the leading fields:
                date_str = field_pop(line_str)
                lat_str = field_pop(line_str)
                lon_str = field_pop(line_str)
                
                ! Check if this is the header line:
                if (date_str == "Date" .or. date_str == "date" .or. date_str == "DATE") then
                  if (debug) write(fates_log(),*) 'The VM event driver header:', trim(line_str)
                  ! Could add additional error checking here.
                else
                  
                  ! Check if this event matches this timestep and location:
                  if (is_now(date_str) .and. is_here(lat_str, lon_str, site)) then
                    if (debug) write(fates_log(), *) 'VM event matches date and location.'
                    
                    ! Load the event:
                    call the_event%load(line_str)
                    
                    ! Is the event a mortality inducing or generative (fecundity) event?
                    ! Only allow one mortality but more than one planting?
                    if (the_event%is_generative()) then
                      if (vm_generative_event%code /= vm_event_null) then
                        write(fates_log(),*) 'A VM generative event already exits for this time step.'
                        ! Note: In future multiple generative events will be allowed to co-occur.
                        call endrun(msg = errMsg(__FILE__, __LINE__))
                      else
                        vm_generative_event = the_event
                        if (debug) then
                          write(fates_log(),*) 'VM generative event loaded.'
                          call vm_generative_event%dump()
                        endif
                      endif
                    else ! Mortality event:
                      if (vm_mortality_event%code /= vm_event_null) then
                        write(fates_log(),*) 'Only one VM mortality event per time step is currently allowed per site.'
                        call endrun(msg = errMsg(__FILE__, __LINE__))
                      else
                        vm_mortality_event = the_event
                        if (debug) then
                          write(fates_log(),*) 'VM mortailty event loaded.'
                          call vm_mortality_event%dump()
                        endif
                      endif ! vm_mortality_event%code /= vm_event_null)
                    endif ! (the_event%is_generative())
                  endif ! (is_now(date_str) ...
                endif ! (date_str == "Date"...
              endif ! (len_trim(line_str) == 0)
              
              ! Read the next line:
              read(driver_file_unit, fmt='(A)',iostat = io_status) line_str
            end do ! while (io_status == 0)
            
          endif ! (io_status /= 0)
          
          ! Close and clean up:
          close(driver_file_unit, iostat = io_status)
          if (io_status /= 0) then
            write(fates_log(),*) 'The vegetation management driver file has failed to close.'
            call endrun(msg = errMsg(__FILE__, __LINE__))
          endif
          call shr_file_freeUnit(driver_file_unit)
          
        endif ! (driver_file_open)
        
      else
        write(fates_log(),*) 'A vegetation management prescribed event driver file was not found.'
      endif ! (driver_file_exists)
    endif ! (hlm_use_vm_driver_file)
    
    if (debug) write(fates_log(), *) 'load_prescribed_events() exiting.'
  end subroutine load_prescribed_events

  !=================================================================================================

  function field_pop(line_str) result(field_str) ! field_pop_str()?
    ! ----------------------------------------------------------------------------------------------
    ! Pop the first field off the front (left) of a line of fields encoded as text.
    ! Returns the field as the return value and removes it from the line string passed in.
    !
    ! The delimiter between fields is one or more spaces.
    ! Note: We could add an optional delimiter argument to expand the utility of this routine.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    character(len = *), intent(inout) :: line_str
    
    ! Locals:
    character(len = line_strlen) :: field_str ! Return value
    integer :: delim_index
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Trim any leading (and trailing) whitespace:
    line_str = adjustl(line_str) ! This pushes whitespace to the end of the line, which we remove.
    line_str = trim(line_str)
    
    ! Find the the following delimiting block of 1 or more spaces:
    delim_index = index(line_str, " ")
    
    ! If 0 this could be the last field, make sure something is there:
    if (delim_index == 0) then
      if (len(line_str) > 0) then
        field_str = line_str
        line_str = ''
      else
        field_str = ''
        line_str = ''
        ! This may not be expected, warn:
        write(fates_log(),*) 'field_pop(): No fields remain.'
      endif
      
    else
      field_str = line_str(:delim_index - 1)
      line_str = adjustl(line_str(delim_index:)) ! Remove the popped field from the line.
    endif
    
  end function field_pop

  !=================================================================================================

  function field_pop_int(line_str) result(field_int)
    ! ----------------------------------------------------------------------------------------------
    ! Pop the first field off the front (left) of a line and return it as an integer.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    character(len = *), intent(inout) :: line_str
    
    ! Locals:
    character(len = line_strlen) :: field_str ! Intermediate
    integer :: field_int ! Return value
    
    ! ----------------------------------------------------------------------------------------------
    
    field_str = field_pop(line_str)
    read(field_str, *) field_int
    
  end function field_pop_int

  !=================================================================================================

  function field_pop_real(line_str) result(field_real)
    ! ----------------------------------------------------------------------------------------------
    ! Pop the first field off the front (left) of a line and return it as a real.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    character(len = *), intent(inout) :: line_str
    
    ! Locals:
    character(len = line_strlen) :: field_str ! Intermediate
    real(r8) :: field_real ! Return value
    
    ! ----------------------------------------------------------------------------------------------
    
    field_str = field_pop(line_str)
    read(field_str, *) field_real
    
  end function field_pop_real

  !=================================================================================================

  function is_now(date_string)
    ! ----------------------------------------------------------------------------------------------
    ! Parse the string passed in as a date and determine if it matches the current time step.
    !
    ! The date string format is somewhat flexible. We expect that the event file may be written
    ! manually in some cases or exported from a spread sheet, etc.  We don't want execution to fail
    ! because a line is not formated in particular way or doesn't match the across lines.
    !
    ! Rules:----------------------------------------------------------------------------------------
    ! Field Order:
    ! - Year first, then month, then day.
    ! Field Separators:
    ! - Currently we handle slash, dash, and period (/,-,.) but not space, since it is our parental
    ! (data line) field delimiter.
    ! Digits:
    ! - We allow 1 or 2 digit month and day strings.  Leading 0s Are allow but not required.
    ! - We allow years of up to 10 digits (in case you have a huge allocation).
    ! Thousands separator:
    ! - No commas or periods to separate 1000s in years are allowed.
    ! Valid examples:
    ! - YYYY-MM-DD, YYYY/MM/DD, 00YY.0M.0D, YY-M-D, YYYYYYYY/MM/DD
    ! - Even this would work: 0000YYYY/M.DD
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use FatesInterfaceTypesMod, only : hlm_current_year, hlm_current_month, hlm_current_day
    
    ! Arguments:
    character(len = *), intent(in) :: date_string
    
    ! Locals:
    logical :: is_now ! Return value
    character(len = len(date_string)) :: work_string ! Modifiable local copy of date_string
    character(len = *), parameter :: delim_chars = "/-."
    integer :: delim_index
    character(len = 10) :: year_str
    character(len = 2) :: month_str, day_str ! Should we increase this to allow errors like 0DD?
    integer :: year, month, day ! Event date components
    
    ! ----------------------------------------------------------------------------------------------
    
    work_string = date_string ! Subsequent operations modify the string in place so make a copy.
    work_string = trim(work_string)! This shouldn't be necessary currently but doesn't hurt.
    
    ! Add error checking...
    ! Check for length and invalid characters!!!!!
    
    ! Parse the date string:
    delim_index = scan(work_string, delim_chars)
    year_str = work_string(:delim_index-1)
    read(year_str, *) year
    
    ! Month:
    work_string = work_string(delim_index+1:)
    delim_index = scan(work_string, delim_chars)
    month_str = work_string(:delim_index-1)
    read(month_str, *) month
    
    ! Day:
    day_str = work_string(delim_index+1:)
    read(day_str, *) day
    
    if (debug) then
      write(fates_log(), *) 'Date string: ', trim(date_string)
      write(fates_log(), *) 'Year:  ', year
      write(fates_log(), *) 'Month: ', month
      write(fates_log(), *) 'Day:   ', day
    endif
    
    if (hlm_current_year == year .and. hlm_current_month == month .and. &
        hlm_current_day == day) then
      is_now = .true.
    else
      is_now = .false.
    endif
  end function is_now

  !=================================================================================================

  function is_here(lat_string, lon_string, site) ! event_is_here()?
    ! ----------------------------------------------------------------------------------------------
    ! Determine if the current site matches the specified coordinates.
    !
    ! Currently the coordinates are matched with a small tolerance.  It would be better if we knew
    ! the actual bounds of the grid cell but those are not available currently.  The tolerance is
    ! arbitrary.  We assume that in most cases FATES will be run at a coarse resolution.
    ! 
    ! Rules:----------------------------------------------------------------------------------------
    ! Format:
    ! - Decimal degrees only.  No sexagesimal degrees (minutes and seconds).
    ! - Bare numbers, no degree symbols.
    ! - No cardinal directions (i.e. N,S,E,W), use sign instead.
    ! - No spaces.
    ! Valid value ranges:
    ! - Latitude: -90 to 90
    ! - Longitude may be specified as -180 to 180 (common) or 0 to 360 (netCDF / model world).
    ! Special values:
    ! - Passing -999 for both latitude and longitude specifies everywhere.
    ! 
    ! ToDo ranges:
    ! It would be pretty simple to allow lat and long ranges to specify a simple rectangular region.
    ! Because of negative degrees ranges would have to be specified in some other way, perhaps
    ! with commas "34.5, -70.1" or with Fortran array specification style "[34.5:-70.1]".
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use EDTypesMod, only : ed_site_type
    
    ! Arguments:
    character(len = *), intent(in) :: lat_string
    character(len = *), intent(in) :: lon_string
    type(ed_site_type), intent(in), target :: site ! The current site object.
    
    ! Locals:
    logical :: is_here ! Return value
    real(r8) :: latitude, longitude_in, longitude_360
    real(r8), parameter :: degree_tolerance = 0.1_r8 ! Arbitrary degree tolerance for matching.
    
    ! ----------------------------------------------------------------------------------------------
    
    ! If ranges were allowed we would check for those first.
    ! e.g. if (index(lat_string, ',') /= 0) then
    
    ! Convert to numbers:
    read(lat_string, *) latitude
    read(lon_string, *) longitude_in
    
    if (debug) then
      write(fates_log(), *) 'Latitude string: ', trim(lat_string)
      write(fates_log(), *) 'Latitude string: ', trim(lon_string)
      write(fates_log(), *) 'Latitude  #:     ', latitude
      write(fates_log(), *) 'Longitude #:     ', longitude_in
    endif
    
    ! Check to see if the location of the event is "everywhere":
    if (latitude == -999.0_r8 .and. longitude_in == -999.0_r8) then
      is_here = .true.
      if (debug) write(fates_log(), *) 'VM event specified as occurring everywhere.'
    else
      ! Validate the coordinates and convert to model coordinates if necessary.
      
      ! Note: Checks are similar to code in FatesInventoryInitMod: assess_inventory_sites():
      if (latitude < -90.0_r8 .or. latitude > 90.0_r8) then
        write(fates_log(),*) 'Vegetation management event has invalid latitude: ', latitude
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
      
      ! Convert the longitude to 0-360 degree notation if necessary.
      if (longitude_in < 0.0_r8) then
        longitude_360 = longitude_in + 360.0_r8
      else
        longitude_360 = longitude_in
      endif
      
      ! Check that the longitude is valid:
      if (longitude_360 < 0.0_r8 .or. longitude_360 > 360.0_r8) then
        write(fates_log(),*) 'Vegetation management event has invalid longitude (converted): ', longitude_in
        call endrun(msg = errMsg(__FILE__, __LINE__))
      endif
      
      ! Compare the latitude and longitude to that of the site being simulated:
      if (abs(site%lat - latitude) < degree_tolerance .and. &
          abs(site%lon - longitude_360) < degree_tolerance) then
        is_here = .true.
        if (debug) write(fates_log(), *) 'VM event matches this location.'
      else
        is_here = .false.
        if (debug) write(fates_log(), *) "VM event doesn't match this location."
      endif
      
    endif
    
  end function is_here

  !=================================================================================================
  ! vm_event Type Member Routines:
  !=================================================================================================

  ! vm_event%zero():
  subroutine zero(this) ! blank()?  Init?
    ! ----------------------------------------------------------------------------------------------
    ! vm_event Type Bound Procedure:
    !   Set data members to default 'blank' (missing or default) values.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    class(vm_event), intent(inout) :: this ! Self reference.  Routine modifies self.
    
    ! Locals: NA
    
    ! ----------------------------------------------------------------------------------------------
    
    this%code = vm_event_null
    
    ! We need sensible defaults for these parameters or a value that indicates the parameter is
    ! 'missing'.  These value are temporary:
    this%pfts = vm_empty_integer ! Change to array (16 in length?)!!!!!
    this%density = vm_empty_real
    this%dbh = vm_empty_real
    this%height = vm_empty_real
    this%row_fraction = vm_empty_real
    this%final_basal_area = vm_empty_real
    this%thin_fraction = vm_empty_real
    this%dbh_min = vm_empty_real
    this%ht_min = vm_empty_real
    this%patch_fraction = vm_empty_real
    
  end subroutine zero

  !=================================================================================================

  ! vm_event%load():
  subroutine load(this, event_str) ! load_from_string(), initialize()
    ! ----------------------------------------------------------------------------------------------
    ! vm_event Type Bound Procedure:
    !   Initialize the event object from an event specification string (from a VM driver file).
    !
    ! Event Specification String Format:------------------------------------------------------------
    !   The event driver specification is formated like a Fortran function call.
    !
    ! Examples:
    ! event_name(arg1 = val1, arg2 = [val2.1, val2.2, val2.3], arg3 = val3)
    ! do_something(dbh_min = 3.0, dbh_max = 55.5, z = 13.1)
    !
    ! Event Type or Name String:
    !   The event specification starts with an event name (i.e. the "function name").
    !   This event name must match an existing event type.  Event specifications map to site level
    ! subroutines that execute the event.  These subroutines declare their interface: the arguments
    ! and values that they take.
    !
    ! Arguments:
    !   Arguments follow the event name, enclosed in parentheses, and separated by commas.
    !   Arguments may be optional, depending on the event, but all arguments that are present must
    ! be passed by name.  Name value pairs are separated by '='.
    !   The order of arguments does not matter but it is best to keep it same as the order in the
    ! event's interface declaration.
    !   The names and meanings of arguments are standardized across the set of available site level
    ! events.  Shared argument names should mean the same thing across events.
    !   Some arguments can take arrays of values.  These are specified using Fortan 2003 style
    ! square brackets, i.e. name = [value1, value2].  Single values will be accepted for array
    ! arguments.
    !
    ! White Space:
    !   Spaces are allowed between elements in event specifications for readability but are not
    ! required.  Other whitespace should be avoided.
    !   The style suggestion is to include spaces between arguments and name values pairs but not
    ! next to parentheses.
    !   Harder to read:  do_something ( dbh_min=3.0,dbh_max=55.5,z=13.1 )
    !   Earlier to read: do_something(dbh_min = 3.0, dbh_max = 55.5, z = 13.1)
    !
    ! Case:
    !   Currently event and argument names are case sensitive and follow FATES conventions of using
    ! lower_case_separated_by_underscores.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    class(vm_event), intent(inout) :: this ! Self reference.  Routine modifies self.
    character(len = *), intent(in) :: event_str
    
    ! Locals:
    character(len = len(event_str)) :: event_type_str, arguments_string, param_string, param_name, param_value, end_string
    integer :: i ! Iterator
    integer :: delim_index

    ! ----------------------------------------------------------------------------------------------
    
    ! Some of the following could all be simplified with field_pop() if it took a delimiter.
    
    ! Parse the event type / name string (i.e the "function name" of the specification string):
    delim_index = index(event_str, '(')
    if (delim_index == 0) then
      write(fates_log(),*) 'Improperly formed VM event specification. Missing opening parenthesis: ', event_str
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    event_type_str = event_str(:delim_index-1)
    ! Remove any surrounding whitespace:
    event_type_str = adjustl(event_type_str)
    event_type_str = trim(event_type_str)
    
    !if (debug) write(fates_log(),*) 'event_type_str: ', trim(event_type_str) ! Temporary!!!!!
    
    select case (event_type_str)
      case ('plant')
        this%code = vm_event_plant
      case ('thin_row_low') ! Name will probably change!!!!!
        this%code = vm_event_thin_test1
      case ('thin_proportional')
        this%code = vm_event_thin_proportional
      case ('thin_low_perfect')
        this%code = vm_event_thin_low_perfect
      case ('thin_low_probabilistic')
        this%code = vm_event_thin_low_probabilistic
      case ('clearcut')
        this%code = vm_event_clearcut
      !case ()
      case default
        write(fates_log(),*) 'VM event name is not recognised:', event_type_str
        call endrun(msg = errMsg(__FILE__, __LINE__))
    end select
    
    ! Get the string segment containing all the arguments (from '(' to ')'):
    arguments_string = event_str(delim_index+1:)
    delim_index = index(arguments_string, ')')
    if (delim_index == 0) then
      write(fates_log(),*) 'Improperly formed VM event specification. Missing closing parenthesis: ', arguments_string
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    arguments_string = arguments_string(:delim_index-1)
    ! Remove any whitespace inside the parentheses:
    arguments_string = adjustl(arguments_string)
    arguments_string = trim(arguments_string)
    
    !if (debug) write(fates_log(),*) 'arguments_string: ', trim(arguments_string) ! Temporary!!!!!
    
    ! There should be no content after the closing parenthesis (any comments have already been removed):
    end_string = arguments_string(delim_index+1:)
    if (len_trim(end_string) > 0) then
      write(fates_log(),*) 'Improperly formed VM event specification. Text after closing parenthesis: ', end_string
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
    ! Parse arguments / parameters:
    ! The number of parameters is variable and may be zero.  Remove each one in turn:
    do while (len_trim(arguments_string) /= 0)
      
      ! Get the next name value pair:
      delim_index = index(arguments_string, ',')
      ! If there is no comma then we are on the last argument:
      if (delim_index /= 0) then
        param_string = arguments_string(:delim_index-1)
        arguments_string = arguments_string(delim_index+1:) ! Remove the processed argument.
      else
        param_string = arguments_string
        arguments_string = ''
      endif
      
      !if (debug) write(fates_log(),*) 'param_string: ', trim(param_string) ! Temporary!!!!!
      
      ! Parse:
      ! The argument name value pairs are separated by equals signs.  A variable amount of
      ! whitespace may also present between these elements.  Arrays may be passed as values using square brackets.
      ! name = value or name = [value1, value2]
      delim_index = index(param_string, '=')
      param_name = param_string(:delim_index-1)
      param_name = adjustl(param_name)
      param_name = trim(param_name)
      
      !if (debug) write(fates_log(),*) 'param_name: ', trim(param_name) ! Temporary!!!!!
      
      param_value = param_string(delim_index+1:)
      ! The following may not be needed since Fortran interprets numeric values pretty robustly:
      param_value = adjustl(param_value)
      param_value = trim(param_value)
      
      !if (debug) write(fates_log(),*) 'param_value: ', trim(param_value) ! Temporary!!!!!
      
      ! The meaning of argument names must be consistent (at least in terms of type) across all routines that use them.
      select case (param_name)
        case ('pfts') ! Many...
          ! This may be a single value or array but for now we assume there is only one (array of 1).
          ! We also should allow group names.
          read(param_value, *) this%pfts
          !pfts = parse_array(param_value)
        
        case ('density') ! plant()
          read(param_value, *) this%density  
        case ('dbh') ! plant()
          read(param_value, *) this%dbh
        case ('height') ! plant()
          read(param_value, *) this%height
        case ('row_fraction') ! thin_row_low()
          read(param_value, *) this%row_fraction
        case ('final_basal_area') ! thin_row_low()
          read(param_value, *) this%final_basal_area
        case ('thin_fraction') ! thin_proportional(), thin_patch_low_probabilistic()
          read(param_value, *) this%thin_fraction
        case ('dbh_min') ! clearcut()
          read(param_value, *) this%dbh_min
        case ('ht_min') ! clearcut()
          read(param_value, *) this%ht_min
        case ('patch_fraction') ! clearcut()
          read(param_value, *) this%patch_fraction
        !More to come:
        !case ('')
        !case ('where')
        case default
          write(fates_log(),*) 'VM parameter name is not recognised.'
          call endrun(msg = errMsg(__FILE__, __LINE__))
      end select
    enddo ! (len_trim(arguments_string) /= 0)
    
    if (debug) then
      write(fates_log(), *) 'Loaded VM event from driver file:'
      call this%dump()
    endif
    
  end subroutine load

  !=================================================================================================

  ! vm_event%is_generative():
  function is_generative(this) ! This doesn't roll off the tongue. is_mortality()?????
    ! ----------------------------------------------------------------------------------------------
    ! vm_event Type Bound Procedure:
    !   Is this event a generative (e.g. planting) or mortality (e.g. harvest) event?
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    class(vm_event), intent(in) :: this ! Self reference
    logical :: is_generative ! Return value
    
    ! Locals: NA
    
    ! ----------------------------------------------------------------------------------------------
    
    if (this%code > vm_event_null .and. this%code <= vm_event_generative_max) then
      is_generative = .true.
    else if (this%code <= vm_event_mortality_max) then
      is_generative = .false.
    else
      write(fates_log(),*) 'Unrecognized or null event code:', this%code
      call endrun(msg = errMsg(__FILE__, __LINE__))
    endif
    
  end function is_generative

  !=================================================================================================

  ! vm_event%dump():
  subroutine dump(this)
    ! ----------------------------------------------------------------------------------------------
    ! vm_event Type Bound Procedure:
    !   Dump the object contents to the log.
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses: NA
    
    ! Arguments:
    class(vm_event), intent(in) :: this ! Self reference.
    
    ! Locals: NA
    
    ! ----------------------------------------------------------------------------------------------
    
    write(fates_log(), *) 'Dumping Vegetation Management Event:--------------'
    write(fates_log(), *) 'Event code:       ', this%code
    
    write(fates_log(), *) 'Parameters:'
    write(fates_log(), *) 'pfts:             ', this%pfts
    write(fates_log(), *) 'density:          ', this%density
    write(fates_log(), *) 'dbh:              ', this%dbh
    write(fates_log(), *) 'height:           ', this%height
    write(fates_log(), *) 'row_fraction:     ', this%row_fraction
    write(fates_log(), *) 'final_basal_area: ', this%final_basal_area
    write(fates_log(), *) 'thin_fraction:    ', this%thin_fraction
    write(fates_log(), *) 'dbh_min:          ', this%dbh_min
    write(fates_log(), *) 'ht_min:           ', this%ht_min
    write(fates_log(), *) 'patch_fraction:   ', this%patch_fraction
    
  end subroutine dump

end module FatesVegetationManagementMod
