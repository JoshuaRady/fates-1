!===================================================================================================
! FatesVegetationManagementMod.F90
! Implemented by: Joshua M. Rady
! Started: 8/26/2020
!
! This module contains subroutines to simulate aspects of human management of vegetation.  While
! written to implement forest management the code contains generic behaviors that could be used to
! implement other activities, e.g. agriculture.
!
! Note: These notes should be revised to conform to FATES coding style, to the extent that exists.
! The coding style is a bit vague with a mix of 2 and 3 space identing.  I will use 2 since they
! indent starting at the module scope.
!
!===================================================================================================

module FatesVegetationManagementMod
  
  use EDTypesMod, only : ed_site_type, ed_patch_type, ed_cohort_type
  use FatesAllometryMod, only : h2d_allom, h_allom
  use FatesConstantsMod, only : r8 => fates_r8
  ! Using generic integers for PFT numbers should be fine but we follow FatesAllometryMod in
  ! explicitly sizing them.
  use FatesConstantsMod, only : i4 => fates_int
  use FatesConstantsMod, only : itrue, ifalse
  use FatesInterfaceTypesMod, only : bc_in_type
  use PRTGenericMod, only : prt_vartypes
  
  ! Logging and error reporting:
  use shr_log_mod, only : errMsg => shr_log_errMsg
  use FatesGlobals, only : fates_log
  
  ! Enforce explicit type declarations:
  implicit none
  private
  
  ! public :: managed_mortality
  public :: managed_fecundity
  
  public :: plant
  public :: init_temporary_cohort ! Currently doen't need to be public.
  
  public :: anthro_disturbance_rate
  public :: anthro_mortality_rate
  public :: management_fluxes
  
  ! character(len=*), parameter, private :: sourcefile = __FILE__
  
  ! Debugging flag for module:
  logical, parameter :: debug = .true.
  
  !=================================================================================================
  
contains
  
  ! Placeholder:
  !subroutine managed_mortality
    ! ----------------------------------------------------------------------------------------------
  !end subroutine
  
  subroutine managed_fecundity(site, bc_in) ! managed_planting?
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
    use FatesInterfaceTypesMod , only : hlm_current_year ! For testing.
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: site
    type(bc_in_type), intent(in) :: bc_in
    
    ! Locals:
    type(ed_patch_type), pointer :: thisPatch
    
    ! ----------------------------------------------------------------------------------------------
    
    if (debug) then
      write(fates_log(), *) 'managed_fecundity() entering.'
    end if
    
    ! ToDo: 
    ! Determine which patch or patches to plant into:
    ! Considerations include the area that needs to be planted, the composition of the patches: i.e.
    ! what species are there, is it a bare patch, and history: is it secondary or managed land.
    !
    ! For initial testing the patch doesn't matter much but we want to make sure it is big enough
    ! that any new cohorts are not eliminated.
    
    ! For testing we use the logging event code to trigger a planting event:
    if (logging_time) then
      
      if (debug) then
        write(fates_log(), *) 'Planting triggered by logging event (for testing).'
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
        call plant(site = site, patch = thisPatch, bc_in = bc_in, pft_index = 2)
        
    end if ! if (logging_time)
    
    if (debug) then
      write(fates_log(), *) 'managed_fecundity() exiting.'
    end if
    
  end subroutine
  
  !=================================================================================================
  ! Planting Subroutines:
  !=================================================================================================
  
  subroutine plant(site, patch, bc_in, pft_index, density, dbh, height) !plant_sapling()?
    ! ----------------------------------------------------------------------------------------------
    ! Plant a new cohort with the PFT, density, and size specified into an existing patch.
    ! Planting differs from recruitment in that the new plants come from somewhere else, currently
    ! not specified, rather than the seed bank.  There are currently no constrains of the size of
    ! the plantings, however the expected use will be primarily to add small plants, i.e. seedlings.
    ! ----------------------------------------------------------------------------------------------
    
    !use EDTypesMod.F90, only : ed_patch_type
    use EDCohortDynamicsMod, only : create_cohort, zero_cohort, InitPRTObject
    use EDPftvarcon, only : EDPftvarcon_inst
    use EDTypesMod, only : num_elements, site_massbal_type, element_list
    use PRTGenericMod, only : num_organ_types
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: site
    ! type(ed_patch_type), intent(inout), target :: patch
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
    
    ! Note: It appears that the compliler on Cheyenne does not support the value keyword.
    ! real(r8), intent(in), optional, value :: density ! The planting density in plants / m^2 or hectare??????
    real(r8), intent(in), optional :: density ! The planting density (plants / m^2)
    real(r8), intent(in), optional :: dbh ! Sapling diameter at breast height (cm)
    real(r8), intent(in), optional :: height ! Sapling height (m)
    
    ! Locals:
    real(r8) :: the_density ! The actual planting density (plants / m^2)
    real(r8) :: the_dbh ! The actual DBH (cm)
    real(r8) :: the_height ! The actual height (m)
    real(r8) :: area ! Patch / cohort area (m^2)
    logical :: use_dbh ! Specify cohort size by DBH, otherwise use height
    type(ed_cohort_type), pointer :: planting_cohort   ! Temporary cohort for calculations
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
    if (present(density)) then
      the_density = density
    else
      the_density = EDPftvarcon_inst%initd(pft_index)
    end if
    
    ! Size:
    ! If DBH or height if provided (in that order), otherwise use the default height.
    if (present(dbh)) then
      the_dbh = dbh
      use_dbh = .true.
    else
      if (present(height)) then
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
    
    ! Add the cohort to the patch...
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
    
    !use EDCohortDynamicsMod, only : zero_cohort, InitPRTObject  [Moved to calling subroutine.]
    use EDPftvarcon, only : EDPftvarcon_inst
    use EDTypesMod, only : num_elements, element_list ! site_massbal_type
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
    ! use PRTGenericMod, only : prt_vartypes !  [Moved to module level.]
    use PRTGenericMod, only : SetState
    use PRTGenericMod, only : carbon12_element, nitrogen_element,  phosphorus_element
    use PRTGenericMod, only : struct_organ, leaf_organ, fnrt_organ, sapw_organ, store_organ
    use PRTGenericMod, only : repro_organ
    use PRTGenericMod, only : prt_carbon_allom_hyp, prt_cnp_flex_allom_hyp
    
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
    !class(prt_vartypes), pointer :: prt_obj ! PARTEH object to hold elemental states.  [-> Agument]
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
        ! call endrun(msg=errMsg(sourcefile, __LINE__))
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
    ! Why not use 0?  I need to confirm that trim is the same oposite of spread
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
    cohort_obj%laimemory = 0._r8
    cohort_obj%sapwmemory = 0._r8
    cohort_obj%structmemory = 0._r8	 
    !cstatus = leaves_on
    cohort_obj%status_coh = leaves_on

    stem_drop_fraction = EDPftvarcon_inst%phen_stem_drop_fraction(cohort_obj%pft)

     if (EDPftvarcon_inst%season_decid(cohort_obj%pft) == itrue .and. &
         any(site%cstatus == [phen_cstat_nevercold,phen_cstat_iscold])) then
      cohort_obj%laimemory = b_leaf
      cohort_obj%sapwmemory = b_sapw * stem_drop_fraction
      cohort_obj%structmemory = b_struct * stem_drop_fraction	    
      b_leaf  = 0._r8
      b_sapw = (1._r8 - stem_drop_fraction) * b_sapw
      b_struct  = (1._r8 - stem_drop_fraction) * b_struct
      !cstatus = leaves_off
      cohort_obj%status_coh = leaves_off
    endif

    if (EDPftvarcon_inst%stress_decid(cohort_obj%pft) == itrue .and. &
        any(site%dstatus == [phen_dstat_timeoff,phen_dstat_moistoff])) then
      cohort_obj%laimemory = b_leaf
      cohort_obj%sapwmemory = b_sapw * stem_drop_fraction
      cohort_obj%structmemory = b_struct * stem_drop_fraction	    
      b_leaf  = 0._r8
      b_sapw = (1._r8 - stem_drop_fraction) * b_sapw
      b_struct  = (1._r8 - stem_drop_fraction) * b_struct	    
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
    
    ! Initialize the PARTEH object:
    ! prt_obj => null()
    ! call InitPRTObject(prt_obj)
    ! The PARTEH object is now passed in.

    do el = 1,num_elements

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
        m_repro  = 0._r8

      case(nitrogen_element)

        m_struct = b_struct * EDPftvarcon_inst%prt_nitr_stoich_p1(cohort_obj%pft, struct_organ)
        m_leaf   = b_leaf * EDPftvarcon_inst%prt_nitr_stoich_p1(cohort_obj%pft, leaf_organ)
        m_fnrt   = b_fnrt * EDPftvarcon_inst%prt_nitr_stoich_p1(cohort_obj%pft, fnrt_organ)
        m_sapw   = b_sapw * EDPftvarcon_inst%prt_nitr_stoich_p1(cohort_obj%pft, sapw_organ)
        m_store  = b_store * EDPftvarcon_inst%prt_nitr_stoich_p1(cohort_obj%pft, store_organ)
        m_repro  = 0._r8

      case(phosphorus_element)

        m_struct = b_struct * EDPftvarcon_inst%prt_phos_stoich_p1(cohort_obj%pft, struct_organ)
        m_leaf   = b_leaf * EDPftvarcon_inst%prt_phos_stoich_p1(cohort_obj%pft, leaf_organ)
        m_fnrt   = b_fnrt * EDPftvarcon_inst%prt_phos_stoich_p1(cohort_obj%pft, fnrt_organ)
        m_sapw   = b_sapw * EDPftvarcon_inst%prt_phos_stoich_p1(cohort_obj%pft, sapw_organ)
        m_store  = b_store * EDPftvarcon_inst%prt_phos_stoich_p1(cohort_obj%pft, store_organ)
        m_repro  = 0._r8
      end select

      select case(hlm_parteh_mode)
      case (prt_carbon_allom_hyp,prt_cnp_flex_allom_hyp )

        ! Put all of the leaf mass into the first bin
        call SetState(cohort_parteh, leaf_organ, element_id, m_leaf,1)
        do iage = 2,nleafage
          call SetState(cohort_parteh, leaf_organ, element_id, 0._r8, iage)
        end do

        call SetState(cohort_parteh, fnrt_organ, element_id, m_fnrt)
        call SetState(cohort_parteh, sapw_organ, element_id, m_sapw)
        call SetState(cohort_parteh, store_organ, element_id, m_store)
        call SetState(cohort_parteh, struct_organ, element_id, m_struct)
        call SetState(cohort_parteh, repro_organ, element_id, m_repro)

      case default
        write(fates_log(),*) 'Unspecified PARTEH module during inventory intitialization'
        ! call endrun(msg=errMsg(sourcefile, __LINE__))
        call endrun(msg = errMsg(__FILE__, __LINE__))
      end select

      end do

    call cohort_parteh%CheckInitialConditions()
    
    ! End Part 4.

  end subroutine init_temporary_cohort

  !=================================================================================================
  ! Managed Mortality Subroutines:
  !=================================================================================================
  
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
  
  subroutine anthro_disturbance_rate(site_in, bc_in, frac_site_primary)
    ! ----------------------------------------------------------------------------------------------
    ! Calculate the extent of any disturbance resulting from potential human vegetation management
    ! at the current time step.  The disturbance, at a patch level, may or may not subsequently
    ! occur depending on whether it is the dominant disturbance.
    !
    ! This will result in update of an update of:
    ! - The site level harvest_carbon_flux value member.  [Not sure why is this calculated here!]
    ! - The patch level disturbance_rates(dtype_ilog) member.
    ! - The cohort level lmort_direct, lmort_collateral, lmort_infra, and l_degrad members.
    !
    ! This routine is called from EDPatchDynamicsMod.F90: disturbance_rates() in the FATES event
    ! loop. This initial version is based on logging code extracted from there but will be expanded
    ! to include other vegetation management.
    !
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    ! use EDPatchDynamicsMod, only : get_frac_site_primary
    ! use EDMortalityFunctionsMod, only : mortality_rates
    ! use FatesAllometryMod, only : carea_allom
    use EDLoggingMortalityMod, only : get_harvest_rate_area
    use EDLoggingMortalityMod, only : logging_time
    use EDLoggingMortalityMod, only : LoggingMortality_frac
    use EDParamsMod, only : logging_export_frac
    use EDPftvarcon, only : EDPftvarcon_inst
    use EDTypesMod, only : dtype_ilog
    use FatesConstantsMod, only : fates_tiny
    use FatesConstantsMod, only : nearzero
    use FatesLitterMod, only : ncwd
    use PRTGenericMod, only : all_carbon_elements
    use PRTGenericMod, only : sapw_organ
    use PRTGenericMod, only : struct_organ
    use SFParamsMod, only : SF_VAL_CWD_FRAC
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: site_in
    type(bc_in_type), intent(in) :: bc_in
    ! Would like to get this with get_frac_site_primary() but that creates a circular dependency:
    real(r8), intent(in) :: frac_site_primary
    
    ! Locals:
    type (ed_patch_type) , pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort

!     real(r8) :: cmort
!     real(r8) :: bmort
!     real(r8) :: hmort
!     real(r8) :: frmort
!     real(r8) :: smort
!     real(r8) :: asmort

    real(r8) :: lmort_direct
    real(r8) :: lmort_collateral
    real(r8) :: lmort_infra
    real(r8) :: l_degrad         ! fraction of trees that are not killed but suffer from forest 
                                 ! degradation (i.e. they are moved to newly-anthro-disturbed 
                                 ! secondary forest patch)
    real(r8) :: dist_rate_ldist_notharvested
!    integer  :: threshold_sizeclass
!    integer  :: i_dist
!    real(r8) :: frac_site_primary
    real(r8) :: harvest_rate
    
    ! ----------------------------------------------------------------------------------------------
    
    ! The following is the based on logging code extracted from EDPatchDynamicsMod.F90:
    ! disturbance_rates().
    
    ! First calculate the fraction of the site that is primary land:
    ! call get_frac_site_primary(site_in, frac_site_primary)
 
    ! An estimate of the harvest flux will be made and stored:
    site_in%harvest_carbon_flux = 0._r8
    
    ! ----------------------------------------------------------------------------------------------
    ! Calculate Mortality Rates:
    ! ----------------------------------------------------------------------------------------------
    
    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))
      currentCohort => currentPatch%shortest
      
      do while(associated(currentCohort))
        currentCohort%patchptr => currentPatch ! Why is this necessary?????
        
        ! Mortality for trees in the understorey.   [What is this about????]
        !call mortality_rates(currentCohort,bc_in,cmort,hmort,bmort,frmort,smort,asmort)
        !currentCohort%dmort  = cmort+hmort+bmort+frmort+smort+asmort

      ! The crown area is used to calculate other disturbance types as well so we assume it has
      ! has been updated already, since growth was applied.
      !call carea_allom(currentCohort%dbh,currentCohort%n,site_in%spread,currentCohort%pft, &
      !currentCohort%c_area)

! Initialize diagnostic mortality rates
! currentCohort%cmort = cmort
! currentCohort%bmort = bmort
! currentCohort%hmort = hmort
! currentCohort%frmort = frmort
! currentCohort%smort = smort
! currentCohort%asmort = asmort

        call LoggingMortality_frac(currentCohort%pft, currentCohort%dbh, &
                                   currentCohort%canopy_layer, lmort_direct, lmort_collateral, &
                                   lmort_infra, l_degrad, bc_in%hlm_harvest_rates, &
                                   bc_in%hlm_harvest_catnames, bc_in%hlm_harvest_units, &
                                   currentPatch%anthro_disturbance_label, &
                                   currentPatch%age_since_anthro_disturbance, frac_site_primary)

        currentCohort%lmort_direct     = lmort_direct
        currentCohort%lmort_collateral = lmort_collateral
        currentCohort%lmort_infra      = lmort_infra
        currentCohort%l_degrad         = l_degrad

        ! Estimate the wood product (trunk_product_site)
        if (currentCohort%canopy_layer>=1) then
          site_in%harvest_carbon_flux = site_in%harvest_carbon_flux + currentCohort%lmort_direct * &
          currentCohort%n * (currentCohort%prt%GetState(sapw_organ, all_carbon_elements) + &
          currentCohort%prt%GetState(struct_organ, all_carbon_elements)) * &
          EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * SF_val_CWD_frac(ncwd) * &
          logging_export_frac
        endif

        currentCohort => currentCohort%taller
      end do ! Cohort loop.
      
      !currentPatch%disturbance_mode = fates_unset_int
      currentPatch => currentPatch%younger
    end do ! Patch loop.
    
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
    
    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))
      currentPatch%total_canopy_area = 0._r8
      currentCohort => currentPatch%shortest
      do while(associated(currentCohort))   
        if(currentCohort%canopy_layer==1)then
             currentPatch%total_canopy_area = currentPatch%total_canopy_area + currentCohort%c_area
        end if
        currentCohort => currentCohort%taller
      end do
      currentPatch => currentPatch%younger
    end do
    
    ! ----------------------------------------------------------------------------------------------
    ! Calculate Disturbance Rates based on the mortality rates just calculated:
    ! ----------------------------------------------------------------------------------------------
    
    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))   

      !currentPatch%disturbance_rates(dtype_ifall) = 0.0_r8
      currentPatch%disturbance_rates(dtype_ilog)  = 0.0_r8
      !currentPatch%disturbance_rates(dtype_ifire) = 0.0_r8

      dist_rate_ldist_notharvested = 0.0_r8

      currentCohort => currentPatch%shortest
      do while(associated(currentCohort))   

        if(currentCohort%canopy_layer == 1)then

! Treefall Disturbance Rate
! currentPatch%disturbance_rates(dtype_ifall) = currentPatch%disturbance_rates(dtype_ifall) + &
! fates_mortality_disturbance_fraction * &
! min(1.0_r8,currentCohort%dmort)*hlm_freq_day*currentCohort%c_area/currentPatch%area

          ! Logging Disturbance Rate:
          currentPatch%disturbance_rates(dtype_ilog) = &
                              currentPatch%disturbance_rates(dtype_ilog) + &
                              min(1.0_r8, currentCohort%lmort_direct + &
                              currentCohort%lmort_collateral + &
                              currentCohort%lmort_infra + &
                              currentCohort%l_degrad ) * &
                              currentCohort%c_area/currentPatch%area

          ! Non-harvested part of the logging disturbance rate:
          dist_rate_ldist_notharvested = dist_rate_ldist_notharvested + currentCohort%l_degrad * &
                                         currentCohort%c_area/currentPatch%area
        endif
        
        currentCohort => currentCohort%taller
      enddo ! Cohort loop

      ! for non-closed-canopy areas subject to logging, add an additional increment of area disturbed
      ! equivalent to the fradction loged to account for transfer of interstitial ground area to new secondary lands
      if (logging_time .and. &
         (currentPatch%area - currentPatch%total_canopy_area) .gt. fates_tiny ) then
        ! The canopy is NOT closed. 

        call get_harvest_rate_area(currentPatch%anthro_disturbance_label, &
                                   bc_in%hlm_harvest_catnames, bc_in%hlm_harvest_rates, &
                                   frac_site_primary, currentPatch%age_since_anthro_disturbance, &
                                   harvest_rate)

        currentPatch%disturbance_rates(dtype_ilog) = currentPatch%disturbance_rates(dtype_ilog) + &
            (currentPatch%area - currentPatch%total_canopy_area) * harvest_rate / currentPatch%area

        ! Non-harvested part of the logging disturbance rate
        dist_rate_ldist_notharvested = dist_rate_ldist_notharvested + &
                                       (currentPatch%area - currentPatch%total_canopy_area) * &
                                       harvest_rate / currentPatch%area
      endif

      ! fraction of the logging disturbance rate that is non-harvested
      if (currentPatch%disturbance_rates(dtype_ilog) .gt. nearzero) then
        currentPatch%fract_ldist_not_harvested = dist_rate_ldist_notharvested / &
                                                 currentPatch%disturbance_rates(dtype_ilog)
      endif

      currentPatch => currentPatch%younger
    enddo ! Patch loop 
    
  end subroutine anthro_disturbance_rate
  
  !=================================================================================================
  
  function anthro_mortality_rate(cohort, bc_in, frac_site_primary) result(dndt_logging) ! anthro_mortality non-disturbing?????
    ! ----------------------------------------------------------------------------------------------
    ! Calculate mortality resulting from human vegetation management at the cohort level.
    ! Mortality is returned as a change in number (density?) per unit time (year).
    ! This routine is called from EDMortalityFunctionsMod: Mortality_Derivative().
    !
    ! Human induced mortality includes logging but other actives as well such as thinning,
    ! understory clearing, agricultural harvest, etc.
    ! The existing code excludes disturbance inducing mortality resulting from logging.  We need to
    ! decide if and when we want to follow suit.
    !
    ! This subroutine is designed to encapsulate the management specific logic so the calling code
    ! does not have to be aware of it.
    ! This code currently extracts the logging specific code from EDMortalityFunctionsMod.F90:
    ! Mortality_Derivative() with only formating and name changes.
    ! However, there is potential problem with that code that I'm working to resolve.
    ! Logic for more activities will follow.
    !
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    !use EDPatchDynamicsMod, only : get_frac_site_primary
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
    else
       ! Mortality from logging in the canopy is ONLY disturbance generating, don't
       ! update number densities via non-disturbance inducing death
       dndt_logging = 0.0_r8
    endif
    
  end function anthro_mortality_rate

  !=================================================================================================

  subroutine management_fluxes(current_site, current_patch, new_patch, patch_site_areadis)
    ! ----------------------------------------------------------------------------------------------
    ! Perform the the fluxes resulting from all the management activities performed during this
    ! timestep.
    ! Called from EDPatchDynamicsMod.F90: spawn_patches().
    !
    ! Currently this is just a wrapper for the existing logging code and will be expanded later.
    ! 
    ! ----------------------------------------------------------------------------------------------
    
    ! Uses:
    use EDLoggingMortalityMod, only : logging_litter_fluxes
    use EDTypesMod, only : dtype_ilog
    
    ! Arguments:
    type(ed_site_type), intent(inout), target :: current_site ! Possibly unnecessary, see below.
    type (ed_patch_type), intent(inout), target :: current_patch ! patch_in?
    type (ed_patch_type), intent(inout), target :: new_patch ! The...
    ! This can be calculated from data in the cohort so doesn't really need to be passed in:
    real(r8), intent(in) :: patch_site_areadis ! total area disturbed in m2 per patch per day
    
    ! Locals:
    ! Do we have to pass in the site to modify it or can we just get it from patch_in?
    !type (ed_site_type), pointer :: current_site
    
    ! ----------------------------------------------------------------------------------------------
    
    ! Checking the disturbance mode here rather than in the calling code allows for more flexibility
    ! but is not really necessary now:
    if (current_patch%disturbance_mode .eq. dtype_ilog) then
      call logging_litter_fluxes(current_site, current_patch, new_patch, patch_site_areadis)
    endif
    
  end subroutine management_fluxes

  !=================================================================================================
  ! Patch and Cohort Search Subroutines:
  !=================================================================================================


end module FatesVegetationManagementMod
