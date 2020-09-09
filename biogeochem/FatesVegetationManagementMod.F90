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

end module FatesVegetationManagementMod
