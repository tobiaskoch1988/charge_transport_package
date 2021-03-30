!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! START module Molecule_class  !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Molecule_class
	type subhoppingsite_t
		integer           ::id
		real, dimension(3)::koord		
	end type subhoppingsite_t

	type hoppingsite_t
		real, dimension(3)                    :: COM
		integer				                  :: N_subhoppingsites
		type (subhoppingsite_t), allocatable  :: subhoppingsite(:)
		logical:: hasCOM=.false.
		logical:: hasSubhoppingsite=.false.
	end type hoppingsite_t

	type Molecule_t	
		!! counts all the molecules in the system
		integer :: molid 
		!! resid in gro format
		integer :: resid
		!! number of atoms in the molecule
		integer :: NAtoms
		character(len=5)   :: residuename
		character(len=500) :: molname
		!! koordinates (1:NAtoms,1:3)
		real, 	  allocatable, dimension(:,:)             :: koordinates
		!! elementnames (1:NAtoms) in the molecule
		character(len=5), allocatable, dimension(:)       :: elementnames	
		!! center of mass of the molecule
		real, dimension(3)   :: COM
		logical::hasCOM=.false.
		type(hoppingsite_t)  :: hoppingsite 
		logical              :: hasHoppingsite=.false.
	end type Molecule_t

	type COM_t
		REAL, DIMENSION(3) :: COM
		REAL :: COM_x,COM_y,COM_z
		logical::hasCOM=.false.
	end type COM_t

	type atom_t
		REAL, DIMENSION(3)   :: koordinate
		CHARACTER(Len=5)     :: elementname
		REAL:: mass
	end type

	type hoppingpair_t
		integer::pairid
		integer::molAid
		integer::molBid
		integer::residA,residB
		real, dimension(3) ::dist_vec
		real::norm_dist_vec
		logical ::has_dist_vec=.false.
		real::S_AB_h,S_AB_e
		real::J_AB_h,J_AB_e
		real::rateAB_h,rateAB_e,rateBA_h,rateBA_e
		real::l_out_h,l_out_e
		real::E_out_h,E_out_e
		real::l_out_h_32mer,l_out_e_32mer
		real::E_out_h_32mer,E_out_e_32mer		
		real::theta_intra
		character(len=10)::seq	
		logical::has_rates_h=.false.,has_rates_e=.false.,has_intramolecular=.false.

	end type hoppingpair_t

	type system_t
		type (Molecule_t),     allocatable :: molecules(:)
		type (hoppingpair_t), allocatable :: hoppingpairs(:)
		integer::N_Resids,N_molecules,N_hopping_pairs
		integer::N_atoms_all
		integer::N_DIPBI,N_PPDI,N_P3HT,N_PBDT_TS1,N_Species
	end type system_t



	contains
	
!!!!!!!!!!!!!!!! SYSTEM START
	subroutine init_system(system)
		type (system_t) :: system  
		type (Molecule_t), allocatable :: molecules(:)
		type (hoppingpair_t), allocatable :: hoppingpairs(:)
		integer::N_Resids,N_molecules,N_hopping_pairs
		integer::N_atoms_all
		integer::N_DIPBI,N_PPDI,N_P3HT,N_PBDT_TS1

		system%N_Resids=0
		system%N_molecules=0
		system%N_atoms_all=0
		system%N_DIPBI=0
		system%N_PPDI=0
		system%N_P3HT=0
		system%N_PBDT_TS1=0
		system%N_Species=0
		system%N_hopping_pairs=0
		
	end subroutine init_system


	subroutine make_new_molecule_for_system(system,residA,NAtomsA,resnameA,molnameA,koordA,elementsA)
		implicit none
		type (system_t)   		 :: system 
		!type (Molecule_t), allocatable  :: molecules(:)
		type (Molecule_t) 		 :: newmolecule
		!! newmolecule with data for molecule A
		integer 			      , intent(in) :: residA,NAtomsA
		character(LEN=*)                      , intent(in) :: resnameA,molnameA
		character(Len=5), dimension(NAtomsA)  , intent(in) :: elementsA
		real, dimension(NAtomsA,3) 	          , intent(in) :: koordA
		integer::molidA
		!!call null_Molecule(newmolecule)
		
		!! increase molid in system
		IF(.not. system_has_molecules(system)) THEN
				molidA=1	
		ELSE	
				molidA=system%N_molecules+1
		END IF
				
		!! initialize newmolecule with data for molecule A
		call init_Molecule(newmolecule,molidA,residA,NAtomsA,resnameA,molnameA,koordA,elementsA)		

		!!! increase number of molecules
		IF( INDEX(resnameA,'DIPBI') /= 0) THEN
			system%N_DIPBI=system%N_DIPBI+1
		ELSE IF( INDEX(resnameA,'PPDI') /= 0) THEN
			system%N_PPDI=system%N_PPDI+1
		ELSE IF( (INDEX(resnameA,'P3') /= 0) .OR. (INDEX(resnameA,'PH') /= 0) .OR. (INDEX(resnameA,'THP') /= 0)) THEN
			system%N_P3HT=system%N_P3HT+1	
		ELSE IF ( (INDEX(resnameA,'8poly') /= 0) .OR. (INDEX(resnameA(1:1),'A') /= 0) &
							&.OR. (INDEX(resnameA(1:1),'B') /= 0) )THEN
			system%N_PBDT_TS1=system%N_PBDT_TS1+1
		ELSE 
			!!! other species, not defined yet / default
			system%N_Species=system%N_Species+1
		ENDIF
		!! general counter
		system%N_Resids=system%N_Resids+1
		system%N_molecules=system%N_molecules+1

		IF(.not. system_has_molecules(system)) THEN
			allocate(system%molecules(1))
			system%molecules(1)=newmolecule
		ELSE
			call append_molecule_to_system(system%molecules,newmolecule)
		ENDIF
			
		system%N_molecules=system%N_molecules+1

	end subroutine make_new_molecule_for_system	


	subroutine append_molecule_to_system(system,newmolecule)
		!! appends the newmolecule to the list of molecules in system
		implicit none
		type (Molecule_t), intent(in)::newmolecule
		type (Molecule_t), ALLOCATABLE, INTENT(INOUT) :: system(:)
		type (Molecule_t), ALLOCATABLE :: new_system(:) 
		integer::nsys

		nsys=size(system)
		allocate(new_system(nsys+1))

		new_system(1:nsys)=system(1:nsys)

		new_system(nsys+1) = newmolecule

		deallocate(system)
		allocate(system(nsys+1))
		!! transfer molecule
		system(:) = new_system(:)

		deallocate(new_system)
 	end subroutine append_molecule_to_system


	!!! prints the coordinates of a molecule "this" in xyz format on the screen, e.g. for debugging
	subroutine print_molecule_xyz(this)
			implicit none
			type (Molecule_t), intent(in)::this
			integer::i
			!Character(len=*), PARAMETER::form2='(2X,A15,1X,A2,1X,3F16.6,1X,A2)'
			WRITE(*,*) this%NAtoms
			WRITE(*,*)
			do i=1,this%NAtoms
				write(*,*) this%elementnames(i),this%koordinates(i,1),this%koordinates(i,2),this%koordinates(i,3)
			end do
 	end subroutine print_molecule_xyz


   	!! Return .true. if the system has a molecules.
   	pure function system_has_molecules( this ) result( hasmolecules )
	      	!! Input: molecule 
	      	type(system_t), intent(in) :: this
	      	!! Output: whether the molecule has molecules
	      	logical                      :: hasmolecules
	      	!!
	      	hasmolecules = Allocated( this%molecules )
   	end function system_has_molecules


   	!! Return .true. if the system has a hoppingpair.
   	pure function system_has_hoppingpairs( this ) result( hashoppingpairs )
	      	!! Input: molecule 
	      	type(system_t), intent(in) :: this
	      	!! Output: whether the molecule has hoppingpairs
	      	logical                      :: hashoppingpairs
	      	!!
	      	hashoppingpairs = Allocated( this%hoppingpairs )
   	end function system_has_hoppingpairs


!!!!!!!!!!!!!!!! SYSTEM END


!!!!!!!!!!!!!!!!! hopping_pair START
	subroutine init_hoppingpair_in_system(this,system)
		implicit none
		type(hoppingpair_t), intent(out) :: this
		type(system_t), intent(inout) :: system
		integer::pairid
		integer::molAid
		integer::molBid
		integer::residA,residB
		real, dimension(3) ::dist_vec
		real::norm_dist_vec
		logical ::has_dist_vec=.false.
		real::S_AB_h,S_AB_e
		real::J_AB_h,J_AB_e
		real::rateAB_h,rateAB_e,rateBA_h,rateBA_e
		real::l_out_h,l_out_e
		real::E_out_h,E_out_e
		real::l_out_h_32mer,l_out_e_32mer
		real::E_out_h_32mer,E_out_e_32mer	
		real::theta_intra
		character(len=10)::seq	
		logical::has_rates_h=.false.,has_rates_e=.false.,has_intramolecular=.false.
		
		this%pairid=system%N_hopping_pairs+1
		system%N_hopping_pairs=system%N_hopping_pairs+1
		
		this%molAid=0
		this%molBid=0
		this%residA=0
		this%residB=0
		this%dist_vec=(/0.0,0.0,0.0/)
		this%norm_dist_vec=0.0
		this%has_dist_vec=.false.
		this%S_AB_h=0.0
		this%S_AB_e=0.0
		this%J_AB_h=0.0
		this%J_AB_e=0.0
		this%rateAB_h=0.0
		this%rateAB_e=0.0
		this%rateBA_h=0.0
		this%rateBA_e=0.0
		this%l_out_h=0.0
		this%l_out_e=0.0
		this%E_out_h=0.0
		this%E_out_e=0.0
		this%l_out_h_32mer=0.0
		this%l_out_e_32mer=0.0
		this%E_out_h_32mer=0.0
		this%E_out_e_32mer=0.0		
		this%has_rates_h=.false.
		this%has_rates_e=.false.
		this%has_intramolecular=.false.
		this%theta_intra=0.0
		this%seq=''
	end subroutine init_hoppingpair_in_system

!!!!!!!!!!!!!!!!! hopping_pair END



!!!!!!!!!!!!!! SUBHOPPINGSITE START
	subroutine make_new_subhoppingsite_for_hoppingsite(hoppingsite,koordA)
		implicit none
		type (hoppingsite_t) ,  intent(inout) ::hoppingsite				
		real, dimension(3),     intent(in)    ::koordA
		type (subhoppingsite_t)  	          ::newsubhoppingsite
		integer::id

		id=hoppingsite%N_subhoppingsites+1
		call set_subhoppingsite(newsubhoppingsite,id,koordA(1),koordA(2),koordA(3))
		
		IF(.not. Hoppingsite_has_SubHoppingsite(hoppingsite))THEN
			allocate(hoppingsite%subhoppingsite(1))
			hoppingsite%subhoppingsite(1)=newsubhoppingsite
		ELSE
			call append_append_subhoppingsite_to_hoppingsitelist(hoppingsite%subhoppingsite,newsubhoppingsite)
		END IF

		IF(hoppingsite%N_subhoppingsites .le. 1) THEN
			hoppingsite%hasSubhoppingsite=.true.
		END IF
		hoppingsite%N_subhoppingsites=hoppingsite%N_subhoppingsites+1
	end subroutine make_new_subhoppingsite_for_hoppingsite


	subroutine append_append_subhoppingsite_to_hoppingsitelist(hoppingsitelist,newsubhoppingsite)
		!! appends the newhoppingsite to the list of hoppingsites in molecule
		implicit none
		type (subhoppingsite_t), intent(in)::newsubhoppingsite
		type (subhoppingsite_t), allocatable, intent(inout) :: hoppingsitelist(:)
		type (subhoppingsite_t), allocatable :: new_hoppingsitelist(:) 
		integer::nsys

		nsys=size(hoppingsitelist)
		allocate(new_hoppingsitelist(nsys+1))

		new_hoppingsitelist(1:nsys)=hoppingsitelist(1:nsys)

		new_hoppingsitelist(nsys+1) = newsubhoppingsite

		deallocate(hoppingsitelist)
		allocate(hoppingsitelist(nsys+1))
		!! transfer newsubhoppingsite
		hoppingsitelist(:) = new_hoppingsitelist(:)

		deallocate(new_hoppingsitelist)
 	end subroutine append_append_subhoppingsite_to_hoppingsitelist


	subroutine set_subhoppingsite(this,id,koordx,koordy,koordz)
		type(subhoppingsite_t), intent(inout)   :: this
		integer 	            :: id
		real, intent(in)        ::koordx,koordy,koordz
		this%id=id
		this%koord(1)=koordx
		this%koord(2)=koordy
		this%koord(3)=koordz
	end subroutine set_subhoppingsite


   	!! Return .true. if the Hoppingsite has a SubHoppingsite.
   	pure function Hoppingsite_has_SubHoppingsite( this ) result( hasSubHoppingsite )
	      	!! Input: Hoppingsite 
	      	type(Hoppingsite_t), intent(in) :: this
	      	!! Output: whether the molecule has a gradient
	      	logical                      :: hasSubHoppingsite
	      	!!
	      	hasSubHoppingsite = Allocated( this%subhoppingsite )
   	end function Hoppingsite_has_SubHoppingsite

	subroutine set_COM_in_Hoppingsite(this,COMx,COMy,COMz)
		IMPLICIT NONE
		type(hoppingsite_t), intent(inout) :: this
		REAL::COMx,COMy,COMz
		this%COM(1)=COMx
		this%COM(2)=COMy
		this%COM(3)=COMz
		this%hasCOM=.true.
	end subroutine set_COM_in_Hoppingsite
!!!!!!!!!!!!!! SUBHOPPINGSITE END	


!!!!!!!!!!! ATOM START 
	subroutine null_atom(this)
		IMPLICIT NONE
		type (atom_t), intent (out) :: this
		REAL, DIMENSION(3)   :: koordinate
		CHARACTER(Len=5)     :: elementname
		REAL:: mass
		this%elementname=' '
		this%koordinate=(/ 0.0, 0.0, 0.0 /)
		this%mass=0.0
	end subroutine null_atom

	subroutine init_atom(this,elementname,koordx,koordy,koordz)
		IMPLICIT NONE
		type (atom_t), intent (out) :: this
		CHARACTER(Len=*)          :: elementname
		REAL, INTENT(IN)::koordx,koordy,koordz
		
		this%elementname=TRIM(elementname)
		this%koordinate=(/ koordx, koordy, koordz/)
		this%mass=0.0
	end subroutine 

!!!!!!!!!!! ATOM END


!!!!!!!!!  MOLECULE START

   	!! Return .true. if the molecule has a koordinates.
   	pure function Mol_has_koordinates( this ) result( haskoordinates )
	      	!! Input: molecule 
	      	type(Molecule_t), intent(in) :: this
	      	!! Output: whether the molecule has a gradient
	      	logical                      :: haskoordinates
	      	!!
	      	haskoordinates = Allocated( this%koordinates )
   	end function Mol_has_koordinates

   	!! Return .true. if the molecule has a elementnames.
   	pure function Mol_has_elementnames( this ) result( haselementnames )
	      	!! Input: molecule 
	      	type(Molecule_t), intent(in) :: this
	      	!! Output: whether the molecule has a gradient
	      	logical                      :: haselementnames
	      	!!
	      	haselementnames = Allocated( this%elementnames )
   	end function Mol_has_elementnames


	subroutine init_Molecule(this,molidA,residA,NAtomsA,resnameA,molnameA,koordA,elementsA)
		IMPLICIT NONE
		! Constructor
		type (Molecule_t), intent (out) :: this
		integer, intent (in) :: molidA,residA,NAtomsA
		character(LEN=*)                      , intent(in) :: resnameA,molnameA
		character(Len=5), dimension(NAtomsA)  , intent(in) :: elementsA
		real, dimension(NAtomsA,3) 	      , intent(in) :: koordA
		integer::i
			this%molid  = molidA
			this%resid  = residA
			this%NAtoms = NAtomsA
			this%residuename=TRIM(resnameA)
			this%molname=TRIM(molnameA)
			allocate(this%koordinates(this%NAtoms,3))
			this%koordinates=0
			allocate(this%elementnames(this%NAtoms))
			this%elementnames=' '
			DO i=1,this%NAtoms
				this%koordinates(i,1:3)=koordA(i,1:3)
				this%elementnames(i)=elementsA(i)
			END DO
			this%COM=(/ 0.0, 0.0, 0.0 /)
	end subroutine init_Molecule


	function get_resid_Molecule(this) result(resid)
		IMPLICIT NONE
		type (Molecule_t), intent (in) :: this
		integer :: resid
		resid = this%resid
	end function get_resid_Molecule

	function get_NAtoms_Molecule(this) result(NAtoms)
		IMPLICIT NONE
		type (Molecule_t), intent (in) :: this
		integer :: NAtoms
		NAtoms = this%NAtoms
	end function get_NAtoms_Molecule

	subroutine set_COM_in_Molecule(this,COMx,COMy,COMz)
		IMPLICIT NONE
		type(Molecule_t), intent(inout) :: this
		REAL::COMx,COMy,COMz
		this%COM(1)=COMx
		this%COM(2)=COMy
		this%COM(3)=COMz
		this%hasCOM=.true.
	end subroutine set_COM_in_Molecule


	subroutine set_COM(this,COMx,COMy,COMz)
		IMPLICIT NONE
		type(COM_t), intent(inout) :: this
		REAL::COMx,COMy,COMz
		this%COM_x=COMx
		this%COM_y=COMy
		this%COM_z=COMz
		this%hasCOM=.true.
	end subroutine set_COM
	
	subroutine print_Molecule(this,printresid)
		IMPLICIT NONE
		type (Molecule_t), intent (in) :: this
		logical, optional, intent (in) :: printresid
		if (present(printresid)) then
			if (printresid) write (*,'(i2,a2)',advance='no') this%resid,': '

		endif
		WRITE(*,*)  this%resid,TRIM(this%molname)
	end subroutine print_Molecule

	subroutine strcpy(s,c)
		IMPLICIT NONE
		character, dimension (:), intent (out) :: s
		character*(*), intent (in) :: c
		integer::i
		do i = 1, max(size(s),len(c))
			s(i) = c(i:i)
		enddo
	end subroutine strcpy

end module Molecule_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! END Molecule_class  !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! START  enumerated_gro_to_VOTCA__DIPBI_P3HT_prepare  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main
	use Molecule_class
	IMPLICIT NONE
	CALL enumerated_gro_to_VOTCA__DIPBI_P3HT_prepare()
CONTAINS

!! PROGRAM to Setup the inputfiles for a votca KMC simulation on DIPBI/P3HT blend
!! Reads gro file, which should be preprocessed to contain the information about the divition of the polymer chain into sites, via adapted resids.
!! Reads Data for the parametization of segments from a lambda_in_all_inputfile. (FORMAT differs!)
!! Adds H atoms to the molecules from an FF to QM representation.
!! Partitioning in DIPBI_K (with chain) P3HT/P3MT with a START segment S, middle M, or terminal segment E. 
!! (eg. Converntion PH06S labels a Polythiophine with 6 units at the chain start. PM01M is a P3MT with a single ring in the middle of a chain.)
!! creates a VOTCA mapping file
!! creates a VOTCA fake_topology file
!! molpol_optionsfiles
!! Can shift the resid index form 0 to 1 to suite VOTCA convention.
!! creates a new gro_outputfile with the new nomenklature
SUBROUTINE enumerated_gro_to_VOTCA__DIPBI_P3HT_prepare()
type (system_t) :: system
!! Einlesen aus der gro_inputfile.gro Datei und Verschieben der Resid um I_shift in gro_outputfile.gro 
Character(5), ALLOCATABLE, Dimension(:)    :: atomsorteA,residue_names,residue_name_list
Real, ALLOCATABLE, DIMENSION(:,:)          :: koordA,vel
!! reorganization energy
Real*8, ALLOCATABLE, DIMENSION(:,:)        ::lambda_array,lambda_array_all
!! Resid und Index , Resid und anzahl der Schwefelatome
INTEGER, ALLOCATABLE, DIMENSION(:,:)     :: R_u_I,R_u_S 
INTEGER, ALLOCATABLE, DIMENSION(:)       :: residue_numbers,residue_numbersH,atom_numbers,atom_numbersH,iHOMO_list
Character(5000)             :: gro_inputfile,gro_outputfile,top_all_outputfile,map_filename,neighbourlist_constrains_filename
Character(5000)             :: lambda_in_all_inputfile,votca_xml_options_filename
INTEGER ::N_Resids,Gesamtanzahl
REAL    ::x_box,y_box,z_box
Character(1000) ::dummyname,method,line
CHARACTER(10000)::bashline
Character(5) :: residue_name,lio
Character(3) :: atom
INTEGER::i,j,k,k_ind,l,ierror,residue_number,residue_number_old,atom_number,start_index,end_index,N_seg_types_all
REAL::x,y,z,v_x,v_y,v_z
INTEGER::I_shift
LOGICAL::DIPBI_P3HT_calc,DIPBI_P3MT_calc
LOGICAL::Datei_vorhanden,DEBUG=.false.,modify_mol_h_atoms,KETTE,print_coordH_to_xyz,nm_to_Ang,use_I_shift=.false.
LOGICAL::create_mapping_file=.false.,create_fake_topology=.false.,orbs_ZINDO,create_votca_boxfile,read_lambda_files=.false.
LOGICAL::add_virtual_H_atoms_to_gro=.false.,exclude_virtual_H_atoms=.false.
LOGICAL::create_molpol_optionsfiles=.false.,create_gromp_mdp_file=.false.
LOGICAL::create_geo_orbs_mps_g09=.false.
!!!! MODIFIY H ATOMS 
INTEGER ::ind1,GesamtanzahlH,GesamtanzahlH_precalc,Natoms1H,R_offset,i_offset ! index offset
REAL, ALLOCATABLE, DIMENSION(:,:) 	:: koordH,koord1H,velH,q_n_el_lo
INTEGER, ALLOCATABLE, Dimension(:,:)    :: R_u_I_H ! Resid,Index_fuer_molekuelstart
Character(5), ALLOCATABLE, Dimension(:) :: atomsorteH,atomsorte1H,residue_namesH,elementsH
Character(500), ALLOCATABLE, Dimension(:) :: MolnameH
INTEGER::N_DIPBI,N_P3HT,N_P3HT_segments,N_Species,N_Schwefel
INTEGER::i_start,i_ende,iHOMO,N_types_lambda
Character(5)::residue_namesH_dummy
Character(500) :: Molname,extra_files_modus
LOGICAL, ALLOCATABLE, DIMENSION (:) ::molecule_da
!!! Variables PPDI, PBDT-TS1
LOGICAL:: PPDI_PBDT_TS1_calc,PPDI_Propyl=.false.,mol_start_da,mol_ende_da
!! count the subunits
INTEGER::N_CarboxyS2,N_PBDT,N_PPDI,N_PBDT_TS1,N_PBDT_TS1_segments,N_molecules_acceptor
!! Intramolecular data files for J_AB and dE_in fuer AB und BA Reihenfolge der Segmente
Character(LEN=500)::PBDT_TS1_AB_dE_intr_in_el_file,PBDT_TS1_BA_dE_intr_in_el_file
Character(LEN=500)::PBDT_TS1_AB_dE_intr_in_lo_file,PBDT_TS1_BA_dE_intr_in_lo_file
Character(LEN=500)::PBDT_TS1_AB_J_intra_in_el_file,PBDT_TS1_BA_J_intra_in_el_file
Character(LEN=500)::PBDT_TS1_AB_J_intra_in_lo_file,PBDT_TS1_BA_J_intra_in_lo_file
Character(LEN=1) :: segtype,fin
CHARACTER(LEN=50):: sequence 
CHARACTER(LEN=50), ALLOCATABLE, DIMENSION(:)::sequence_list
Character(5),ALLOCATABLE, DIMENSION(:)::residue_name_site
!!
INTEGER::trajstep,Numcharged,NCPUS,Mem
!! neighbourlist daten
INTEGER, ALLOCATABLE, Dimension(:,:) ::neighbour_list_Res
INTEGER                              ::N_Neighbours
Character(LEN=500)::neighbourlist_filename
LOGICAL::use_min_d_COM_vector_list=.false.


!! SECTION WHERE THE USER CAN SELECT THE OPTIONS APPLIED IN THE SCRIPT
!! Modify for user
DEBUG=.false.
KETTE=.true.
!! System
DIPBI_P3HT_calc=.true.
DIPBI_P3MT_calc=.false.
PPDI_PBDT_TS1_calc=.false.


modify_mol_h_atoms=.true.
create_mapping_file=.true.
create_fake_topology=.true.
create_votca_boxfile=.true.
create_molpol_optionsfiles=.false.
create_gromp_mdp_file=.true.
print_coordH_to_xyz=.true.
nm_to_Ang=.true.
!! Use the HOMO and LUMO levels for ZINDO orbitals 
!! Adapt this, if other J_AB than VOTCA MOO/iZINDO method is used
orbs_ZINDO=.true.
read_lambda_files=.true.
!! Adds H atoms at the intrachain P3HT junctions between neighbouring segments, for the treatment of each segement as an entire molecule / site 
!! This is relevant for the concordance of the number of atoms in the QM files in *.mps files.
add_virtual_H_atoms_to_gro=.true.
!! exclude_virtual_H_atoms=.true.   sets virtual H atoms to 1 in <virtual_mps> section => Excludes H in Emultipole calculation
!! exclude_virtual_H_atoms=.false.  sets virtual H atoms to 0 in <virtual_mps> section => Includes H in Emultipole calculation
exclude_virtual_H_atoms=.true.


!!  make_geo_orbs_mps_g09 
!! extra_files_modus = [false, all, make_g09_for_mps, make_ZINDO, make_xyz, make_folders]    applied in "make_geo_orbs_mps_g09"
!! make_g09_for_mps: creates g09 inputfiles for the geometries to evaluate with data for *.mps files 
!! make_ZINDO: creates g09 input files for the ZINDO calculations, which are needed for *.orbs
!! make_xyz: creates name1_n.xyz for the QC_FILE folder as needed for VOTCA calculations.
!! make_folders: creates folders QC_FILES and MP_FILES, if they are not present in the current directory.
create_geo_orbs_mps_g09=.true.
extra_files_modus="false"


N_Neighbours=12
use_min_d_COM_vector_list=.false.

!!! ggf modifizieren Datei mit allen Daten in eV
!! Format Residue_name  OptS0_n_monomerA	R2=OptS0_el_monomerA	R3=OptS0_lo_monomerA	R7=SCF_el_OptS0_n_monomerA	R8=SCF_lo_OptS0_n_monomerA	R9=SCF_n_OptS0_el_monomerA	R10=SCF_n_OptS0_lo_monomerA
IF(DIPBI_P3HT_calc) THEN
	lambda_in_all_inputfile='data/lambda_in_all_DIPBI_P3MT_B3LYP_6_311Gss_eV.dat'
	votca_xml_options_filename='../options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml'
	neighbourlist_filename='data/new_sorted_neighbours_equi_2000ps_500k_def2_theta_75_G0.ngh'
ELSE IF (PPDI_PBDT_TS1_calc) THEN
	lambda_in_all_inputfile='data/lambda_in_all_PPDI_PBDT-TS1_B3LYP_6_31Gs_eV_V1.dat'
	votca_xml_options_filename='options_cpt_run_PPDI_P3DT_TS1.xml'
	neighbourlist_filename='sim_data_12N/new_sorted_neighbours_mischbox_900k_30ns-34ns_60deg_50NN.ngh'
END IF




!!! change the default setting
trajstep=0
Numcharged=0
NCPUS=72
Mem=64     !in GB



!!! Inicialisation
!!!  Shift
I_shift=0
use_I_shift=.false.
IF (command_argument_count() >= 4) THEN	                
        CALL get_command_argument(1,method)
        IF (TRIM(method) == 'gro_to_VOTCA') THEN
		CONTINUE 
	ELSE
		WRITE(*,*) 'USE: gro_to_VOTCA'
		CALL EXIT (1)
	END IF
        CALL get_command_argument(2,gro_inputfile)
        inquire(file=TRIM(gro_inputfile),exist=Datei_vorhanden) ! Abfrage ob gro-Datei vorhanden ist.
	IF (check_filetermination(gro_inputfile,'.gro')) THEN
        	IF( .NOT. Datei_vorhanden) THEN 
            		WRITE(*,*) ' Fehler: Die_gro_Datei_ist_nicht_vorhanden. '//TRIM(gro_inputfile)
            		WRITE(*,*) ' ENDE'
            		CALL EXIT (1)
        	END IF ! gro-file vorhanden  
	ELSE
		WRITE(*,*) 'Fehler: Datei_hat_keine_gro_Endung'//TRIM(gro_inputfile)
		CALL EXIT (1)
	END IF ! check if *.gro
        CALL get_command_argument(3,gro_outputfile)
        inquire(file=TRIM(gro_outputfile),exist=Datei_vorhanden) ! Abfrage ob gro-Datei vorhanden ist.
	IF (check_filetermination(gro_outputfile,'.gro')) THEN
        	IF(Datei_vorhanden) THEN 
            		WRITE(*,*) ' Fehler: Die_gro_Datei_is_bereits_vorhanden. '//TRIM(gro_outputfile)
			WRITE(*,*) ' Stellen_Sie_sicher_dass_diese_nicht_ueberschrieben_wird! '
            		WRITE(*,*) ' ENDE'
            		CALL EXIT (1)
        	END IF ! gro-file vorhanden  
	ELSE
		WRITE(*,*) 'Fehler: Datei_hat_keine_gro_Endung'//TRIM(gro_outputfile)
		CALL EXIT (1)
	END IF ! check if *.gro
	    CALL get_command_argument(4,dummyname)  
		READ(dummyname,'(I10)') I_shift       
		WRITE(*,*) TRIM(dummyname),' I_shift:',I_shift
 		IF( command_argument_count() .ge. 4) THEN
			use_I_shift=.true.
		ELSE
			use_I_shift=.false.
		END IF
		
		!!! read neighbourlist_filename  
		IF( command_argument_count() >= 5) THEN
			CALL get_command_argument(5,neighbourlist_filename)
			IF (check_filetermination(neighbourlist_filename,'.ngh')) THEN
					IF(Datei_vorhanden) THEN 
							WRITE(*,*) ' Fehler: Die_gro_Datei_is_bereits_vorhanden. '//TRIM(neighbourlist_filename)
					WRITE(*,*) ' Stellen_Sie_sicher_dass_diese_nicht_ueberschrieben_wird! '
							WRITE(*,*) ' ENDE'
							CALL EXIT (1)
					END IF ! neighbourlist_filename vorhanden  
			ELSE
				WRITE(*,*) 'Fehler: Datei_hat_keine_ngh_Endung'//TRIM(neighbourlist_filename)
				CALL EXIT (1)
			END IF ! check if *.ngh	
		END IF !check_filetermination
		
ELSE
	WRITE(*,*) ' mod_resid_index_gro um den index der Resids in der gro datei zu verschieben um I_shift. '
	WRITE(*,*) ' '
	WRITE(*,*) '      mod_resid_index_gro:        method     enumerated_inputfilename.gro    outputfilename.gro  I_shift   '&
	&//'    [ neighbourlist_filename.ngh ] ' 
	WRITE(*,*) ' Use: mod_resid_index_gro   gro_to_VOTCA                	  in.gro               out.gro         1        '&
	&//'    new_sorted_neighbourlist.ngh'
	CALL EXIT (0)
END IF

WRITE(*,*) 'Use a RESID shift: ',I_shift


top_all_outputfile='fake_topology_'//TRIM(gro_inputfile(:LEN_TRIM(gro_inputfile)-4))//'.top'
inquire(file=TRIM(top_all_outputfile),exist=Datei_vorhanden) ! Abfrage ob gro-Datei vorhanden ist.
IF(Datei_vorhanden) THEN 
            		WRITE(*,*) ' Fehler: Die_top_Datei_is_bereits_vorhanden. '//TRIM(top_all_outputfile)
			WRITE(*,*) ' Stellen_Sie_sicher_dass_diese_nicht_ueberschrieben_wird! '
            		WRITE(*,*) ' ENDE'
            		CALL EXIT (1)
END IF ! top-file vorhanden  


map_filename='map_'//TRIM(gro_inputfile(:LEN_TRIM(gro_inputfile)-4))//'.xml'
inquire(file=TRIM(map_filename),exist=Datei_vorhanden) ! Abfrage ob gro-Datei vorhanden ist.
IF(Datei_vorhanden) THEN 
            		WRITE(*,*) ' Fehler: Die_VOTCA_mapping_Datei_is_bereits_vorhanden. '//TRIM(map_filename)
			WRITE(*,*) ' Stellen_Sie_sicher_dass_diese_nicht_ueberschrieben_wird! '
            		WRITE(*,*) ' ENDE'
            		CALL EXIT (1)
END IF ! map-file vorhanden  

neighbourlist_constrains_filename='neighbours_constrained_VOTCA__'//TRIM(gro_inputfile(:LEN_TRIM(gro_inputfile)-4))//'.xml'
inquire(file=TRIM(neighbourlist_constrains_filename),exist=Datei_vorhanden) ! Abfrage ob gro-Datei vorhanden ist.
IF(Datei_vorhanden) THEN 
            		WRITE(*,*) ' Fehler: Die_neighbourlist_constrains_file_is_bereits_vorhanden. '&
						&//TRIM(neighbourlist_constrains_filename)
			WRITE(*,*) ' Stellen_Sie_sicher_dass_diese_nicht_ueberschrieben_wird! '
            		WRITE(*,*) ' ENDE'
            		CALL EXIT (1)
END IF ! neighbourlist_constrains_filename-file vorhanden  



IF(read_lambda_files) THEN
	WRITE(*,*) 'use: '//TRIM(lambda_in_all_inputfile)
END IF
WRITE(*,*) 'uses    gro: ',TRIM(gro_inputfile)
WRITE(*,*) 'creates gro: ',TRIM(gro_outputfile)
WRITE(*,*) 'creates topology: ',TRIM(top_all_outputfile)
WRITE(*,*) 'creates map: ',TRIM(map_filename)



! Einlesen
OPEN(UNIT=17,FILE=TRIM(gro_inputfile),STATUS='OLD',IOSTAT=ierror)
Read(17,*) dummyname                             !Einlesen des Infokopfes im File
Read(17,*) Gesamtanzahl                          !Einlesen der Atomanzahl im File
    !!! Erstelle String lio=Durcheinander, um die Rechnungen mit Hilfe von Zufallszahlen einzigartig zu machen, sodass keine Interferenzen mit anderen Rechnungen ggf. im gleichen Ordner auftreten können. 
    lio=TRIM(make_lio())
    !! Einlesen der Anzahl der Resids=N_Mol aus gro-file
    Bashline="start_index=$(sed -n '3p' "//TRIM(gro_inputfile)//' | cut -c 1-5) ; end_index=$(tail -2 '//TRIM(gro_inputfile)//&
    &' | head -1 | cut -c 1-5) ; echo ${start_index}  ${end_index} >> s2e_'//TRIM(lio)//'.dat'  ! start 2 end 
    CALL execute_command_line(TRIM(Bashline)) 
    OPEN(UNIT=32,FILE='s2e_'//TRIM(lio)//'.dat',STATUS='unknown',IOSTAT=ierror,action='read')
    READ(32,*) start_index,end_index
    CLOSE(32)
    Bashline=' rm s2e_'//TRIM(lio)//'.dat '
    CALL execute_command_line(TRIM(Bashline)) 

    N_Resids=end_index-start_index+1
    WRITE(*,*)'  N_Resids= ',N_Resids,' in ',TRIM(gro_inputfile)
    WRITE(*,*)'  N_Atoms= ',Gesamtanzahl,' in ',TRIM(gro_inputfile)
ALLOCATE(atomsorteA(Gesamtanzahl))
ALLOCATE(koordA(Gesamtanzahl,3))
ALLOCATE(vel(Gesamtanzahl,3))           ! velocities(:,v_x,v_y,v_z)
ALLOCATE(R_u_I(N_Resids+1,2)) ! einen Eintrag mehr, damit an der letzten Stelle, der R_u_I: Resid_und_start_index(N_Resids,Gesamtanzahl+1,2) stehen kann, um die Schleifenbegrenzungen nicht anpassen zu muessen

ALLOCATE(residue_names(Gesamtanzahl))
ALLOCATE(residue_numbers(Gesamtanzahl))
ALLOCATE(atom_numbers(Gesamtanzahl))


!!! initialize system
call init_system(system)

koordA=0.0
vel=0.0
R_u_I=0
R_u_I(N_Resids+1,1)=-100
R_u_I(N_Resids+1,2)=Gesamtanzahl+1
residue_number_old=-1
i=0
N_DIPBI=0
N_PPDI=0
N_P3HT=0
N_P3HT_segments=0
N_PBDT_TS1=0
N_PBDT_TS1_segments=0
N_Species=0
DO j=1,Gesamtanzahl
    	READ(17, "(a)",IOSTAT=ierror) line
		BACKSPACE(UNIT=17)
        IF(ierror < 0) EXIT
        IF(ierror > 0) THEN
            WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(gro_inputfile),' aufgetreten.'//TRIM(line)
            WRITE(*,*) 'Beende Einlesen !'
            CALL EXIT (1)
        END IF
	!WRITE(*,*) 'laenge',LEN(TRIM(line)),' line: ',TRIM(line)
	IF (LEN(TRIM(line)) .le. 45 ) THEN ! check if no velocities are written
		!WRITE(*,*) 'case 1'
		!READ(17,*) dummyname
		!WRITE(*,*) dummyname
		!READ(17,'(i5,2a5,i5,3f8.3)') residue_number,residue_name,dummyname,i,x,y,z
		!WRITE(*,*) residue_number,residue_name,trim(dummyname),i,x,y,z
    		READ(17,'(i5,2a5,i5,3f8.3)',IOSTAT=ierror) residue_numbers(j),residue_names(j),atomsorteA(j),atom_numbers(j),&
                                            &koordA(j,1),koordA(j,2),koordA(j,3)
    		WRITE(*,'(i5,2a5,i5,3f8.3)',IOSTAT=ierror) residue_numbers(j),residue_names(j),atomsorteA(j),atom_numbers(j),&
                                           &koordA(j,1),koordA(j,2),koordA(j,3)

	ELSE ! ordinary *.gro file
    		Read(17,'(i5,2a5,i5,3f8.3,3f8.4)',IOSTAT=ierror) residue_numbers(j),residue_names(j),atomsorteA(j),atom_numbers(j),&
                                           	       &koordA(j,1),koordA(j,2),koordA(j,3),vel(j,1),vel(j,2),vel(j,3)
    		!WRITE(*,'(i5,2a5,i5,3f8.3,3f8.4)',IOSTAT=ierror) residue_numbers(j),residue_names(j),atomsorteA(j),atom_numbers(j),&
    		!                                      &koordA(j,1),koordA(j,2),koordA(j,3),vel(j,1),vel(j,2),vel(j,3)
	END IF
    IF(ierror < 0) EXIT
    IF(ierror > 0) THEN
        WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(gro_inputfile),' aufgetreten.'
        WRITE(*,*) 'Beende Einlesen !'
        CALL EXIT (1)
    END IF
    residue_number=residue_numbers(j)
    IF (residue_number /= residue_number_old ) THEN
        i=i+1
        R_u_I(i,1)=residue_number
        R_u_I(i,2)=j
        residue_number_old=residue_number
		IF((residue_names(j) == 'DIPBI') .OR. (residue_names(j) == 'DPBIK'))THEN
			N_DIPBI=N_DIPBI+1
		ELSE IF (residue_names(j)(1:4) == 'PPDI') THEN
			N_PPDI=N_PPDI+1
		ELSE IF (residue_names(j) == '8poly') THEN
			N_PBDT_TS1_segments=N_PBDT_TS1_segments+1
			IF( TRIM(adjustl(atomsorteA(j))) == 'CH3') THEN ! new PBDT_TS1 chain-Start with segment A and a CH3 group.
				N_PBDT_TS1=N_PBDT_TS1+1			
			END IF ! new PBDT_TS1 chain-Start
		ELSE IF ((residue_names(j)(1:3) == 'THP') .OR. (residue_names(j)(1:2) == 'PH') .OR. (residue_names(j)(1:2) == 'PM'))THEN
			N_P3HT_segments=N_P3HT_segments+1
			IF( (residue_names(j)(1:5) == 'THP1A') .OR. (residue_names(j)(5:5) == 'S') )THEN ! 32mer start
				N_P3HT=N_P3HT+1
			END IF ! 32mer start
		ELSE 
			N_Species=N_Species+1
		END IF
    END IF
    IF( i == N_Resids+1) EXIT
    IF (i==0 .and. j==Gesamtanzahl ) THEN
       write(*,*) ' Error: Es ist ein Fehler beim Einlesen ausgetreten: '//TRIM(gro_inputfile)
       write(*,*) ' ENDE '
       CALL EXIT (1)
   END IF  
END DO
Read(17,*,IOSTAT=ierror) x_box,y_box,z_box
IF( ierror > 0) STOP 'Error: Reading_the_boxsize'
CLOSE(17)


!!!! lambda data
IF(read_lambda_files)THEN
	CALL read_lambda_in_all_inputfile(N_types_lambda,lambda_array,lambda_in_all_inputfile)
	IF(DEBUG) THEN
		WRITE(*,*) 'Array with lambda data in eV!'
		DO i=1,N_types_lambda
			WRITE(*,*) lambda_array(i,:)
		END DO
	END IF
ELSE 
	IF (Allocated(lambda_array))  Deallocate(lambda_array)
	N_types_lambda=95
	ALLOCATE(lambda_array(N_types_lambda,7))
	WRITE(*,*) 'WARNING: Set lambda_data to zero!'
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Ende Einlesen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



R_offset=1-residue_numbers(1)
i_offset=1
WRITE(*,*) 'modify_mol_h_atoms:',modify_mol_h_atoms 
WRITE(*,*) 'Resid offset:', R_offset
WRITE(*,*) 'i_offset:',i_offset,' WARINING HARD CODED!'

WRITE(*,*) 'N_Resids: ',N_Resids
IF(N_DIPBI /=0) WRITE(*,*) 'N_DIPBI: ',N_DIPBI
IF(N_PPDI  /=0) WRITE(*,*) 'N_PPDI:  ',N_PPDI
IF(N_P3HT_segments /=0) WRITE(*,*) 'N_P3HT:    ',N_P3HT,' N_P3HT_segments:',N_P3HT_segments      !  Expected  N_P3HT=416   !416  ! Anzahl 32mere
IF(N_PBDT_TS1_segments /=0)  WRITE(*,*) 'N_PBDT_TS1:',N_PBDT_TS1,' N_PBDT_TS1_segments:',N_PBDT_TS1_segments   !N_PBDT_TS1= 464  !464= Anzahl ketten 
IF(N_Species /=0) WRITE(*,*)' N_Species:',N_Species


IF(nm_to_Ang) THEN
	koordA=koordA*10.
	x_box=x_box*10.
	y_box=y_box*10.
	z_box=z_box*10.
END IF



!!! Zaehle die Anzahl der Schwefelatome pro resid
IF ( Allocated(R_u_S))  Deallocate(R_u_S)
!! Residund Schwefel R_u_S(i,j)  i=1,N_Resids   j=1 => Resid ; j=2 N_Schwefel, j=3 Startsegment ==1 sonst 0 ; j==4 Endsegment sonst 0
ALLOCATE(R_u_S(N_Resids+1,4))
!! sequence list for polymers like PBDT-TS1
ALLOCATE(sequence_list(N_Resids))
ALLOCATE(residue_name_site(N_Resids))
residue_name_site=''
sequence_list=''
R_u_S=0
R_u_S(N_Resids+1,1)=N_Resids+1

IF(DIPBI_P3HT_calc) THEN
	DO i=1,N_Resids
		N_Schwefel=0
		DO j=R_u_I(i,2),R_u_I(i+1,2)
		       if (TRIM(adjustl(atomsorteA(j)))=='S') N_Schwefel=N_Schwefel+1
		       if (TRIM(adjustl(atomsorteA(j)))=='S-S_R') N_Schwefel=N_Schwefel+1
		END DO
		R_u_S(i,1)=R_u_I(i,1)
		R_u_S(i,2)=N_Schwefel
		!!! Startsegment 
		IF( (TRIM(adjustl(atomsorteA(R_u_I(i,2)))) == 'HC') .AND. (TRIM(adjustl(atomsorteA(R_u_I(i,2)+1))) == 'CA') ) THEN
			R_u_S(i,3)=1
		ELSE
			R_u_S(i,3)=0
		END IF

		!!! Endsegment
		IF( TRIM(adjustl(atomsorteA(R_u_I(i+1,2)-8))) == 'HC' ) THEN
			R_u_S(i,4)=1
		ELSE
			R_u_S(i,4)=0
		END IF
		!WRITE(*,*) R_u_S(i,1),R_u_S(i,2),R_u_S(i,3),R_u_S(i,4),atomsorteA(R_u_I(i+1,2)-8)
	END DO

ELSE IF (PPDI_PBDT_TS1_calc) THEN
     
	!! system2
	sequence=''
	DO i=1,N_Resids

		IF ( TRIM(residue_names(R_u_I(i,2))) == 'PPDI_' ) THEN 
			R_u_S(i,1)=R_u_I(i,1)			
			R_u_S(i,2)=1
			R_u_S(i,3)=0
			R_u_S(i,4)=0
			residue_name_site(i)='PPDI_'
			sequence_list(i)='PPDI'
			CYCLE
		END IF
		mol_start_da=.false.
		mol_ende_da=.false.
		N_Schwefel=0
		N_CarboxyS2=0
		N_PBDT=0
		! polymer
		DO j=R_u_I(i,2),R_u_I(i+1,2)
		       	!IF (TRIM(adjustl(atomsorteA(j)))=='S') N_Schwefel=N_Schwefel+1
		       	!IF (TRIM(adjustl(atomsorteA(j)))=='S-S_R') N_Schwefel=N_Schwefel+1
			IF ( residue_names(j) == '8poly' ) THEN 
	    			IF (TRIM(adjustl(atomsorteA(j)))=='F')  N_CarboxyS2=N_CarboxyS2+1
	    			IF (TRIM(adjustl(atomsorteA(j)))=='HC') N_PBDT=N_PBDT+1
				IF( (TRIM(adjustl(atomsorteA(j-1)))=='CH3') .and. (TRIM(adjustl(atomsorteA(j+3)))=='F') &
									   &.and. (TRIM(adjustl(atomsorteA(j+7)))=='OE')  )THEN  ! Start CarboxyS2 unit
					mol_start_da=.true.
					segtype='A'
					sequence=TRIM(sequence)//TRIM(segtype)
				ELSE IF( (TRIM(adjustl(atomsorteA(j)))=='CH3') .and. (TRIM(adjustl(atomsorteA(j-1)))=='HC') )THEN ! End chain 
					mol_ende_da=.true.
				ELSE IF( TRIM(adjustl(atomsorteA(j+3)))=='F' .and. TRIM(adjustl(atomsorteA(j-1)))=='HC') THEN ! Mittlel CarboxyS2
					segtype='A'
					sequence=TRIM(sequence)//TRIM(segtype)		
				ELSE IF (       (TRIM(adjustl(atomsorteA(j+1)))=='S') .and. (TRIM(adjustl(atomsorteA(j+10)))=='HC') ) THEN  ! PBDT segment
						IF  (   (TRIM(adjustl(atomsorteA(j+1))) =='S') .and. (TRIM(adjustl(atomsorteA(j))) =='C')  &    ! angle theta atoms
						& .and. (TRIM(adjustl(atomsorteA(j-2 ))) =='C') .and. (TRIM(adjustl(atomsorteA(j-1))) == 'S')) THEN 
							segtype='B'
							sequence=TRIM(sequence)//TRIM(segtype)
						END IF ! angle theta atoms
				END IF ! PBDT segment	
			END IF ! 8poly
		END DO !j
		R_u_S(i,1)=R_u_I(i,1)

		!!! Startsegment 
		IF(mol_start_da) THEN
			R_u_S(i,3)=1
		ELSE
			R_u_S(i,3)=0
		END IF

		IF( .NOT. mol_ende_da) THEN
				sequence=TRIM(sequence(:LEN_TRIM(sequence)-1)) 
				R_u_S(i,2)=PBDT_TS1_sequence_to_index(sequence) 
				sequence_list(i)=TRIM(sequence)
				sequence='' 
				R_u_S(i,4)=0
		ELSE !!! Startsegment 
				R_u_S(i,2)=PBDT_TS1_sequence_to_index(sequence) 
				sequence_list(i)=TRIM(sequence)  
				sequence=''
				R_u_S(i,4)=1	
		END IF 		

          	IF( mol_start_da .and. mol_ende_da ) THEN
				fin='A'
            	ELSE IF (mol_start_da) THEN
				fin='S'
		ELSE IF (mol_ende_da) THEN
				fin='E'
		ELSE
				fin='M'
		END IF

		IF( TRIM(sequence_list(i)(1:1)) == 'A' ) THEN
				residue_name_site(i)='A'//TRIM(str(N_CarboxyS2))//'B'//TRIM(str(N_PBDT/6))//TRIM(fin)
		ELSE IF( TRIM(sequence_list(i)(1:1)) == 'B' ) THEN
				residue_name_site(i)='B'//TRIM(str(N_PBDT/6))//'A'//TRIM(str(N_CarboxyS2))//TRIM(fin)
		END IF 

		 
		IF(DEBUG) WRITE(*,*) R_u_S(i,1),residue_name_site(i),R_u_S(i,2),R_u_S(i,3),R_u_S(i,4),TRIM(sequence_list(i))
	END DO !i=1,N_Resids
ELSE
	WRITE(*,*) 'Error: selected system is not supported.'
	CALL EXIT (1)
END IF



IF (modify_mol_h_atoms) THEN
	!! Check how much memory needs to be allocated => Result in GesamtanzahlH
	GesamtanzahlH=0
	!CALL how_many_h_atoms_to_mol(R_u_I,atomsorteA,residue_names,Gesamtanzahl,N_Resids,GesamtanzahlH,KETTE) 
	IF( DIPBI_P3HT_calc .OR. DIPBI_P3MT_calc) THEN
		IF(KETTE) THEN
		    GesamtanzahlH=GesamtanzahlH+N_DIPBI*218
		ELSE
		    GesamtanzahlH=GesamtanzahlH+N_DIPBI*86 ! Anzahl mit Methylgruppen
		END IF
	    
		IF(KETTE) THEN ! P3HT
			GesamtanzahlH=GesamtanzahlH+N_P3HT*802 ! Komplette Ketten
		ELSE IF (DIPBI_P3MT_calc) THEN
			GesamtanzahlH=GesamtanzahlH+N_P3HT*322  ! 32mer mit Methylgruppen
		END IF
	
		IF (add_virtual_H_atoms_to_gro) THEN ! Fuer jedes Bruchstuecke 2 HAtome mehr
			WRITE(*,*) '        GesamtanzahlH: ',GesamtanzahlH
			GesamtanzahlH=GesamtanzahlH+2*abs(N_Resids-N_DIPBI-N_P3HT)
			WRITE(*,*) 'virtual GesamtanzahlH: ',GesamtanzahlH
		END IF! no add_virtual_H_atoms_to_gro
	END IF ! DIPBI_P3HT_calc	

	IF (PPDI_PBDT_TS1_calc) THEN
		IF(KETTE) THEN !
		    GesamtanzahlH=GesamtanzahlH+N_PPDI*144
		ELSE IF(PPDI_Propyl) THEN
		    GesamtanzahlH=GesamtanzahlH+N_PPDI*96
		ELSE 
		    GesamtanzahlH=GesamtanzahlH+N_PPDI*78 ! Anzahl mit H-gruppen
		END IF !PPDI
		
		N_CarboxyS2=8
		N_PBDT=8
		GesamtanzahlH=GesamtanzahlH+N_PBDT_TS1*(38*N_PBDT+22*N_CarboxyS2+2)

		IF (add_virtual_H_atoms_to_gro) THEN ! Fuer jedes Bruchstuecke 2 HAtome mehr
			WRITE(*,*) '        GesamtanzahlH: ',GesamtanzahlH
			GesamtanzahlH=GesamtanzahlH+2*abs(N_Resids-N_PPDI-N_PBDT_TS1)
			WRITE(*,*) 'virtual GesamtanzahlH: ',GesamtanzahlH
		END IF ! no add_virtual_H_atoms_to_gro

	END IF ! PPDI_PBDT_TS1_calc
	GesamtanzahlH_precalc=GesamtanzahlH
ELSE
	GesamtanzahlH=Gesamtanzahl
	GesamtanzahlH_precalc=Gesamtanzahl
END IF  
 
!! Allocate Arrays
IF ( Allocated(R_u_I_H))  Deallocate(R_u_I_H)
IF ( Allocated(residue_numbersH))  Deallocate(residue_numbersH)
IF ( Allocated(residue_namesH))  Deallocate(residue_namesH)
IF ( Allocated(atomsorteH))  Deallocate(atomsorteH)
IF ( Allocated(koordH))  Deallocate(koordH)
IF ( Allocated(MolnameH))  Deallocate(MolnameH)
IF ( Allocated(elementsH))  Deallocate(elementsH)
IF ( Allocated(atom_numbersH))  Deallocate(atom_numbersH)
IF ( Allocated(velH))  Deallocate(velH)



ALLOCATE(R_u_I_H(N_Resids+1,2))
ALLOCATE(residue_numbersH(GesamtanzahlH))
ALLOCATE(residue_namesH(GesamtanzahlH))
ALLOCATE(atomsorteH(GesamtanzahlH))
ALLOCATE(koordH(GesamtanzahlH,3))
ALLOCATE(MolnameH(GesamtanzahlH))
ALLOCATE(elementsH(GesamtanzahlH))
ALLOCATE(atom_numbersH(GesamtanzahlH))
ALLOCATE(velH(GesamtanzahlH,3))


R_u_I_H=0
atomsorteH=''
koordH=0
GesamtanzahlH=0
MolnameH=''
residue_namesH=''
elementsH=''
atom_numbersH=0
velH=0
residue_numbersH=0

! Uebergabe der aller Resids
R_u_I_H(:,1)=R_u_I(:,1)


IF (modify_mol_h_atoms) THEN
    IF(	DIPBI_P3HT_calc .OR. DIPBI_P3MT_calc .OR. PPDI_PBDT_TS1_calc ) THEN
	IF (add_virtual_H_atoms_to_gro) THEN 
		    DO i=1,N_Resids
			ind1=R_u_I(i,1)+i_offset ! offset ggf anpassen !!!
			!WRITE(*,*) ind1,R_u_I(ind1,1),R_u_I(ind1,2)
			CALL add_h_atoms_to_mol(atomsorteA(R_u_I(ind1,2):(R_u_I(ind1+1,2)-1)),koordA(R_u_I(ind1,2):(R_u_I(ind1+1,2)-1),1:3),&
				    &(R_u_I(ind1+1,2)-R_u_I(ind1,2)),residue_names(R_u_I(ind1,2)),KETTE,atomsorte1H,koord1H,NAtoms1H,Molname)
			R_u_I_H(ind1,2)=GesamtanzahlH+1
			DO j=1,NAtoms1H            
				atomsorteH(GesamtanzahlH+j)=atomsorte1H(j)
				koordH(GesamtanzahlH+j,:)=koord1H(j,:)
				residue_namesH(GesamtanzahlH+j)=residue_names(ind1)
				MolnameH(GesamtanzahlH+j)=Molname
			END DO !j
			GesamtanzahlH=GesamtanzahlH+NAtoms1H
		    END DO ! N_Resids
		    R_u_I_H(N_Resids+1,2)=GesamtanzahlH+1
		    R_u_I_H(N_Resids+1,1)=-100
		    WRITE(*,*) 'GesamtanzahlH: ',GesamtanzahlH,' Gesamtanzahl: ',Gesamtanzahl
		    IF(GesamtanzahlH_precalc /= GesamtanzahlH) THEN	
			WRITE(*,*) 'WARNING: Number of estimatesd needed atoms including extra H atoms',GesamtanzahlH_precalc,&
			&',does not match the added number',GesamtanzahlH,'! Revise implementation or input data!'
			CALL EXIT (1)
		    END IF
	ELSE ! DO NOT add_virtual_H_atoms_to_gro => unite polymer chain segments to one entire chain and add H atoms afterwards
 
	  IF(DIPBI_P3HT_calc) THEN !!! DIPBI / PPDI section
		N_molecules_acceptor=N_DIPBI
		WRITE(*,*) 'Start to prepare DIPBI'
	  ELSE IF (PPDI_PBDT_TS1_calc) THEN	
	    	WRITE(*,*) 'Start to prepare PPDI'
		N_molecules_acceptor=N_PPDI
	  ELSE
		N_molecules_acceptor=0
	  END IF ! system 

	    !!! DIPBI / PPDI section
	    DO i=1,N_molecules_acceptor
		ind1=R_u_I(i,1)+i_offset ! offset ggf anpassen !!!
		!WRITE(*,*) ind1,R_u_I(ind1,1),R_u_I(ind1,2)
		CALL add_h_atoms_to_mol(atomsorteA(R_u_I(ind1,2):(R_u_I(ind1+1,2)-1)),koordA(R_u_I(ind1,2):(R_u_I(ind1+1,2)-1),1:3),&
		            &(R_u_I(ind1+1,2)-R_u_I(ind1,2)),residue_names(R_u_I(ind1,2)),KETTE,atomsorte1H,koord1H,NAtoms1H,Molname)
		R_u_I_H(ind1,2)=GesamtanzahlH+1		
		DO j=1,NAtoms1H            
		        atomsorteH(GesamtanzahlH+j)=atomsorte1H(j)
		        koordH(GesamtanzahlH+j,:)=koord1H(j,:)
			residue_namesH(GesamtanzahlH+j)=residue_names(ind1)
			MolnameH(GesamtanzahlH+j)=Molname
		END DO !j
		GesamtanzahlH=GesamtanzahlH+NAtoms1H
	    END DO ! N_Resids N_molecules_acceptor
	

	  	IF(print_coordH_to_xyz) THEN ! ok
			dummyname='xyz_Acceptor_H.xyz'
			INQUIRE(file=TRIM(dummyname),exist=Datei_vorhanden) 
			IF( .NOT. Datei_vorhanden) THEN
				OPEN(UNIT=8,FILE=TRIM(dummyname),STATUS='NEW',IOSTAT=ierror)
				WRITE(8,*) GesamtanzahlH
				WRITE(8,*)
				DO i=1,GesamtanzahlH
					WRITE(8,*) atomsorteH(i),koordH(i,1),koordH(i,2),koordH(i,3)
				END DO 
			ELSE
				WRITE(*,*) 'Warning: File already exists, so it will not be replaced:'//TRIM(dummyname)
			END IF 
		END IF ! print_coordH_to_xyz

		i_ende=R_u_I(N_molecules_acceptor+1,2)-1  ! offset ggf anpassen !!!
		WRITE(*,*) 'i_ende:',i_ende

	IF(DIPBI_P3HT_calc) THEN !!! unite P3HT or PBDT-TS1 units 
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		    !!!! Section unite all P3HT units to one entire molecule and add H atoms to 32mer  !!!!!!
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	 	    WRITE(*,*) 'Start to prepare P3HT'
		    DO i=1,N_P3HT !N_Resids
			i_start=i_ende+1
			i_ende=i_start+386-1
			!WRITE(*,*) 'tmp R and I',R_u_I(i,1),R_u_I(i,2)
			!ind1=R_u_I(i,1)+i_offset ! offset ggf anpassen !!!
			!WRITE(*,*) ind1,R_u_I(ind1,1),R_u_I(ind1,2)
			!i_max32_mer=R_u_I(ind1,2)+386 !+386+i_offset  !max index nach 32mer in gro file
			!WRITE(*,*) atomsorteA(R_u_I(1,2)),koordA(R_u_I(1,2),1:3)
			!WRITE(*,*) 'START: ',atomsorteA(i_start),koordA(i_start,1:3)
			!WRITE(*,*) 'ENDE: ',atomsorteA(i_ende),koordA(i_ende,1:3)
			!WRITE(*,*) ind1,R_u_I(ind1,1),R_u_I(ind1,2),ind1,i_max32_mer,residue_names(ind1),residue_names(ind1+1)

			CALL add_h_atoms_to_mol(atomsorteA(i_start:i_ende),koordA(i_start:i_ende,1:3),&
				    &(i_ende-i_start+1),residue_names(i_start),KETTE,atomsorte1H,koord1H,NAtoms1H,Molname)
			!!! HIER EINTEILUNG IN 32mere => wird unten korrigiert
			ind1=ind1+1
			R_u_I_H(ind1,2)=GesamtanzahlH+1
			DO j=1,NAtoms1H            
				atomsorteH(GesamtanzahlH+j)=atomsorte1H(j)
				koordH(GesamtanzahlH+j,:)=koord1H(j,:)
				residue_namesH(GesamtanzahlH+j)=residue_names(ind1)
				MolnameH(GesamtanzahlH+j)=Molname
				!WRITE(*,*)	atomsorte1H(j),koord1H(j,1),koord1H(j,2),koord1H(j,3)
			END DO !j
			GesamtanzahlH=GesamtanzahlH+NAtoms1H
		    END DO ! N_Resids


		    R_u_I_H(N_Resids+1,2)=GesamtanzahlH+1
		    R_u_I_H(N_Resids+1,1)=-100


		   ! WRITE(*,*) 'GesamtanzahlH: ',GesamtanzahlH,' Gesamtanzahl: ',Gesamtanzahl
		   IF(print_coordH_to_xyz) THEN !ok
				WRITE(*,*) ' print_coordH_to_xyz ' 
				dummyname='xyz_output_all_H.xyz'
				INQUIRE(file=TRIM(dummyname),exist=Datei_vorhanden) 
				IF( .NOT. Datei_vorhanden) THEN
					OPEN(UNIT=9,FILE=TRIM(dummyname),STATUS='REPLACE',IOSTAT=ierror)
					WRITE(9,*) GesamtanzahlH
					WRITE(9,*)
					DO i=1,GesamtanzahlH
						WRITE(9,*) atomsorteH(i),koordH(i,1),koordH(i,2),koordH(i,3)
					END DO 
					CLOSE(9)
				ELSE
					WRITE(*,*) 'Warning: File already exists, so it will not be replaced:'//TRIM(dummyname)
				END IF 
		   END IF
		   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		   !!!! Ende Section unite all P3HT !!!!!!!!!!!!!
		   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ELSE IF(PPDI_PBDT_TS1_calc) THEN 
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		    !!!! Section unite all PBDT_TS1 units to one entire molecule and add H atoms to A8B8A chain ABABABABABABABAB  !!!!!!
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	 	    WRITE(*,*) 'Start to prepare P3HT'
		    DO i=1,N_PBDT_TS1  !N_Resids
			i_start=i_ende+1
			i_ende=i_start+378-1
			!WRITE(*,*) 'tmp R and I',R_u_I(i,1),R_u_I(i,2)
			!ind1=R_u_I(i,1)+i_offset ! offset ggf anpassen !!!
			!WRITE(*,*) ind1,R_u_I(ind1,1),R_u_I(ind1,2)
			!i_max32_mer=R_u_I(ind1,2)+378 !+378+i_offset  !max index nach 32mer in gro file
			!WRITE(*,*) atomsorteA(R_u_I(1,2)),koordA(R_u_I(1,2),1:3)
			!WRITE(*,*) 'START: ',atomsorteA(i_start),koordA(i_start,1:3)
			!WRITE(*,*) 'ENDE: ',atomsorteA(i_ende),koordA(i_ende,1:3)
			!WRITE(*,*) ind1,R_u_I(ind1,1),R_u_I(ind1,2),ind1,i_max32_mer,residue_names(ind1),residue_names(ind1+1)

			CALL add_h_atoms_to_mol(atomsorteA(i_start:i_ende),koordA(i_start:i_ende,1:3),&
				    &(i_ende-i_start+1),residue_names(i_start),KETTE,atomsorte1H,koord1H,NAtoms1H,Molname)
			!!! HIER EINTEILUNG IN A8B8A 
			ind1=ind1+1
			R_u_I_H(ind1,2)=GesamtanzahlH+1
			DO j=1,NAtoms1H            
				atomsorteH(GesamtanzahlH+j)=atomsorte1H(j)
				koordH(GesamtanzahlH+j,:)=koord1H(j,:)
				residue_namesH(GesamtanzahlH+j)=residue_names(ind1)
				MolnameH(GesamtanzahlH+j)=Molname
				!WRITE(*,*)	atomsorte1H(j),koord1H(j,1),koord1H(j,2),koord1H(j,3)
			END DO !j
			GesamtanzahlH=GesamtanzahlH+NAtoms1H
		    END DO ! N_Resids


		    R_u_I_H(N_Resids+1,2)=GesamtanzahlH+1
		    R_u_I_H(N_Resids+1,1)=-100


		   ! WRITE(*,*) 'GesamtanzahlH: ',GesamtanzahlH,' Gesamtanzahl: ',Gesamtanzahl
		   IF(print_coordH_to_xyz) THEN !ok
				WRITE(*,*) ' print_coordH_to_xyz ' 
				dummyname='xyz_output_all_H.xyz'
				INQUIRE(file=TRIM(dummyname),exist=Datei_vorhanden) 
				IF( .NOT. Datei_vorhanden) THEN
					OPEN(UNIT=90,FILE=TRIM(dummyname),STATUS='REPLACE',IOSTAT=ierror)
					WRITE(90,*) GesamtanzahlH
					WRITE(90,*)
					DO i=1,GesamtanzahlH
						WRITE(90,*) atomsorteH(i),koordH(i,1),koordH(i,2),koordH(i,3)
					END DO 
					CLOSE(90)
				ELSE
					WRITE(*,*) 'Warning: File already exists, so it will not be replaced:'//TRIM(dummyname)
				END IF 
		   END IF
	   	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   	    !!!! Ende Section unite all PBDT_TS1 !!!!!!!!!!!!!
	   	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	   END IF !!! Unite polymer segments to an entire chain.
	END IF !! add virtual H atoms to gro ??!
		
    END IF ! DIPBI_P3HT_calc .OR. DIPBI_P3MT_calc .OR. PPDI_PBDT_TS1_calc

    IF(GesamtanzahlH_precalc /= GesamtanzahlH) THEN	
	WRITE(*,*) 'WARNING: Number of estimated needed atoms including extra H atoms',GesamtanzahlH_precalc,&
	&',does not match the added number',GesamtanzahlH,'! Revise implementation or input data!'
	CALL EXIT (1)
    END IF
 ELSE ! modify_mol_h_atoms
    ! Uebergabe ohne H-Veränderung
    Do i=1,Gesamtanzahl
        atomsorteH(i)=atomsorteA(i)
        koordH(i,:)=koordA(i,:)
    END DO
    DO i=1,N_Resids
        R_u_I_H(i,:)=R_u_I(i,:)
        residue_namesH(i)=residue_names(i)
	MolnameH(i)=residue_names(i)
    END DO
    GesamtanzahlH=Gesamtanzahl
    R_u_I_H(N_Resids+1,2)=GesamtanzahlH+1
    R_u_I_H(N_Resids+1,1)=-100
 END IF ! modify_mol_h_atoms




!!! Ubergabe aller elemente an das elementsH array, da die atomsorteH unten modifiziert wird
elementsH=atomsorteH

!!!!!!!!  Section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Residue name list for all possible names !!!!!!!!!!!!!!!!!!!!!!!!!!

IF(DIPBI_P3HT_calc .OR. DIPBI_P3MT_calc) THEN ! DIPBI_P3HT: 95
	N_seg_types_all=95
ELSE IF (PPDI_PBDT_TS1_calc) THEN  ! PPDI_PBDT_TS1: 92
	N_seg_types_all=92
END IF


IF ( Allocated(residue_name_list))  Deallocate(residue_name_list)
IF ( Allocated(molecule_da))  Deallocate(molecule_da)
IF ( Allocated(iHOMO_list))  Deallocate(iHOMO_list)
IF ( Allocated(lambda_array_all))  Deallocate(lambda_array_all)
ALLOCATE(residue_name_list(N_seg_types_all))
ALLOCATE(molecule_da(N_seg_types_all))
ALLOCATE(iHOMO_list(N_seg_types_all))
ALLOCATE(lambda_array_all(N_seg_types_all,7))

residue_name_list=''
molecule_da=.false.
iHOMO_list=0
lambda_array_all=0


IF (PPDI_PBDT_TS1_calc) THEN  !!! Start Residue name list PPDI_PBDT_TS1
	!!! PPDI
	k=1
	residue_name_list(k)='PPDI_'
	iHOMO_list(k)=205
	lambda_array_all(k,:)=lambda_array(1,:)

	!!! PBDT_TS1
	DO j=1,3  !!! Startsegment, Mittelsegment, Ensegment
		IF (j==1) fin='M'	
		IF (j==2) fin='S'
		IF (j==3) fin='E'	
		
		!!! chains starting with segment A
		sequence=''
		iHOMO=0
		DO i=1,8
			k=k+1
			segtype='A'
			sequence=TRIM(sequence)//TRIM(segtype)	
			residue_name_list(k)='A'//TRIM(str(i))//'B'//TRIM(str(i-1))//TRIM(fin)
			IF(i==1) THEN
				iHOMO=40
			ELSE 
				iHOMO=iHOMO+39
			END IF		
			iHOMO_list(k)=iHOMO
			lambda_array_all(k,:)=lambda_array(PBDT_TS1_sequence_to_index(sequence),:) 

			IF (i==8) CYCLE ! exclude case for complete chain
			k=k+1
			segtype='B'
			sequence=TRIM(sequence)//TRIM(segtype)	
			residue_name_list(k)='A'//TRIM(str(i))//'B'//TRIM(str(i))//TRIM(fin)
			iHOMO=iHOMO+64
			iHOMO_list(k)=iHOMO
			lambda_array_all(k,:)=lambda_array(PBDT_TS1_sequence_to_index(sequence),:)
		END DO !i segment A
		
		!!! chains starting with segment B
		sequence=''
		iHOMO=0
		DO i=1,8
			k=k+1
			segtype='B'
			sequence=TRIM(sequence)//TRIM(segtype)	
			residue_name_list(k)='B'//TRIM(str(i))//'A'//TRIM(str(i-1))//TRIM(fin)
			IF(i==1)THEN
				iHOMO=65
			ELSE
				iHOMO=iHOMO+64
			END IF
			iHOMO_list(k)=iHOMO
			lambda_array_all(k,:)=lambda_array(PBDT_TS1_sequence_to_index(sequence),:) 

			IF (i==8) CYCLE ! exclude case for complete chain
			k=k+1
			segtype='A'
			sequence=TRIM(sequence)//TRIM(segtype)	
			residue_name_list(k)='B'//TRIM(str(i))//'A'//TRIM(str(i))//TRIM(fin)
			iHOMO=iHOMO+39
			iHOMO_list(k)=iHOMO
			lambda_array_all(k,:)=lambda_array(PBDT_TS1_sequence_to_index(sequence),:)
		END DO !i segment B
	END DO !j fin= S,M,E
	
	!!! Complete chain
	k=k+1
	i=8
	segtype='A'
	sequence=TRIM(segtype)//TRIM(sequence)	
	residue_name_list(k)='A'//TRIM(str(i))//'B'//TRIM(str(i))//'A'
	iHOMO_list(k)=825
	lambda_array_all(k,:)=lambda_array(PBDT_TS1_sequence_to_index(sequence),:) 
	IF(DEBUG) THEN
		DO k=1,N_seg_types_all
			WRITE(*,*) residue_name_list(k),iHOMO_list(k),lambda_array_all(k,:)
		END DO
	END IF 

      !!! Ende Residue name list PPDI_PBDT_TS1
      
ELSE  !!! START Residue name list for all possible names in DIPBI P3HT/P3MT
	k=0
	DO i=1,32  ! MITTE
		k=k+1
		IF(DIPBI_P3HT_calc) THEN
	 		residue_namesH_dummy='PH'
			iHOMO_list(k)=30*i+1
		ELSE IF (DIPBI_P3MT_calc) THEN
			residue_namesH_dummy='PM'
			iHOMO_list(k)=15*i+1
		END IF
		IF (i .LE. 9) THEN
			residue_namesH_dummy=TRIM(ADJUSTL(residue_namesH_dummy))//TRIM(str(0))//TRIM(str(i))//'M'
		ELSE
			residue_namesH_dummy=TRIM(ADJUSTL(residue_namesH_dummy))//TRIM(str(i))//'M'
		END IF
		lambda_array_all(k,:)=lambda_array(i,:)
		residue_name_list(k)=residue_namesH_dummy
	END DO

	DO i=1,31  ! START
		k=k+1
		IF(DIPBI_P3HT_calc) THEN
	 		residue_namesH_dummy='PH'
			iHOMO_list(k)=30*i+1
		ELSE IF (DIPBI_P3MT_calc) THEN
			residue_namesH_dummy='PM'
			iHOMO_list(k)=15*i+1
		END IF
		IF (i .LE. 9) THEN
			residue_namesH_dummy=TRIM(ADJUSTL(residue_namesH_dummy))//TRIM(str(0))//TRIM(str(i))//'S'
		ELSE
			residue_namesH_dummy=TRIM(ADJUSTL(residue_namesH_dummy))//TRIM(str(i))//'S'
		END IF
		lambda_array_all(k,:)=lambda_array(i,:)
		residue_name_list(k)=residue_namesH_dummy
	END DO

	DO i=1,31  ! ENDE
		k=k+1
		IF(DIPBI_P3HT_calc) THEN
	 		residue_namesH_dummy='PH'
			iHOMO_list(k)=30*i+1
		ELSE IF (DIPBI_P3MT_calc) THEN
			residue_namesH_dummy='PM'
			iHOMO_list(k)=15*i+1
		END IF
		IF (i .LE. 9) THEN
			residue_namesH_dummy=TRIM(ADJUSTL(residue_namesH_dummy))//TRIM(str(0))//TRIM(str(i))//'E'
		ELSE
			residue_namesH_dummy=TRIM(ADJUSTL(residue_namesH_dummy))//TRIM(str(i))//'E'
		END IF
		lambda_array_all(k,:)=lambda_array(i,:)
		residue_name_list(k)=residue_namesH_dummy
	END DO

	IF(DIPBI_P3HT_calc) THEN!!! DIPBI an residue_name_list(95) mit KETTE DPBIK
		residue_name_list(32)='P3HTA'
		residue_name_list(95)='DPBIK'
		lambda_array_all(95,:)=lambda_array(33,:)
		iHOMO_list(95)=293 !PM3 DIPBI Hexyl chains
	ELSE IF (DIPBI_P3MT_calc) THEN
		residue_name_list(32)='P3MTA'
		residue_name_list(95)='DIPBI'
		lambda_array_all(95,:)=lambda_array(33,:)
		iHOMO_list(95)=161 !PM3 DIPBI Methy groups
	END IF
END IF!!! END  Residue name list for all possible names in DIPBI P3HT/P3MT





IF (DEBUG) THEN
	WRITE(*,*) 'List with possible residue names'
	DO i=1,N_seg_types_all
		WRITE(*,*) i,TRIM(residue_name_list(i)),'   iHOMO:',iHOMO_list(i)
	END DO 
END IF


!IF (add_virtual_H_atoms_to_gro) THEN 
!	DO i=1,N_Resids
!		!!! Uebergabe der Resid und Nummerrierung der Elementenamen
!		k=0
!		DO j=R_u_I_H(i,2),R_u_I_H(i+1,2)-1
!			residue_namesH(j)=residue_namesH_dummy
!			k=k+1
!			atomsorteH(j)=TRIM(ADJUSTL(atomsorteH(j)))//TRIM(str(k))
!			residue_numbersH(j)=R_u_I_H(i,1)
!		END DO
!	END DO ! i=1,N_Resids
!ELSE
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!         START  SECTION         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!! Modify new residue_namesH      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF(DIPBI_P3HT_calc) THEN 
		ind1=R_u_I_H(N_DIPBI+1,2) ! Startindex des ersten Molekuels
		WRITE(*,*) 'START P3HT ind:',ind1
		DO i=1,N_Resids
			IF ( .NOT. add_virtual_H_atoms_to_gro) THEN 
				!!! Modify Resid index R_u_I mit einzelsplitting 
				IF(DIPBI_P3HT_calc .AND. (i .GT. N_DIPBI) ) THEN
					!Startindex fuer das nachste Element ind1 + 25* N_Schwefel + H_Start + H_Ende	
					!WRITE(*,*) R_u_S(i,1),R_u_S(i,2),R_u_S(i,3),R_u_S(i,4)			
					R_u_I_H(i+1,2)=ind1+25*R_u_S(i,2)+R_u_S(i,3)+R_u_S(i,4)
					ind1=R_u_I_H(i+1,2)
				END IF
			END IF !add_virtual_H_atoms_to_gro

			!!! Modify residue_namesH
			IF (R_u_S(i,2) == 0) THEN ! DIPBI
				residue_namesH_dummy=residue_name_list(95)
				IF( .NOT. molecule_da(95)) molecule_da(95)=.true.
			ELSE IF(R_u_S(i,3)==1 .AND. R_u_S(i,4)==1) THEN  !32mer
				residue_namesH_dummy=residue_name_list(32)
				IF( .NOT. molecule_da(32)) molecule_da(32)=.true.
			ELSE
				!!! Schwefel
				IF(R_u_S(i,3)==1) THEN !Start
					residue_namesH_dummy=residue_name_list(R_u_S(i,2)+32)
					IF( .NOT. molecule_da(R_u_S(i,2)+32)) molecule_da(R_u_S(i,2)+32)=.true.
				ELSE IF (R_u_S(i,4)==1) THEN !Ende
					residue_namesH_dummy=residue_name_list(R_u_S(i,2)+63)
					IF( .NOT. molecule_da(R_u_S(i,2)+63)) molecule_da(R_u_S(i,2)+63)=.true.
				ELSE ! Hauptteil Mittelstueck
				 	residue_namesH_dummy=residue_name_list(R_u_S(i,2))
					IF( .NOT. molecule_da(R_u_S(i,2))) molecule_da(R_u_S(i,2))=.true.
				END IF		
			END IF
	
			!!! Uebergabe der Resid und Nummerrierung der Elementenamen
			k=0
			DO j=R_u_I_H(i,2),R_u_I_H(i+1,2)-1
				residue_namesH(j)=residue_namesH_dummy
				k=k+1
				atomsorteH(j)=TRIM(ADJUSTL(atomsorteH(j)))//TRIM(str(k))
				residue_numbersH(j)=R_u_I_H(i,1)
			END DO
			IF(i == N_Resids) THEN
				WRITE(*,*)  'Schwefel: ',R_u_S(i,1),R_u_S(i,2),R_u_S(i,3),R_u_S(i,4)	
				WRITE(*,*)  'last Resid:   ',i,R_u_I_H(i,1),residue_numbersH(R_u_I_H(i,1)),TRIM(residue_namesH_dummy)
			END IF
		END DO
	ELSE IF	(PPDI_PBDT_TS1_calc) THEN
		ind1=R_u_I_H(N_PPDI+1,2) ! Startindex des ersten Molekuels
		WRITE(*,*) 'START PBDT_TS1 ind:',ind1
		DO i=1,N_Resids
				IF ( .NOT. add_virtual_H_atoms_to_gro) THEN 
					!!! Modify Resid index R_u_I mit einzelsplitting 
					IF(PPDI_PBDT_TS1_calc .AND. (i .GT. N_PPDI) ) THEN
						!!!   22*N_CarboxyS2 +38 N_PBDT + H_Start + H_Ende
						k=22*countsubstring(sequence_list(i), 'A')+38*countsubstring(sequence_list(i), 'B')	
						R_u_I_H(i+1,2)=ind1+k+R_u_S(i,3)+R_u_S(i,4)
						ind1=R_u_I_H(i+1,2)
					END IF
				END IF !add_virtual_H_atoms_to_gro

				!!! Modify residue_namesH
				IF (R_u_S(i,2) == 1) THEN ! PPDI
					residue_namesH_dummy=residue_name_list(1)
					IF( .NOT. molecule_da(1)) molecule_da(1)=.true.
				ELSE IF(R_u_S(i,3)==1 .AND. R_u_S(i,4)==1) THEN  ! PBDT_TS1 complete chain A8B8A
					residue_namesH_dummy=residue_name_list(92)
					IF( .NOT. molecule_da(92)) molecule_da(92)=.true.
				ELSE
				
				!!! Schwefel
				IF(R_u_S(i,3)==1) THEN !Start
					k_ind=seq_to_k_ind(sequence_list(i),'S')
					residue_namesH_dummy=residue_name_list(k_ind)
					IF( .NOT. molecule_da(k_ind)) molecule_da(k_ind)=.true.
					!WRITE(*,*) 'TMP:',"  ",molecule_da(k_ind),k_ind,residue_namesH_dummy,'  ',TRIM(sequence_list(i))
					!CALL EXIT
				ELSE IF (R_u_S(i,4)==1) THEN !Ende
					k_ind=seq_to_k_ind(sequence_list(i),'E')
					residue_namesH_dummy=residue_name_list(k_ind)
					IF( .NOT. molecule_da(k_ind)) molecule_da(k_ind)=.true.
				ELSE ! Hauptteil Mittelstueck
					k_ind=seq_to_k_ind(sequence_list(i),'M')
				 	residue_namesH_dummy=residue_name_list(k_ind)
					IF( .NOT. molecule_da(k_ind)) molecule_da(k_ind)=.true.
				END IF		
			END IF

			IF(DEBUG) WRITE(*,*) i,residue_namesH_dummy
			!!! Uebergabe der Resid und Nummerrierung der Elementenamen
			k=0
			DO j=R_u_I_H(i,2),R_u_I_H(i+1,2)-1
				residue_namesH(j)=residue_namesH_dummy
				k=k+1
				atomsorteH(j)=TRIM(ADJUSTL(atomsorteH(j)))//TRIM(str(k))
				residue_numbersH(j)=R_u_I_H(i,1)
				
			END DO !j
			IF(i == N_Resids) THEN
				WRITE(*,*)  'Schwefel:',R_u_S(i,1),R_u_S(i,2),R_u_S(i,3),R_u_S(i,4)	
				WRITE(*,*)  'Resid: ',i,R_u_I_H(i,1),residue_numbersH(R_u_I_H(i,1)),TRIM(residue_namesH_dummy)
				!CALL EXIT (1)
			END IF
		END DO ! i=1,N_Resids

	END IF ! DIPBI_P3HT_calc or  PPDI_PBDT_TS1_calc

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
	!!!  		Modify new residue_namesH   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!    END SECTION   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!							
!END IF !add_virtual_H_atoms_to_gro



IF(DEBUG) THEN
		WRITE(*,*) GesamtanzahlH
		WRITE(*,*) ' '
		DO j=1,GesamtanzahlH
			WRITE(*,*) atomsorteH(j),koordH(j,1),koordH(j,2),koordH(j,3)
		END DO 
END IF




!!! SETUP data in system_d
IF( use_min_d_COM_vector_list) THEN
	WRITE(*,*) 'WARNING: Implementation not finished!  Set use_min_d_COM_vector_list=.false.'
	call make_min_d_COM_vector_list(system,N_Resids,R_u_I_H,GesamtanzahlH,residue_numbersH,residue_namesH,elementsH,koordH,&
		&neighbourlist_filename,gro_inputfile,N_DIPBI,N_P3HT,N_PPDI,N_PBDT_TS1,N_Species,N_Neighbours,neighbour_list_Res,&
		&x_box,y_box,z_box)
ENDIF




IF(use_I_shift) THEN !!! Shift ggf the index as VOTCA does not allow to start with resid 0
	R_u_I(:,1)=R_u_I(:,1)+I_shift
	R_u_I_H(:,1)=R_u_I_H(:,1)+I_shift
	residue_numbersH(:)=residue_numbersH(:)+I_shift
	WRITE(*,*) 'Use I_shift: ',I_shift
ELSE
	IF( R_u_I(1,1) == 0) THEN
		WRITE(*,*) 'WARNING: Votca needs to start with resid 1, apply I_shift=1!!'
	END IF
END IF







IF(create_fake_topology) THEN
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! START create fake topology files !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	OPEN(UNIT=20,FILE=TRIM(top_all_outputfile),STATUS='REPLACE',IOSTAT=ierror)

	WRITE(20,*) '[ defaults ]                                                    '
	WRITE(20,*) '; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeLJ '
	WRITE(20,*) '       1                3             yes           0.5     0.5 '
	WRITE(20,*) '                                                                '
	WRITE(20,*) '[ atomtypes ]                                                   '
	WRITE(20,*) '; name       mass      charge    ptype        sigma          eps'
	!!! Element section
	WRITE(20,*) '     C     12.000       0.000        A  0.00000e+00  0.00000e+00'
	WRITE(20,*) '     N     14.000       0.000        A  0.00000e+00  0.00000e+00'
	WRITE(20,*) '     O     16.000       0.000        A  0.00000e+00  0.00000e+00'
	WRITE(20,*) '     H      1.000       0.000        A  0.00000e+00  0.00000e+00'
	WRITE(20,*) '     S     32.000       0.000        A  0.00000e+00  0.00000e+00'
	IF ( N_DIPBI .GT. 0 ) THEN
		WRITE(20,*) '    Cl     35.000       0.000        A  0.00000e+00  0.00000e+00'
	END IF
	IF ( (N_PBDT_TS1 .GT. 0 ) .AND. PPDI_PBDT_TS1_calc) THEN
		WRITE(20,*) '     F     19.000       0.000        A  0.00000e+00  0.00000e+00'
	END IF

	WRITE(20,*) '                                                                '
	DO i=1,size(residue_name_list)
		IF(molecule_da(i)) THEN
			WRITE(20,*) '#include "'//TRIM(residue_name_list(i))//'_raw.itp"'
		END IF
	END DO 
	WRITE(20,*) '                                                                '
	WRITE(20,*) '[ system ]                                                      '
	WRITE(20,*) '; Name                                                          '
	WRITE(20,*) 'Protein                                                         '
	WRITE(20,*) '                                                                '
	WRITE(20,*) '[ molecules ]                                                   '
	WRITE(20,*) '; Compound        #mols                                         '

	IF(DIPBI_P3HT_calc .OR. DIPBI_P3MT_calc) THEN
		IF (N_DIPBI .GT. 0 ) THEN
			WRITE(20,*) TRIM(residue_name_list(95))//'   ',N_DIPBI
		END IF
		DO i=N_DIPBI+1,N_Resids
			WRITE(20,*) TRIM(residue_namesH(R_u_I_H(i,2)))//'   ',1
		END DO
	ELSE IF (PPDI_PBDT_TS1_calc) THEN
		IF (N_PPDI .GT. 0 ) THEN
			WRITE(20,*) TRIM(residue_name_list(1))//'   ',N_PPDI
		END IF
		DO i=N_PPDI+1,N_Resids
			WRITE(20,*) TRIM(residue_namesH(R_u_I_H(i,2)))//'   ',1
		END DO
	END IF
	WRITE(20,*) '                                                                '
	CLOSE(20)
END IF ! create_fake_topology



!! Create VOTCA mapping file
IF(create_mapping_file) THEN
	OPEN(UNIT=22,FILE=TRIM(map_filename),STATUS='replace', access = 'append')
	WRITE(22,'(g0)') '<topology>'
	WRITE(22,'(g0)') '	<molecules>'
	!WRITE(22,'(g0)') '		<molecule>'
	CLOSE(22)
END IF !create_mapping_file



IF(create_fake_topology .OR. create_mapping_file .OR. create_geo_orbs_mps_g09) THEN
	!!! Make fake *.itp Dateien erzeugen
	DO j=1,size(molecule_da)
		IF (molecule_da(j)) THEN
			DO i=1,N_Resids
				IF( residue_name_list(j) == residue_namesH(R_u_I_H(i,2)) ) THEN
					IF(create_fake_topology) THEN
			 			CALL make_fake_itp(R_u_I_H(i+1,2)-R_u_I_H(i,2),residue_name_list(j),&
						&elementsH(R_u_I_H(i,2):R_u_I_H(i+1,2)-1),atomsorteH(R_u_I_H(i,2):R_u_I_H(i+1,2)-1))
						WRITE(*,*) 'ADDED: ',residue_name_list(j)//'_raw.itp'
					END IF !create_fake_topology				

					IF(create_mapping_file) THEN
						CALL make_votca_kmc_mapfile(map_filename,R_u_I_H(i+1,2)-R_u_I_H(i,2),&
						&residue_namesH(R_u_I_H(i,2):R_u_I_H(i+1,2)-1),&
						&elementsH(R_u_I_H(i,2):R_u_I_H(i+1,2)-1),atomsorteH(R_u_I_H(i,2):R_u_I_H(i+1,2)-1),&
						&iHOMO_list(j),lambda_array_all(j,:),orbs_ZINDO,add_virtual_H_atoms_to_gro,&
						&exclude_virtual_H_atoms)
					END IF ! create_mapping_file
	
					IF(create_geo_orbs_mps_g09) THEN ! make inputfiles if selected in extra_files_modus
						CALL  make_geo_orbs_mps_g09(extra_files_modus,R_u_I_H(i+1,2)-R_u_I_H(i,2),&
							& residue_namesH(R_u_I_H(i,2)),elementsH(R_u_I_H(i,2):R_u_I_H(i+1,2)-1),&
							& koordH(R_u_I_H(i,2):R_u_I_H(i+1,2)-1,1:3),&
							& trajstep,Numcharged,NCPUS,Mem)
					END IF ! create_geo_orbs_mps_g09

					EXIT 
				END IF
			END DO
			!WRITE(*,*) 'Checked Molecule: ',TRIM(residue_name_list(j))
		END IF 
	END DO
END IF !create_fake_topology OR:1 create_mapping_file
!!!HIER WEITER!!! 


IF(create_mapping_file) THEN
	OPEN(UNIT=22,FILE=TRIM(map_filename),STATUS='unknown', access = 'append')
	!WRITE(22,'(g0)') '		</molecule>'
	WRITE(22,'(g0)') '	</molecules>'
	WRITE(22,'(g0)') '</topology>'
	CLOSE(22)
	WRITE(*,*) 'Votca mapping file generated: ',TRIM(map_filename)
END IF !create_mapping_file


!!!! ENDE create fake topology file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(create_gromp_mdp_file) THEN
	CALL make_gromp_mdp_file()
END IF

IF(create_mapping_file) THEN
	! Creates a file for the constrain requirements, needs to be included into the VOTCA options.xml file 
	! to build the neighbourlist in VOTCA
	CALL make_votca_neighbourlist_constrains(N_seg_types_all,residue_name_list,molecule_da,neighbourlist_constrains_filename)
END IF



!!! Create molpol options files for:    ctp_tools -e molpol -o options_molpol.xml !!!
IF(create_molpol_optionsfiles) THEN
	DO i=1,size(molecule_da)
		IF(molecule_da(i)) THEN
			residue_name=TRIM(residue_name_list(i))			
			IF ((residue_name(:2) == 'PM') .OR. (residue_name(:2) == 'PH')) THEN ! P3MT oder P3HT to All "A"
				residue_name=TRIM(residue_name(1:4))//'A'
			ELSE IF((residue_name(1:1) == 'A') .OR. (residue_name(1:1) == 'B')) THEN  ! PBDT-TS1
				residue_name=TRIM(residue_name(1:4))//'A'
			END IF
			!WRITE(*,*) TRIM(residue_name)
			CALL make_votca_molpol_optionsfiles(residue_name)
		END IF
	END DO
	WRITE(*,*) 'Use:  ctp_tools -e molpol -o options_molpol_RESNAME_C.xml '
	WRITE(*,*) 'make_votca_molpol_optionsfiles: done'
END IF



!!! Convert back to nm
IF(nm_to_Ang) THEN ! Zurueck zu nm
	koordA=koordA*0.1
	koordH=koordH*0.1
	x_box=x_box*0.1
	y_box=y_box*0.1
	z_box=z_box*0.1
	IF(create_votca_boxfile) THEN ! in nm!
		CALL make_votca_boxfile(x_box,y_box,z_box,gro_inputfile)
	END IF
END IF





WRITE(*,*) 'N_Resids: ',N_Resids
!! Output GROMACS DATEI *.gro
OPEN(UNIT=18,FILE=TRIM(gro_outputfile),STATUS='REPLACE',IOSTAT=ierror)
WRITE(18,*) 'enumerated_gro_to_VOTCA__'//TRIM(gro_inputfile)
WRITE(18,*) GesamtanzahlH 
DO j=1,GesamtanzahlH
	atom_numbersH(j)=mod(j,100000)
    WRITE(18,'(i5,2a5,i5,3f8.3,3f8.4)',IOSTAT=ierror) residue_numbersH(j),residue_namesH(j),adjustr(atomsorteH(j)),&
					   &atom_numbersH(j),&
                                           &koordH(j,1),koordH(j,2),koordH(j,3),velH(j,1),velH(j,2),velH(j,3)
    IF(ierror > 0) THEN
        WRITE(*,*) 'Error: Es ist eine Fehler beim Schreiben von ',TRIM(gro_outputfile),' aufgetreten.'
    END IF
END DO
WRITE(18,'(3f8.3)',IOSTAT=ierror) x_box,y_box,z_box
IF(ierror > 0) THEN
        WRITE(*,*) 'Error: Es ist eine Fehler beim Schreiben von ',TRIM(gro_outputfile),' aufgetreten.'
END IF
WRITE(*,*) '#New gro file created: '//TRIM(gro_outputfile)
WRITE(*,*) ' '



WRITE(*,*) '###Next Steps TO DO: '
WRITE(*,*) '#!/bin/bash'
WRITE(*,*) ' '
WRITE(*,*) ' ### 1) create GROMACS *.tpr'
WRITE(*,*) ' /opt/gromacs-5.0.4-d-WP/bin/grompp_d  -c '//TRIM(gro_outputfile)//' -p '//TRIM(top_all_outputfile)//' -f gromp.mdp'
WRITE(*,*) ' ### 2) GROMACS *.tpr'
WRITE(*,*) ' /opt/gromacs-5.0.4-d-WP/bin/mdrun_d  -s topol.tpr -x '//TRIM(gro_inputfile(:LEN_TRIM(gro_inputfile)-4))//'.xtc'
WRITE(*,*) ' '
WRITE(*,*) ' ### 3a) Mapping'
WRITE(*,*) ' ctp_map -t topol.tpr -c '//TRIM(gro_inputfile(:LEN_TRIM(gro_inputfile)-4))//'.xtc -s '//&
									&TRIM(map_filename)//' -f state.sql'
WRITE(*,*) ' ### 3b) System.xml'
WRITE(*,*) ' cp '//TRIM(map_filename)//' system.xml'
WRITE(*,*) ' ### 4) Neighbor list'
WRITE(*,*) ' ctp_run -e neighborlist -o '//TRIM(neighbourlist_constrains_filename)//' -f state.sql'
WRITE(*,*) ' '

WRITE(*,*) ' ### 5)  Site energies'
WRITE(*,*) ' ctp_run -e emultipole -o  '//TRIM(votca_xml_options_filename)//' -f state.sql -t 8'
WRITE(*,*) ' ctp_run -e einternal -o  '//TRIM(votca_xml_options_filename)//' -f state.sql -t 8'
WRITE(*,*) ' ctp_run -e eoutersphere -o  '//TRIM(votca_xml_options_filename)//' -f state.sql -t 8'
WRITE(*,*) ' '
WRITE(*,*) ' ### 6) Transfer integrals'
WRITE(*,*) ' ctp_run -e izindo -o  '//TRIM(votca_xml_options_filename)//' -f state.sql -t 8'
WRITE(*,*) ' '
WRITE(*,*) ' ### 7)  Rates'
WRITE(*,*) ' ctp_run -e rates -o  '//TRIM(votca_xml_options_filename)//' -f state.sql -t 8'
WRITE(*,*) ' '
WRITE(*,*) ' ### 9) Charge Transport'
WRITE(*,*) ' # kmc_run -e kmcmultiple -o  '//TRIM(votca_xml_options_filename)//' -f state.sql '


WRITE(*,*) ' '
WRITE(*,*) ' '
WRITE(*,*) ' ### optional: pewald3D '
WRITE(*,*) ' # where required: create mps files ctp_tools -e log2mps -o options.xml  and  ctp_tools -e molpol -o options.xml '
WRITE(*,*) ' ### 10) stateserver' 
WRITE(*,*) ' # ctp_run -e stateserver -o options_stateserver.xml -f state.sql'
WRITE(*,*) ' ### 11) jobwriter '
WRITE(*,*) ' # ctp_run -e jobwriter -o options_jobwriter.xml -f state.sql '
WRITE(*,*) ' ### 12) ewald3D background polarization  '
WRITE(*,*) ' # ctp_run -e ewdbgpol -o options_ewdbgpol.xml -f state.sql -t 8 '
WRITE(*,*) ' ### 13) pewald3D '
WRITE(*,*) ' # ctp_parallel -e pewald3d -o options_pewald3d.xml -f /absolute/path/to/state.sql -s 0 -t 8 -c 8 '

WRITE(*,*) '##enumerated_gro_to_VOTCA__DIPBI_P3HT_prepare: Done'
END SUBROUTINE enumerated_gro_to_VOTCA__DIPBI_P3HT_prepare




!! create votca mapping file e.g. with virtual sites. For DIPBI/P3HT
SUBROUTINE make_votca_kmc_mapfile(map_filename,N_atoms_in_segment,residue_namesH,elementsH,atomsorteH,&
					&iHOMO_list,lambda_array,orbs_ZINDO,add_virtual_H_atoms_to_gro,&
					&exclude_virtual_H_atoms)
	IMPLICIT NONE		
	INTEGER, INTENT(IN)::N_atoms_in_segment,iHOMO_list
	CHARACTER(5), INTENT(IN), DIMENSION(N_atoms_in_segment)  ::residue_namesH
	CHARACTER(5), INTENT(IN), DIMENSION(N_atoms_in_segment)  :: elementsH,atomsorteH
	CHARACTER(5000), INTENT(IN)::map_filename
	REAL*8, DIMENSION (7), INTENT(IN)::lambda_array
	LOGICAL, INTENT(IN)::orbs_ZINDO,add_virtual_H_atoms_to_gro,exclude_virtual_H_atoms
	INTEGER::i,i_start,i_QM,ierror,iHOMO
	LOGICAL::Datei_vorhanden
	CHARACTER(5) ::mod_residue_namesH,residue_name,mod_mpsfilename
	CHARACTER(8) ::string
 	!! Strings to store the data for the atoms in the given representation
	CHARACTER(100000):: l_mdatoms,l_qmatoms,l_mpoles,l_weights,l_weights_mps,l_virtual_mps,l_seg_map
	LOGICAL::seg_ende
	l_mdatoms=''
	l_qmatoms=''
	l_mpoles=''
	l_weights=''
	l_weights_mps=''
	l_virtual_mps=''
	l_seg_map=''

	i_start=1
	seg_ende=.false.
	! Modify residue name
	mod_residue_namesH=residue_namesH(1)
	residue_name=residue_namesH(1)
	mod_mpsfilename=residue_namesH(1)
	!WRITE(*,*) 'tmp ',mod_residue_namesH,lambda_array
	IF ( (mod_residue_namesH(:2) == 'PH')   .OR. (mod_residue_namesH(:2) == 'PM') .OR. & ! P3HT / P3MT
	   & (mod_residue_namesH(1:1) == 'A')   .OR. (mod_residue_namesH(1:1)== 'B' ) )THEN  ! PBDT_TS1
		IF (mod_residue_namesH(5:) == 'M') THEN
			i_start=2
		ELSE IF (mod_residue_namesH(5:) == 'E') THEN
			i_start=2
			seg_ende=.true.
		END IF
		mod_residue_namesH=TRIM(mod_residue_namesH(:4))//'A'
		!WRITE(*,*) 'Modified:',mod_residue_namesH
	END IF


	!!!! Modify string for additional virtual atoms
	IF (add_virtual_H_atoms_to_gro) THEN 
		i_start=1
		DO i=1,N_atoms_in_segment
			IF(.NOT. exclude_virtual_H_atoms) THEN
				l_virtual_mps=TRIM(l_virtual_mps)//' 0'
			ELSE  !!! TREAT VIRTUAL H-ATOMS = 1
				IF( (i == N_atoms_in_segment) .AND. ( (residue_name(5:) == 'S') ))THEN !!! letztes H virtuell am ersten segment , S= start segment 
					l_virtual_mps=TRIM(l_virtual_mps)//' 1'
				ELSE IF ((residue_name(5:) == 'M') .AND. ((i==N_atoms_in_segment) .OR. (i==1))) THEN ! letzte H Atome am Mittelstueck 
					l_virtual_mps=TRIM(l_virtual_mps)//' 1'
				ELSE IF ((residue_name(5:) == 'E') .AND. (i==1) ) THEN
					l_virtual_mps=TRIM(l_virtual_mps)//' 1'
				ELSE ! Normalfall Aktivierung!
					l_virtual_mps=TRIM(l_virtual_mps)//' 0'
				END IF
			END IF
		END DO
	END IF


	IF(seg_ende) THEN
		i_QM=i_start
		DO i=1,N_atoms_in_segment
			l_mdatoms=TRIM(l_mdatoms)//'  1:'//TRIM(adjustl(residue_namesH(i)))//':'//TRIM(adjustl(atomsorteH(i)))
			IF( mod_residue_namesH(:2) == 'PH' .AND.  (i==N_atoms_in_segment-20) ) THEN ! modify index endsegment P3HT
				IF(add_virtual_H_atoms_to_gro) THEN
					l_qmatoms=TRIM(l_qmatoms)//'  '//TRIM(str(N_atoms_in_segment))//':'//TRIM(adjustl(elementsH(i)))
					l_mpoles=TRIM(l_mpoles)//'  '//TRIM(str(N_atoms_in_segment))//':'//TRIM(adjustl(elementsH(i)))
					l_seg_map=TRIM(l_seg_map)//'  '//TRIM(str(N_atoms_in_segment))
				ELSE
					l_qmatoms=TRIM(l_qmatoms)//'  '//TRIM(str(N_atoms_in_segment+1))//':'//TRIM(adjustl(elementsH(i)))
					l_mpoles=TRIM(l_mpoles)//'  '//TRIM(str(N_atoms_in_segment+1))//':'//TRIM(adjustl(elementsH(i)))
					l_seg_map=TRIM(l_seg_map)//'  '//TRIM(str(N_atoms_in_segment+1))
				END IF
			ELSE IF( mod_residue_namesH(:2) == 'PM' .AND.  (i==N_atoms_in_segment-5) ) THEN ! modify index endsegment P3MT
				IF(add_virtual_H_atoms_to_gro) THEN
					l_seg_map=TRIM(l_seg_map)//'  '//TRIM(str(N_atoms_in_segment)) ! Setzen des letzten Atoms
					l_qmatoms=TRIM(l_qmatoms)//'  '//TRIM(str(N_atoms_in_segment))//':'//TRIM(adjustl(elementsH(i)))
					l_mpoles=TRIM(l_mpoles)//'  '//TRIM(str(N_atoms_in_segment))//':'//TRIM(adjustl(elementsH(i)))
				ELSE
					l_qmatoms=TRIM(l_qmatoms)//'  '//TRIM(str(N_atoms_in_segment+1))//':'//TRIM(adjustl(elementsH(i)))
					l_mpoles=TRIM(l_mpoles)//'  '//TRIM(str(N_atoms_in_segment+1))//':'//TRIM(adjustl(elementsH(i)))
					l_seg_map=TRIM(l_seg_map)//'  '//TRIM(str(N_atoms_in_segment+1))
				END IF			
			ELSE
				l_qmatoms=TRIM(l_qmatoms)//'  '//TRIM(str(i_QM))//':'//TRIM(adjustl(elementsH(i)))
				l_mpoles=TRIM(l_mpoles)//'  '//TRIM(str(i_QM))//':'//TRIM(adjustl(elementsH(i)))
				l_seg_map=TRIM(l_seg_map)//'  '//TRIM(str(i_QM))
				i_QM=i_QM+1
			END IF
			l_weights=TRIM(l_weights)//'  '//TRIM(str(NINT(element_to_mass(elementsH(i),.false.))))
			l_weights_mps=TRIM(l_weights_mps)//'  '//str(NINT(element_to_mass(elementsH(i),.false.)))
			IF (.NOT. add_virtual_H_atoms_to_gro) THEN 
				l_virtual_mps=TRIM(l_virtual_mps)//' 0'
			END IF
		END DO
		! Orbitals
		IF(orbs_ZINDO) THEN
			iHOMO=iHOMO_list
		ELSE
			iHOMO=Calc_N_electrons_in(elementsH,N_atoms_in_segment,.false.)/2
		END IF ! Orbitals

	ELSE ! NORMAL
		i_QM=i_start
		DO i=1,N_atoms_in_segment
			l_mdatoms=TRIM(l_mdatoms)//'  1:'//TRIM(adjustl(residue_namesH(i)))//':'//TRIM(adjustl(atomsorteH(i)))
			l_qmatoms=TRIM(l_qmatoms)//'  '//TRIM(str(i_QM))//':'//TRIM(adjustl(elementsH(i)))
			l_mpoles=TRIM(l_mpoles)//'  '//TRIM(str(i_QM))//':'//TRIM(adjustl(elementsH(i)))
			l_seg_map=TRIM(l_seg_map)//'  '//TRIM(str(i_QM))
			l_weights=TRIM(l_weights)//'  '//TRIM(str(NINT(element_to_mass(elementsH(i),.false.))))
			l_weights_mps=TRIM(l_weights_mps)//'  '//str(NINT(element_to_mass(elementsH(i),.false.)))
			IF (.NOT. add_virtual_H_atoms_to_gro) THEN 
				l_virtual_mps=TRIM(l_virtual_mps)//' 0'
			END IF
			i_QM=i_QM+1
		END DO
		! Orbitals
		IF(orbs_ZINDO) THEN
			iHOMO=iHOMO_list
		ELSE
			iHOMO=Calc_N_electrons_in(elementsH,N_atoms_in_segment,.false.)/2
		END IF ! Orbitals
	END IF

!!! Create mapfile
INQUIRE(file=TRIM(map_filename),exist=Datei_vorhanden) ! Abfrage ob map-Datei vorhanden ist.
IF(.NOT. Datei_vorhanden) THEN 
	WRITE(*,*) 'Error File does not exist:'//TRIM(map_filename)
	CALL EXIT(1)	
ELSE
OPEN(UNIT=21,FILE=TRIM(map_filename),STATUS='unknown', access = 'append')
!WRITE(21,'(g0)') '<topology>'
!WRITE(21,'(g0)') '	<molecules>'
WRITE(21,'(g0)') '		<molecule>'
WRITE(21,'(g0)') '			<name>'//TRIM(adjustl(residue_namesH(1)))//'</name>'
WRITE(21,'(g0)') '			<mdname>'//TRIM(adjustl(residue_namesH(1)))//'</mdname>'
WRITE(21,'(g0)') '			<segments>'
WRITE(21,'(g0)') '				<segment>'
WRITE(21,'(g0)') '					<name>'//TRIM(adjustl(residue_namesH(1)))//'</name>'
WRITE(21,'(g0)') '					<qmcoords>QC_FILES/'//TRIM(adjustl(mod_residue_namesH))//'_n.xyz</qmcoords>'
WRITE(21,'(g0)') '					<orbitals>QC_FILES/'//TRIM(adjustl(mod_residue_namesH))//'.orbs</orbitals>'
WRITE(21,'(g0)') '					<basisset>INDO</basisset>'
WRITE(string,'(1I7)') iHOMO
WRITE(21,*) '					<torbital_h>'//TRIM(adjustl(string))//'</torbital_h>'
WRITE(string,'(1f6.3)') Real(lambda_array(3)-lambda_array(1))
WRITE(21,*) '					<U_cC_nN_h>'//TRIM(adjustl(string))//'</U_cC_nN_h>'
WRITE(string,'(1f6.3)') Real(lambda_array(7)-lambda_array(1))
WRITE(21,*) '					<U_nC_nN_h>'//TRIM(adjustl(string))//'</U_nC_nN_h>'
WRITE(string,'(1f6.3)') Real(lambda_array(5)-lambda_array(3))
WRITE(21,*) '					<U_cN_cC_h>'//TRIM(adjustl(string))//'</U_cN_cC_h>'
WRITE(string,'(1I7)') iHOMO+1
WRITE(21,*) '					<torbital_e>'//TRIM(adjustl(string))//'</torbital_e>'
WRITE(string,'(1f6.3)') Real(lambda_array(2)-lambda_array(1))
WRITE(21,*) '					<U_cC_nN_e>'//TRIM(adjustl(string))//'</U_cC_nN_e>'
WRITE(string,'(1f6.3)') Real(lambda_array(6)-lambda_array(1))
WRITE(21,*) '					<U_nC_nN_e>'//TRIM(adjustl(string))//'</U_nC_nN_e>'
WRITE(string,'(1f6.3)') Real(lambda_array(4)-lambda_array(2))
WRITE(21,*) '					<U_cN_cC_e>'//TRIM(adjustl(string))//'</U_cN_cC_e>'
WRITE(21,'(g0)') '					<multipoles_n>MP_FILES/'//TRIM(adjustl(mod_residue_namesH))//'_n.mps</multipoles_n>'
WRITE(21,'(g0)') '					<multipoles_h>MP_FILES/'//TRIM(adjustl(mod_residue_namesH))//'_h.mps</multipoles_h>'
WRITE(21,'(g0)') '					<multipoles_e>MP_FILES/'//TRIM(adjustl(mod_residue_namesH))//'_e.mps</multipoles_e>'
WRITE(21,'(g0)') '					<map2md>1</map2md>'
WRITE(21,'(g0)') '					<map>'//TRIM(l_seg_map)//' </map>'
WRITE(21,'(g0)') '					<weights>'//TRIM(l_weights)//' </weights>'
WRITE(21,'(g0)') '					<fragments>'
WRITE(21,'(g0)') '						<fragment>'
WRITE(21,'(g0)') '							<name>'//TRIM(adjustl(residue_namesH(1)))//'_1</name>'
WRITE(21,'(g0)') '							<mdatoms>'//TRIM(adjustl(l_mdatoms))//'</mdatoms>'
WRITE(21,'(g0)') '							<qmatoms>'//TRIM(l_qmatoms)//' </qmatoms>'
WRITE(21,'(g0)') '							<mpoles>'//TRIM(l_mpoles)//' </mpoles>'
WRITE(21,'(g0)') '							<weights>'//TRIM(l_weights)//' </weights>'
WRITE(21,'(g0)') '							<localframe>  2   3   4  </localframe>'
WRITE(21,'(g0)') '							<localframe_mps>  2   3   4  </localframe_mps>'
WRITE(21,'(g0)') '							<weights_mps>'//TRIM(l_weights_mps)//' </weights_mps>'
WRITE(21,'(g0)') '							<virtual_mps>'//TRIM(l_virtual_mps)//' </virtual_mps>'
WRITE(21,'(g0)') '						</fragment>'
WRITE(21,'(g0)') '					</fragments>'
WRITE(21,'(g0)') '				</segment>'
WRITE(21,'(g0)') '			</segments>'
WRITE(21,'(g0)') '		</molecule>'
!WRITE(21,'(g0)') '	</molecules>'
!WRITE(21,'(g0)') '</topology>'
END IF 
CLOSE(21)
END SUBROUTINE make_votca_kmc_mapfile


!! make_votca_molpol_optionsfiles
SUBROUTINE make_votca_molpol_optionsfiles(residue_name)
	IMPLICIT NONE
	CHARACTER(5),INTENT(IN)::residue_name
	CHARACTER(5000)::mps_filename,molpol_options_filename
	CHARACTER(5)::charge
	CHARACTER(7)::mol_c ! residue_name_charge
	INTEGER::i
	LOGICAL::Datei_vorhanden
DO i=1,3	
	charge=''
	IF(i==1)THEN
		charge='n'
	ELSE IF(i==2)THEN
		charge='e'
	ELSE IF(i==3)THEN
		charge='h'
	END IF
	mol_c=TRIM(residue_name)//'_'//TRIM(charge)
	mps_filename=TRIM(mol_c)//'.mps'
	INQUIRE(file='MP_FILES'//TRIM(mps_filename),exist=Datei_vorhanden)
	IF( .NOT. Datei_vorhanden)THEN
		WRITE(*,*) 'WARNING: mps file for residue name is not available: MP_FILES/'//TRIM(mps_filename)
	END IF
	molpol_options_filename='options_molpol_'//TRIM(mol_c)//'.xml'
	INQUIRE(file=TRIM(molpol_options_filename),exist=Datei_vorhanden) ! Abfrage ob molpol_options_filename-Datei vorhanden ist.
	IF(Datei_vorhanden) THEN 
		WRITE(*,*) 'File does already exist:'//TRIM(molpol_options_filename)
	ELSE
		OPEN(UNIT=25,FILE=TRIM(molpol_options_filename),STATUS='REPLACE')	
		WRITE(25,*) '<options>'
		WRITE(25,*) '   <log2mps>'
		WRITE(25,*) '      <package>gaussian</package>'
		WRITE(25,*) '      <logfile>MP_FILES/'//TRIM(mol_c)//'.log</logfile>'
		WRITE(25,*) '      <mpsfile>MP_FILES/'//TRIM(mps_filename)//'</mpsfile>'
		WRITE(25,*) '   </log2mps>'
		WRITE(25,*) '   <molpol>'
		WRITE(25,*) '      <mpsfiles>'
		WRITE(25,*) '         <input>MP_FILES/'//TRIM(mps_filename)//'</input>'
		WRITE(25,*) '         <output>output_'//TRIM(mps_filename)//'</output>'
		WRITE(25,*) '         <polar>polar_'//TRIM(mol_c)//'.xml</polar>'
		WRITE(25,*) '      </mpsfiles>'
		WRITE(25,*) '         <induction>'
		WRITE(25,*) '            <expdamp>0.39000</expdamp>'
		WRITE(25,*) '            <wSOR>0.30000</wSOR>'
		WRITE(25,*) '            <maxiter>1024</maxiter>'
		WRITE(25,*) '            <tolerance>0.00001</tolerance>'
		WRITE(25,*) '         </induction>'
		WRITE(25,*) '            <target>'
		WRITE(25,*) '            <optimize>true</optimize>'
		WRITE(25,*) '            <molpol>77 0 0 77 0 77</molpol>'
		WRITE(25,*) '            <tolerance>0.00001</tolerance>'
		WRITE(25,*) '      </target>'
		WRITE(25,*) '   </molpol>'
		WRITE(25,*) '</options>'
		CLOSE(25)
	END IF
END DO ! loop mps files charges n, e, h
END SUBROUTINE make_votca_molpol_optionsfiles


!! create a votca boxfile from the *gro imput data
SUBROUTINE make_votca_boxfile(x_box,y_box,z_box,gro_inputfile)
	! Creates the VOTCA boxfile in formated output.
	IMPLICIT NONE
	REAL ::x_box,y_box,z_box
	CHARACTER(5000)::gro_inputfile
	LOGICAL::Datei_vorhanden
	CHARACTER(5000)::boxfilename
	CHARACTER(7)::formated_boxsize
	boxfilename='box.xml'
	INQUIRE(file=TRIM(boxfilename),exist=Datei_vorhanden) ! Abfrage ob itp-Datei vorhanden ist.
	IF(Datei_vorhanden) THEN 
		WRITE(*,*) 'WARINING: File for votca box file allready exist: '&
				&//TRIM(boxfilename)
		WRITE(*,*) 'NO new file is produced: '//TRIM(boxfilename)
		CONTINUE
	ELSE
		WRITE(*,*) 'Boxsisze:',x_box,y_box,z_box
		OPEN(UNIT=24,FILE=TRIM(boxfilename),STATUS='REPLACE')	
		WRITE(24,*) '<options>'
		WRITE(formated_boxsize,'(1f7.3)') x_box 
		WRITE(24,*) '	<boxX>'//TRIM(adjustl(formated_boxsize))//'</boxX>   <!-- '//TRIM(gro_inputfile)//' -->'
		WRITE(formated_boxsize,'(1f7.3)') y_box 
		WRITE(24,*) '	<boxY>'//TRIM(adjustl(formated_boxsize))//'</boxY>'
		WRITE(formated_boxsize,'(1f7.3)') z_box 
		WRITE(24,*) '	<boxZ>'//TRIM(adjustl(formated_boxsize))//'</boxZ>'
		WRITE(24,*) '	<boxXY>0</boxXY>'
		WRITE(24,*) '	<boxXZ>0</boxXZ>'
		WRITE(24,*) '	<boxYZ>0</boxYZ>'
		WRITE(24,*) '</options>'	
		CLOSE(24)
		WRITE(*,*) 'Votca boxfile was generated: '//TRIM(boxfilename) 
	END IF
END SUBROUTINE make_votca_boxfile


!! Create a neighbourlist file, here you can modify the cutoffs for the pairs in kmc simulations
SUBROUTINE make_votca_neighbourlist_constrains(N_seg_types_all,residue_name_list,molecule_da,neighbourlist_constrains_filename)
		! Creates a file for the constrain requirements, needs to be included into the VOTCA options.xml file 
		! to build the neighbourlist in VOTCA
		IMPLICIT NONE
		INTEGER, INTENT(IN) ::N_seg_types_all
		CHARACTER(5), DIMENSION(N_seg_types_all), INTENT(IN)  :: residue_name_list
		LOGICAL, DIMENSION(N_seg_types_all), INTENT(IN)  :: molecule_da
		CHARACTER(5000)::neighbourlist_constrains_filename
		INTEGER::i,j
		LOGICAL::Datei_vorhanden
		CHARACTER(3)::NNcutoff_1,NNcutoff_2,NNcutoff_3,NNcutoff_4
	!!! Cutoff fuer nearest neighbour list, 1) P3HT /P3MT segment
	NNcutoff_1='1.5' !! 0.9
       !higher value for DPBIK
	NNcutoff_2='1.3'
       ! cutoff  PPDI
	NNcutoff_3='1.3'
       ! cutoff  PBDT-TS1
	NNcutoff_4='3.0'

	INQUIRE(file=TRIM(neighbourlist_constrains_filename),exist=Datei_vorhanden) ! Abfrage ob itp-Datei vorhanden ist.
	IF(Datei_vorhanden) THEN 
		WRITE(*,*) 'WARINING: File for neighbourlist_constrains_filename allready exist: '&
				&//TRIM(neighbourlist_constrains_filename)
		WRITE(*,*) 'NO new file is produced.'
		CONTINUE
	ELSE
	OPEN(UNIT=23,FILE=TRIM(neighbourlist_constrains_filename),STATUS='REPLACE', access = 'append')
	WRITE(23,'(g0)') '<options>'
	WRITE(23,'(g0)') '    <neighborlist>  <!--  Standard <cutoff>0.7</cutoff>   -->'
	DO j=1,N_seg_types_all
		IF(molecule_da(j)) THEN
			DO i=j,N_seg_types_all
				IF(molecule_da(i)) THEN
					WRITE(23,'(g0)') '	     <segments>'
					WRITE(23,'(g0)') '               <type> '//TRIM(adjustl(residue_name_list(j)))//&
									&' '//TRIM(adjustl(residue_name_list(i)))//'</type>'
					!!! Cutoffs for neighbour search 
					IF ((residue_name_list(j) == 'DPBIK') .OR. (residue_name_list(j) == 'DIPBI') .OR. & 
					  & (residue_name_list(i) == 'DPBIK') .OR. (residue_name_list(i) == 'DIPBI')) THEN
					       WRITE(23,'(g0)') '               <cutoff>'//TRIM(adjustl(NNcutoff_2))//'</cutoff>'
					ELSE IF( (residue_name_list(j)(1:4) == 'PPDI') .OR. (residue_name_list(i)(1:4) == 'PPDI')) THEN
						WRITE(23,'(g0)') '               <cutoff>'//TRIM(adjustl(NNcutoff_3))//'</cutoff>'
					ELSE IF( (residue_name_list(j)(1:1) == 'A') .OR. (residue_name_list(j)(1:1) == 'B') .OR. & 
						 (residue_name_list(i)(1:1) == 'A') .OR. (residue_name_list(i)(1:1) == 'B') ) THEN !!! PBDT-TS1 poleymer segment
						WRITE(23,'(g0)') '               <cutoff>'//TRIM(adjustl(NNcutoff_4))//'</cutoff>'
					ELSE ! Standard
						WRITE(23,'(g0)') '               <cutoff>'//TRIM(adjustl(NNcutoff_1))//'</cutoff>'
					END IF

					WRITE(23,'(g0)') '            </segments>'
				END IF! i molecule_da
			END DO !i
		END IF ! j molecule_da
	END DO !j
	WRITE(23,'(g0)') '	 </neighborlist>'
	WRITE(23,'(g0)') '</options>'
	WRITE(*,*) 'votca neighbourlist constrained file was generated: ',TRIM(neighbourlist_constrains_filename)	
	END IF
END SUBROUTINE make_votca_neighbourlist_constrains



!! make fake itp files that will be linked in the topology file to produce a gromacs inputfiles
SUBROUTINE make_fake_itp(N_atoms_in_segment,residue_name,elementsH,atomsorteH)
	IMPLICIT NONE
	INTEGER, INTENT(IN)::N_atoms_in_segment
	CHARACTER(5), INTENT(IN)  ::residue_name
	CHARACTER(5), INTENT(IN), DIMENSION(N_atoms_in_segment) :: elementsH,atomsorteH
	INTEGER::i,ierror
	LOGICAL::Datei_vorhanden
	CHARACTER(5000)::itp_filename
	itp_filename=TRIM(residue_name)//'_raw.itp'	
	INQUIRE(file=TRIM(itp_filename),exist=Datei_vorhanden) ! Abfrage ob itp-Datei vorhanden ist.
	IF(Datei_vorhanden) THEN 
		CONTINUE
	ELSE
		OPEN(UNIT=21,FILE=TRIM(itp_filename),STATUS='REPLACE',IOSTAT=ierror)
		WRITE(21,'(g0)') '[ moleculetype ]'
		WRITE(21,'(g0)') '; Name            nrexcl                                                        '
		WRITE(21,'(g0)') TRIM(residue_name)//'             3'
		WRITE(21,'(g0)') '                                                                                '
		WRITE(21,'(g0)') '[ atoms ]                                                                       '
		WRITE(21,'(g0)') ';       nr      type     resnr   residue      atom      cgnr    charge      mass'
		WRITE(21,'(g0)') '; res '//TRIM(adjustl(residue_name))//'_1'     
		!!!    FORMATIERTE AUSGABE
		!!	WRITE(21,*) '         1         H         1     P3MTA        H1         1     0.000     1.000'
		DO i=1,N_atoms_in_segment
			WRITE(21,'(I10,A10,I10,2A10,I10,2f10.3)') i,ADJUSTR(TRIM(elementsH(i))),1,ADJUSTR(TRIM(residue_name)),&
			&ADJUSTR(TRIM(atomsorteH(i))),i,0.0,1.0
		END DO                                                              
		WRITE(21,'(g0)') ' '
	END IF

END SUBROUTINE make_fake_itp

!! read lambda_in_all_inputfile
SUBROUTINE read_lambda_in_all_inputfile(N_types,array,lambda_in_all_inputfile)
		!! Reads lambda in data in eV!!!!
		!!	R1=OptS0_n_${monomerA}	R2=OptS0_el_${monomerA}	R3=OptS0_lo_${monomerA}	R7=SCF_el_OptS0_n_${monomerA}	R8=SCF_lo_OptS0_n_${monomerA}	R9=SCF_n_OptS0_el_${monomerA}	R10=SCF_n_OptS0_lo_${monomerA}
		IMPLICIT NONE
		CHARACTER(5000), INTENT(IN)::lambda_in_all_inputfile
		REAL*8, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) ::array
		INTEGER, INTENT(OUT)::N_types
		CHARACTER (5):: residue_name
		CHARACTER (500)::dummy
		INTEGER::i,j,N_lines,N_skip_line,ierror
		LOGICAL::Datei_vorhanden

		INQUIRE(file=TRIM(lambda_in_all_inputfile),exist=Datei_vorhanden) ! Abfrage ob Datei vorhanden ist.
		IF( .NOT. Datei_vorhanden) THEN
		      WRITE(*,*) ' Fehler: Die_Datei_ist_nicht_vorhanden. '//TRIM(lambda_in_all_inputfile)
		      WRITE(*,*) ' ENDE'
		      CALL EXIT (1)
		END IF ! Datei-file vorhanden  

		OPEN(UNIT=20,FILE=TRIM(lambda_in_all_inputfile),STATUS='OLD',IOSTAT=ierror)
		WRITE(*,*) 'Reading data from: ',TRIM(lambda_in_all_inputfile)
		ierror=0
		N_skip_line=0
		N_lines=0
		DO WHILE ( ierror == 0 ) 
		    Read(20,*,IOSTAT=ierror) dummy
		    IF(ierror < 0) EXIT
		    IF(ierror > 0) THEN
			WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(lambda_in_all_inputfile),' aufgetreten.'
			WRITE(*,*) 'Beende Einlesen !'
			CALL EXIT (1)
		    END IF
		    IF ( dummy(1:1) == '#' ) THEN ! Check ob erster Buchstabe ein #-Zeichen ist, um kommentare zu lesen
			N_skip_line=N_skip_line+1
		    ELSE
			N_lines=N_lines+1
		    END IF
		END DO       
		CLOSE(20)        
		WRITE(*,*) 'N_lines=',N_lines,'N_skip_line',N_skip_line        

		OPEN(UNIT=20,FILE=TRIM(lambda_in_all_inputfile),STATUS='OLD',IOSTAT=ierror)
		IF (N_skip_line > 0 ) THEN
		    DO i =1,N_skip_line 
			Read(20,*,IOSTAT=ierror) dummy
		    END DO
		END IF
		
		N_types=N_lines
		ALLOCATE(array(N_types,7))		

		DO i=1,N_types
		    Read(20,*,IOSTAT=ierror) residue_name,array(i,1),array(i,2),array(i,3),array(i,4),array(i,5), & 
			    &  array(i,6),array(i,7) !,array(i,8),array(i,9),array(i,10), &
			    !&  array(i,11),array(i,12),array(i,13),array(i,14),array(i,15), &
			    !&  array(i,16),array(i,17),array(i,18),array(i,19),array(i,20), &
			    !&  array(i,21),array(i,22),array(i,23),array(i,24),array(i,25), &
			    !&  array(i,26),array(i,27),array(i,28),array(i,29),array(i,30), &
			    !&  array(i,31),array(i,32),array(i,33)
		   	    !!! !WRITE(*,*) array(i,1),array(i,2),array(i,3),array(i,4),array(i,5)
		    IF(ierror < 0) EXIT
		    IF(ierror > 0) THEN
			WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von '//TRIM(lambda_in_all_inputfile)//' aufgetreten.'
			WRITE(*,*) 'Beende Einlesen !'
			CALL EXIT (1)
		    END IF
		END DO
		CLOSE(20)

END SUBROUTINE read_lambda_in_all_inputfile

!! make gromp mdp file, for an infinitisimal set in gromacs to check if all the inputfiles match the gromacs format.
SUBROUTINE make_gromp_mdp_file()
	IMPLICIT NONE
	CHARACTER (LEN=500)::gromp_filename
	LOGICAL::Datei_vorhanden
	INTEGER::ierror
	!!! Creates a minimal gromp.mdp file for gmx for an infinitisimal step
	gromp_filename='gromp.mdp'
	inquire(file=TRIM(gromp_filename),exist=Datei_vorhanden) ! Abfrage ob gro-Datei vorhanden ist.
	IF(Datei_vorhanden) THEN 
		    		WRITE(*,*) ' Warning: Die_VOTCA_gromp.mdp_Datei_is_bereits_vorhanden. '//TRIM(gromp_filename)
				WRITE(*,*) ' Stellen_Sie_sicher_dass_diese_nicht_ueberschrieben_wird! '
		    		WRITE(*,*) ' ENDE'
		    		CALL EXIT (1)
	ELSE
			! Erstell gromp_filename
		OPEN(UNIT=27,FILE=TRIM(gromp_filename),STATUS='REPLACE',IOSTAT=ierror)
		WRITE(27,*) 'integrator               = md       '
		WRITE(27,*) 'dt                       = 0.0002   '
		WRITE(27,*) 'nsteps                   = 1        '
		WRITE(27,*) 'nstxout                  = 1        '
		WRITE(27,*) 'nstvout                  = 1        '
		WRITE(27,*) 'nstlog                   = 1        '
		WRITE(27,*) 'nstenergy                = 300      '
		WRITE(27,*) 'nstxout-compressed       = 500      '
		WRITE(27,*) 'nstlist                  = 10       '
		WRITE(27,*) 'ns-type                  = grid     '
		WRITE(27,*) 'rlist                    = 0.8      '
		WRITE(27,*) 'coulombtype              = cut-off  '
		WRITE(27,*) 'rcoulomb                 = 1.4      '
		WRITE(27,*) 'rvdw                     = 1.4      '
		WRITE(27,*) 'tcoupl                   = V-rescale'
		WRITE(27,*) 'tc-grps                  = system   '
		WRITE(27,*) 'tau-t                    = 0.1      '
		WRITE(27,*) 'ref-t                    = 300      '
		WRITE(27,*) 'gen-vel                  = yes      '
		WRITE(27,*) 'gen-temp                 = 300      '
		WRITE(27,*) 'gen-seed                 = 173529   '
		WRITE(27,*) 'constraints              = all-bonds'

	END IF ! gromp-file vorhanden  


END SUBROUTINE make_gromp_mdp_file



SUBROUTINE make_geo_orbs_mps_g09(modus,NAtoms,residue_name,elements,koord,&
							& trajstep,Numcharged,NCPUS,Mem)
	!! creates the QM input geometries, which are needed to generate the name1_n.xyz, ZINDO name1.orbs files name1.mps files.
	!! extra_files_modus / modus = "all" selects the evaluation or one method.
	!! modus = [all, make_g09_for_mps, make_ZINDO, make_xyz, make_folders]
	!! make_g09_for_mps: creates g09 inputfiles for the geometries to evaluate with data for *.mps files 
	!! make_ZINDO: creates g09 input files for the ZINDO calculations, which are needed for *.orbs
 	!! make_xyz: creates name1_n.xyz in a QC_FILE folder as needed for VOTCA calculations.
 	!! make_folders: creates folders QC_FILES and MP_FILES, if they are not present in the current directory.
	!!
	!! g09: select: e.g. trajstep=0, Numcharged=0, NCPUS=8, Mem=8
	IMPLICIT NONE
	INTEGER, INTENT (IN):: NAtoms,trajstep,Numcharged,NCPUS,Mem
	CHARACTER(LEN=5), DIMENSION(NAtoms), INTENT (IN)  :: elements
	REAL, DIMENSION(NAtoms,3), INTENT (IN)            :: koord
	CHARACTER(LEN=5)  , INTENT (IN)                   :: residue_name
	CHARACTER(LEN=500), INTENT (IN)                   :: modus
	!!! Local variables
	CHARACTER(LEN=500) :: foldername,xyz_filename,ZINDO_orbs_filename,mps_filename,filename
	CHARACTER(LEN=10000) :: bashline
	Character(500) :: molA,g09_inputline,calculation_method
	Character(5)   ::mod_residue_name
	Character(1)   ::charge_name
	LOGICAL :: folder_vorhanden,file_exists,make_folders,make_g09_inp_for_mps,make_g09_inp_ZINDO,make_xyz
	INTEGER :: i,charge

	folder_vorhanden=.false.
	make_xyz=.false.
	make_folders=.false.
	make_g09_inp_ZINDO=.false.
	make_g09_inp_for_mps=.false.

	IF( (INDEX(TRIM(modus), 'false') /= 0) .OR. (INDEX(TRIM(modus), 'False') /= 0) )THEN  
		RETURN
	ELSE IF( TRIM(modus) == 'all') THEN
		make_xyz=.true.
		make_folders=.true.
		make_g09_inp_ZINDO=.true.
		make_g09_inp_for_mps=.true.
	ELSE IF ( TRIM(modus) == 'make_g09_for_mps') THEN 
		make_g09_inp_for_mps=.true.
	ELSE IF( TRIM(modus) == 'make_ZINDO') THEN 
		make_g09_inp_ZINDO=.true.
	ELSE IF( TRIM(modus) == 'make_xyz') THEN
		make_xyz=.true.
	ELSE IF( TRIM(modus) == 'make_folders') THEN
		make_folders=.true.
	ELSE
		WRITE(*,*) 'ERROR: Modus not found.'//TRIM(modus)
		CALL EXIT (1)
	END IF ! modus ?
	
	mod_residue_name=TRIM(residue_name)
	IF ( (mod_residue_name(:2) == 'PH')   .OR. (mod_residue_name(:2) == 'PM') .OR. & ! P3HT / P3MT
	   & (mod_residue_name(1:1) == 'A')   .OR. (mod_residue_name(1:1)== 'B' ) )THEN  ! PBDT_TS1
		mod_residue_name=mod_residue_name(:4)//'A'
		!WRITE(*,*) 'Modified:',mod_residue_name
	END IF




	IF(make_folders) THEN
		foldername='QC_FILES'
                INQUIRE(file=TRIM(foldername),exist=folder_vorhanden)               
                IF (.NOT. folder_vorhanden) THEN
                        write(bashline,*) 'mkdir '//TRIM(foldername)
                        CALL execute_command_line(TRIM(bashline)) 
			WRITE(*,*) TRIM(bashline)
                END IF ! folder 
                
                foldername='MP_FILES'
                INQUIRE(file=TRIM(foldername),exist=folder_vorhanden)               
                IF (.NOT. folder_vorhanden) THEN
                        write(bashline,*) 'mkdir '//TRIM(foldername)
                        CALL execute_command_line(TRIM(bashline)) 
			WRITE(*,*) TRIM(bashline)
                END IF ! folder      
                make_folders=.false.          
	END IF ! make folders
	
	IF(make_xyz) THEN
		foldername='QC_FILES'
                INQUIRE(file=TRIM(foldername),exist=folder_vorhanden)               
                IF (folder_vorhanden) THEN
					molA=TRIM(mod_residue_name)//'_n'
					xyz_filename=TRIM(foldername)//'/'//TRIM(molA)//'.xyz'
					INQUIRE(file=TRIM(xyz_filename),exist=file_exists)
					IF( .NOT. file_exists) THEN
						call make_coord_to_xyz_in_dir(elements,koord,NAtoms,molA,foldername)
						WRITE(*,*) "New file created: "//TRIM(xyz_filename)
					END IF 
		ELSE
					WRITE(*,*) 'Warning: QC_FILES folder does not exist. No *.xyz file produced.'
		END IF
	END IF ! make_xyz
	
	IF(make_g09_inp_ZINDO) THEN
		foldername='QC_FILES'
                INQUIRE(file=TRIM(foldername),exist=folder_vorhanden)               
                IF (folder_vorhanden) THEN
					ZINDO_orbs_filename=TRIM(mod_residue_name)//'.orbs'
					filename=TRIM(foldername)//'/'//TRIM(ZINDO_orbs_filename)
					INQUIRE(file=TRIM(filename),exist=file_exists)
					IF( .NOT. file_exists) THEN
						molA=TRIM(mod_residue_name)
						g09_inputline=" ZINDO  punch=mo  iop(7/33=1) SCF(XQC,MaxConventionalCycles=400,MaxCycle=800) "
						g09_inputline=TRIM(g09_inputline)//" pop=minimal  iop(6/80=1)  nosymm IOp(6/7=3)"
						calculation_method="ZINDO"
						CALL make_coord_to_g09_charged_inp(elements,koord,NAtoms,calculation_method,&
							&g09_inputline,molA,0,0,0,NCPUS,Mem)
						WRITE(*,*) "New g09_ZINDO.inp file created for: "//TRIM(molA)
					END IF ! orbs file_exists
                ELSE
					WRITE(*,*) 'Warning: QC_FILES folder does not exist. No g09.inp file produced.'
		END IF !  folder_vorhanden
	END IF ! make_g09_inp_ZINDO
	
	
	IF(make_g09_inp_for_mps) THEN
		foldername='MP_FILES'
                INQUIRE(file=TRIM(foldername),exist=folder_vorhanden)               
                IF (folder_vorhanden) THEN
					DO i=-1,1,1  ! charge states for mps files
						IF(i ==0) charge_name='n'
						IF(i ==1) charge_name='h'
						IF(i ==-1) charge_name='e'
						
						mps_filename=TRIM(mod_residue_name)//'_'//TRIM(charge_name)//'.mps'
						filename=TRIM(foldername)//'/'//TRIM(mps_filename)
						INQUIRE(file=TRIM(filename),exist=file_exists)
						IF( .NOT. file_exists) THEN
							molA=TRIM(mod_residue_name)
							g09_inputline=" PBE1PBE/6-31G*  pop=minimal SCF(XQC,MaxConventionalCycles=400,MaxCycle=800) " 
							g09_inputline=TRIM(g09_inputline)//"   pop=(minimal,CHelpG) nosymm punch=mo IOp(6/7=3) "
							calculation_method="PBE0_6-31Gs"
							CALL make_coord_to_g09_charged_inp(elements,koord,NAtoms,calculation_method,&
								&g09_inputline,molA,0,0,i,NCPUS,Mem)
							WRITE(*,*) "New g09_to_mps.inp file created for: "//TRIM(molA)
						END IF ! orbs file_exists
					END DO  ! charge states for mps files
		ELSE
						WRITE(*,*) 'Warning: MP_FILES folder does not exist. No g09.inp file produced.'					
		END IF !  folder_vorhanden			
	END IF ! make_g09_for_mps
	
	
END SUBROUTINE make_geo_orbs_mps_g09




SUBROUTINE make_coord_to_g09_charged_inp(atomsorte,koord,MAXATM,calculation_method,&
                &g09_inputline,molA,trajstep,Numcharged,charge,NCPUS,Mem)
    !! creates g09.inp files for a given geometry
    IMPLICIT NONE
    Integer, INTENT (IN):: MAXATM,charge,trajstep,Numcharged,NCPUS,Mem
    Real, Dimension(MAXATM,3), INTENT (IN) :: koord
    Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
    Character(500), INTENT (IN) :: molA
    Character(500),INTENT (IN) :: calculation_method,g09_inputline
    Character(500) ::filename
    Character(1024) ::bashline
    Character(2) ::charge_name
    Integer :: i,ierror,multiplicity
    logical::folder_vorhanden
    LOGICAL::file_exists=.false.
    
    IF(charge < 0) THEN
        charge_name='el'
        multiplicity=2
    ELSE IF(charge == 0 ) THEN    
        charge_name='n'
        multiplicity=1
    ELSE
        charge_name='lo'
        multiplicity=2
    END IF
    filename='SCF_'//TRIM(molA)//'_'//TRIM(charge_name)//'_'//TRIM(calculation_method)&
                       &//'_G'//TRIM(str(trajstep))//'_N'//TRIM(str(Numcharged))//'.inp'
    
    INQUIRE(file=TRIM(filename),exist=file_exists)
    IF (.NOT. file_exists) THEN
		OPEN(UNIT=78,FILE=TRIM(filename),STATUS='REPLACE',IOSTAT=ierror)
		  
		write(78,*) '%NProcShared='//TRIM(str(NCPUS))
		write(78,*) '%Mem='//TRIM(str(Mem))//'GB'
		write(78,*) '# '//TRIM(g09_inputline)//' NoSymm ' 
		write(78,*) '# GFINPUT 6D 10F' 
		write(78,*) ' '
		write(78,*) ' with '//TRIM(filename) 
		write(78,*) ' '
		write(78,*) TRIM(str(charge))//'  '//TRIM(str(multiplicity))
				DO i=1,MAXATM
						write(78,'(2X,A3,1X,3F16.6)') atomsorte(i),koord(i,1),koord(i,2),koord(i,3)
				END DO 
		write(78,*) ' '
		CLOSE(78)
	END IF ! file exists

END SUBROUTINE make_coord_to_g09_charged_inp 







!! Checks if the user selected filename has git a usefull termination, else => programm termintates with an error
LOGICAL FUNCTION check_filetermination(filename,termination)
	IMPLICIT NONE
	! Checks if the filename has got the correct termination
	! filename=my_molecule.xyz  and termination=.xyz => CALL EXIT (1) program with an error.
	CHARACTER (LEN=*), INTENT(IN) ::filename,termination
	LOGICAL::termination_ok
	termination_ok=.false.
	IF(index(filename((LEN_TRIM(termination)):),TRIM(termination)) == 0 ) THEN
		WRITE(*,*) 'Fehler: Eine falsche Endung der Datei '//TRIM(filename)//' ist nicht: '//TRIM(termination)
		CALL EXIT (1)
	ELSE
		termination_ok=.true.
	END IF
	check_filetermination=termination_ok
END FUNCTION


!! Integer to string conversion with WRITE statement
character(len=30) function str(k) ! Integer_to_String
IMPLICIT NONE
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

character(len=30) function strReal(k) ! Real_to_String
IMPLICIT NONE
!   "Convert an real to string."
    REAL, intent(in) :: k
    write (strReal,'(F10.6)') k
    strReal = adjustl(strReal)
end function strReal


!! make_lio create an string with RANDOM_NUMBER, in order_to_create_unique string:  lio=chaos in espanol !
CHARACTER(len=5) function make_lio() 
       ! Creates a random string unsing clock-time 
       IMPLICIT NONE
       !   "Convert an real to string."
       character(len=10) :: k
       INTEGER::i
       REAL::HARVEST
       integer, allocatable :: seed(:)
       integer:: size2
       LOGICAL::first_call=.true.
	
	IF (first_call) THEN
    		CALL init_random_seed()
		first_call=.false.
	END IF

    	call random_seed(size=size2)
    	allocate(seed(size2))
    	! set seed(:) somehow
       call random_seed(put=seed)
       call system_clock(i)
       CALL RANDOM_NUMBER(HARVEST)
       !WRITE(*,*) HARVEST,abs(mod(HARVEST,REAL(i)))*1.E5
       write (k,'(I5)') INT(abs(mod(HARVEST,REAL(i)))*1.E5)
       !WRITE(*,*) k
       make_lio=adjustl(TRIM(k))
end function make_lio

!! Random seed 
subroutine init_random_seed()
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
  call random_seed(size = n)
  allocate(seed(n))
  call system_clock(count=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
  deallocate(seed)
end


!! Add H atoms to the chain or cuts the chain if needed.
!! KETTE = .true. selects the full side chain representation.
SUBROUTINE add_h_atoms_to_mol(atomsorte,koord,NAtoms,residue_name,KETTE,atomsorte1H,koord1H,NAtoms1H,Molname) 
    IMPLICIT NONE
    Character(5), Dimension(NAtoms), INTENT(INOUT) :: atomsorte  ! , ALLOCATABLE ??
    Real, Dimension(NAtoms,3), INTENT(IN)       :: koord
    Character(5), INTENT(IN) :: residue_name
    INTEGER, INTENT(IN)::NAtoms
    LOGICAL, INTENT(IN)::KETTE
    Character(5), ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorte1H
    Real, ALLOCATABLE, Dimension(:,:), INTENT(OUT) :: koord1H
    INTEGER, INTENT(OUT) ::NAtoms1H
    !! 
    INTEGER::i,j
    LOGICAL::DEBUG=.false.
    LOGICAL::  DIPBI,DIPBI_KETTE,P3HT,P3MT,P3HT_KETTE,PPDI,PPDI_KETTE,PBDT_TS1,HDI_KETTE,PPDI_Propyl
    Character(500), INTENT(OUT)::Molname
    DIPBI=.false.
    DIPBI_KETTE=.true.
    P3HT=.false. 
    P3MT=.false.
    P3HT_KETTE=.false.
    PPDI=.false.
    PPDI_Propyl=.false.
    PPDI_KETTE=.false.
    PBDT_TS1=.false.
    HDI_KETTE=.false. 

    IF(DEBUG) THEN
	DO i=1,NAtoms
		WRITE(*,*) i,atomsorte(i),koord(i,1),koord(i,2),koord(i,3)
	END DO
     END IF

IF ('THP'==TRIM(residue_name(1:3))) THEN        
    P3HT=.TRUE.
    Molname="P3MT"
    P3MT=.TRUE.
    IF(KETTE) THEN
        Molname="P3HT"
        P3MT=.false.
	P3HT_KETTE=.true.
    END IF    	
ELSE IF( TRIM(residue_name)=='DIPBI' ) THEN
    DIPBI=.TRUE.
    Molname="DIPBI"
    IF(KETTE) THEN
        Molname="DIPBI_KETTE"
        DIPBI_KETTE=.true.
    END IF
ELSE IF ( 'HDI__'==TRIM(residue_name(1:5)) ) THEN
    Molname='HDI'          
ELSE IF ( 'HDI_K'==TRIM(residue_name(1:5)) ) THEN
    Molname='HDI_K' 
    IF(KETTE) THEN
        HDI_KETTE=.true.
        Molname="HDI_KETTE"
    END IF
ELSE IF ( 'TDI__'==TRIM(residue_name(1:5)) ) THEN
    Molname='TDI'
ELSE IF ( 'K1K__'==TRIM(residue_name(1:5)) .OR. 'PtK1K'==TRIM(residue_name(1:5)) ) THEN
    Molname='PtK1K'                     
ELSE IF ( 'K2K__'==TRIM(residue_name(1:5)) .OR. 'PtK2K'==TRIM(residue_name(1:5)) ) THEN
    Molname='PtK2K'                   
ELSE IF ( 'K3K__'==TRIM(residue_name(1:5)) .OR. 'PtK3K'==TRIM(residue_name(1:5)) ) THEN
    Molname='PtK3K'
ELSE IF ( 'K4K__'==TRIM(residue_name(1:5)) .OR. 'PtK4K'==TRIM(residue_name(1:5)) ) THEN
    Molname='PtK4K'                     
ELSE IF ( 'mCP__'==TRIM(residue_name(1:5)) ) THEN
    Molname='mCP'        
ELSE IF ( 'NPB__'==TRIM(residue_name(1:5)) ) THEN
    Molname='NPB'
ELSE IF ( 'PPDI_'==TRIM(residue_name(1:5)) ) THEN
    Molname='PPDI'   
    PPDI=.true.
    IF (KETTE) THEN
        PPDI_KETTE=.true.
    END IF
ELSE IF ( '8poly'==TRIM(residue_name(1:5)) ) THEN
    Molname='PBDT_TS1'
    PBDT_TS1=.true.
ELSE
    Molname=TRIM(residue_name)
END IF 



IF(DIPBI) THEN
    IF(DIPBI_KETTE) THEN
        call h_diPBI_kette(atomsorte,koord,atomsorte1H,koord1H,NAtoms,Molname)
    ELSE
        call h_diPBI(atomsorte,koord,atomsorte1H,koord1H,NAtoms,Molname)
    END IF
ELSE IF (P3HT) THEN
    IF(P3MT) THEN
        CALL h_P3MT(atomsorte,koord,atomsorte1H,koord1H,NAtoms,Molname)
    ELSE IF (P3HT_KETTE) THEN	

        CALL h_P3HT(atomsorte,koord,atomsorte1H,koord1H,NAtoms,Molname)  
	!DO i=1,size(atomsorte1H)
	!	WRITE(*,*) atomsorte1H(i),koord1H(i,1),koord1H(i,2),koord1H(i,3)
	!END DO
	!CALL EXIT (1)	
    END IF
ELSE IF (PPDI) THEN
    IF(PPDI_KETTE) THEN
        CALL h_PPDI_kette(atomsorte,koord,atomsorte1H,koord1H,NAtoms,Molname)
    ELSE IF (PPDI_Propyl) THEN    
        CALL h_PPDI_propyl(atomsorte,koord,atomsorte1H,koord1H,NAtoms,Molname)
    ELSE
        !DO i=1,NAtoms
        !    WRITE(*,*) atomsorte(i),koord(i,:)
        !END DO
        CALL h_PPDI(atomsorte,koord,atomsorte1H,koord1H,NAtoms,Molname)   
    END IF
ELSE IF(PBDT_TS1) THEN
    !DO i=1,NAtoms
    !    WRITE(*,*) atomsorte(i),koord(i,:)
    !END DO
    CALL h_PBDT_TS1(atomsorte,koord,atomsorte1H,koord1H,NAtoms,Molname)    
ELSE IF( (.NOT. HDI_KETTE) .AND. (TRIM(Molname) =='HDI' .OR. TRIM(Molname) =='HDI_K' .OR. TRIM(Molname) =='HDI_KETTE')) THEN
    !DO i=1,NAtoms
        !WRITE(*,*) atomsorte(i),koord(i,:)
    !END DO
    CALL h_cut_chain_HDI(atomsorte,koord,NAtoms,atomsorte1H,koord1H,Molname)
ELSE IF (index(TRIM(Molname),'Alq') .ne. 0) THEN ! H-Atoms are not changed
    CALL alq3_rename_elements(atomsorte,koord,NAtoms,atomsorte1H,koord1H,Molname)
ELSE
    
    ! Uebergabe ohne Bearbeitung der H-Atome
    ALLOCATE(koord1H(NAtoms,3))
    ALLOCATE(atomsorte1H(NAtoms))
    DO i=1,NAtoms
        koord1H(i,:)=koord(i,:)
        atomsorte1H(i)=atomsorte(i)
    END DO    
    write(*,*) ' Keine Bearbeitung der H-Atome'
END IF

NAtoms1H=size(atomsorte1H)
IF (DEBUG) THEN
    !WRITE(*,*) 'Koordinaten molA: atomsorte1H, koord1H(i,:)'
    !DO i=1,NAtoms1H
    !    WRITE(*,*) atomsorte1H(i),koord1H(i,:)
    !END DO
    WRITE(*,*) 'molA: '//TRIM(Molname)//'  ,PBDT_TS1:',PBDT_TS1,' PPDI:',PPDI,' DIPBI:',DIPBI,' P3HT:',P3HT
END IF


END SUBROUTINE add_h_atoms_to_mol



!! Add H Atoms to a P3HT molecule
SUBROUTINE h_P3HT(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel)
IMPLICIT NONE
Integer, INTENT (IN):: MAXATM
Real, Dimension(MAXATM,3), INTENT (IN) :: koord
Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
Character(5), ALLOCATABLE, Dimension(:), INTENT (OUT) :: atomsorteH
Character(5), ALLOCATABLE, Dimension(:)::atomsortedummy
Character(500), INTENT (IN) :: ziel
Integer ::i,j,k,ind,anzahl,anzahlneu,anzahlringe,ind5RingEnde,ierror1,cAtomRingInd,sAtomInd,zaehler
Character(3) :: dummy
Real, Dimension(3) :: cAtom,cAtomRing,hAtom1,hAtom2,hAtom3,hAtomOld,vorigesC,naechstesC,kreuzprod,c1Atom,c2Atom,sAtom
Integer, Dimension(MAXATM) :: ind5Ring
LOGICAL::HStart,HEnde,DEBUG=.false.

IF(DEBUG) THEN
	WRITE(*,*) 'Koordinaten vor h_P3HT ' 
	DO i=1,MAXATM,1
        	write(*,*)  atomsorte(i),(koord(i,j),j=1,3)
	END DO
END IF



HStart=.false.
HEnde=.false.

! Bestimmung ob erstes und letztes H-Atom existiert

IF( TRIM(adjustl(atomsorte(1))) == 'HC' .AND. TRIM(adjustl(atomsorte(2))) == 'CA' ) HStart=.true.
IF( TRIM(adjustl(atomsorte(MAXATM-7))) == 'HC' ) HEnde=.true.

!Indexarray fur 5Ringe und Wasserstoffe erzeugen
anzahl=MAXATM
i =1
ind=1
ind5Ring =0
DO WHILE(ind <= anzahl)
	ind5Ring(i) = ind
	IF(index(atomsorte(ind),'S') >0) THEN !Hexylgruppe ÃŒberspringen
		ind = ind +6 !jetzt: ind letztes Atom der Hexylgruppe
	END IF
	i=i+1
	ind=ind+1
END DO
ind5RingEnde = i-1 !letzter belegter Platz

!Write(*,*) 'HStart',HStart,'HEnde',HEnde

IF( HStart .AND. HEnde ) THEN
    anzahlringe = (ind5RingEnde -2)/6                   !anzahlringe = (ind5RingEnde -2)/6 !fuer ganze Kette 
    anzahlneu = anzahl + anzahlringe*13                 !Gesamtanzahl Atome in neuer Datei
else
    IF(HStart .OR. HEnde) THEN
        anzahlringe = (ind5RingEnde -1)/6  
        anzahlneu = anzahl + anzahlringe*13 + 1 !Gesamtanzahl Atome in neuer Datei +1 fuer Anfangs/Abschluss H-Atom
    ELSE 
        anzahlringe = (ind5RingEnde)/6 
        anzahlneu = anzahl + anzahlringe*13 + 2 !Gesamtanzahl Atome in neuer Datei +2 fuer Anfangs+Abschluss H-Atom
    END IF
END IF
IF(DEBUG) THEN
	WRITE(*,*) "MAXATM",MAXATM
	WRITE(*,*) "anzahlringe",anzahlringe
	WRITE(*,*) "anzahlneu",anzahlneu
	WRITE(*,*) "ind5RingEnde",ind5RingEnde
END IF !DEBUG
ALLOCATE(koordH(anzahlneu,3))
ALLOCATE(atomsorteH(anzahlneu))
ALLocate(atomsortedummy(anzahlneu))
atomsortedummy=' '
ind =1
i=1
DO WHILE(i <= MAXATM )
     dummy = adjustl(atomsorte(i))
     atomsorteH(ind) = dummy(1:1)//'  '
     atomsortedummy(i) = dummy(1:1)//'  '    !Schneidet hintere von CH das H ab / laesst nur noch ersten Buchstaben durch
     IF(index(atomsorte(i),'S') >0) THEN
        ind = ind+19                ! Ueberspringen um 19 Atome =6*C+13*H
        i=i+6
     END IF
    ind=ind+1
    i=i+1
END DO 

k=1       ! Index fuer das Speichern der Koordinaten in KoordH
!neue Koordinaten schreiben
DO i=1,ind5RingEnde,1
	!Write(23,'(2X,A2,1X,3F16.6)') atomsortedummy(ind5Ring(i)),koord(ind5Ring(i),1),koord(ind5Ring(i),2),koord(ind5Ring(i),3)    
    koordH(k,:)=koord(ind5Ring(i),:) !Uebergabe der Koordinaten an koordH
	k=k+1

    IF(index(atomsortedummy(ind5Ring(i)),'S') >0) THEN !H-Atome an Hexylgruppe anfangen
		cAtomRingInd = ind5Ring(i) -2 !Indizes zwischenspeichern
		sAtomInd = ind5Ring(i)
		
		!1. C-Atom mit H bestuecken
		DO ind=1,3,1 !Koordinaten merken
			cAtom(ind) = koord(sAtomInd+1,ind)
			vorigesC(ind) = koord(cAtomRingInd,ind)
			naechstesC(ind) = koord(sAtomInd+2,ind)
		END DO
				
		CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
		hAtom1 = cAtom + kreuzprod
		hAtom2 = cAtom - kreuzprod

		CALL rotate(naechstesC,vorigesC,hAtom1,32.) !Winkel anpassen
		CALL rotate(naechstesC,vorigesC,hAtom2,-32.)

		CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen
		CALL setDistance(cAtom,hAtom2,1.095) !Abstand anpassen
		
		!Write(23,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
        koordH(k,:)=cAtom(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='C ' 
        k=k+1
        !Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
        koordH(k,:)=hAtom1(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='H ' 
        k=k+1
		!Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
        koordH(k,:)=hAtom2(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='H ' 
        k=k+1
		!C-Atome 2 - 5 bestuecken
		DO zaehler=2,5,1
			DO ind=1,3,1 !Koordinaten merken
				cAtom(ind) = koord(sAtomInd+zaehler,ind)
				vorigesC(ind) = koord(sAtomInd+zaehler-1,ind)
				naechstesC(ind) = koord(sAtomInd+zaehler+1,ind)
			END DO
					
			CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
			hAtom1 = cAtom + kreuzprod
			hAtom2 = cAtom - kreuzprod
			
			CALL rotate(naechstesC,vorigesC,hAtom1,32.) !Winkel anpassen
			CALL rotate(naechstesC,vorigesC,hAtom2,-32.)
	
			CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen
			CALL setDistance(cAtom,hAtom2,1.095) !Abstand anpassen

			!Write(23,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
            koordH(k,:)=cAtom(:) !Uebergabe der Koordinaten an koordH
            atomsorteH(k)='C ' 
            k=k+1           
			!Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
            koordH(k,:)=hAtom1(:) !Uebergabe der Koordinaten an koordH
            atomsorteH(k)='H ' 
            k=k+1           
			!Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
            koordH(k,:)=hAtom2(:) !Uebergabe der Koordinaten an koordH
            atomsorteH(k)='H ' 
            k=k+1  
		END DO

		!letztes C-Atom mit 3 H-Atomen bestuecken
		DO ind=1,3,1 !Koordinaten merken
			cAtom(ind) = koord(sAtomInd+6,ind)
			vorigesC(ind) = koord(sAtomInd+5,ind)
			naechstesC(ind) = koord(sAtomInd+4,ind) !vorvoriges C-atom
		END DO		
		
		CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)

		hAtom1 = vorigesC
		CALL setDistance(cAtom,hAtom1,1.095)
		CALL rotate(cAtom,cAtom+kreuzprod,hAtom1,112.)
        	hAtom2 = hAtom1  !schon mal die ungefaehre Position
		hAtom3 = hAtom1

		CALL rotate(vorigesC,cAtom,hAtom2,120.) !in korrekte Position bringen
		CALL rotate(vorigesC,cAtom,hAtom3,240.)
        
        !Rotation um C-C Achse sodass der Diederwinkel H-C-C-H ~60° wird
        ! (Korrektur zur passenden Newman-Projektion)
           hAtom1=hAtom1-cAtom(:)
           call Rotationsmatrix_alpha(hAtom1,60.,vorigesC(:)-cAtom(:),.true.)
           hAtom1(:)=hAtom1(:)+cAtom(:)

           hAtom2=hAtom2-cAtom(:)
           call Rotationsmatrix_alpha(hAtom2,60.,vorigesC(:)-cAtom(:),.true.)
           hAtom2(:)=hAtom2(:)+cAtom(:)

           hAtom3=hAtom3-cAtom(:)
           call Rotationsmatrix_alpha(hAtom3,60.,vorigesC(:)-cAtom(:),.true.)
           hAtom3(:)=hAtom3(:)+cAtom(:)

		!letztes C-Atom ausgeben
		!Write(23,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
        koordH(k,:)=cAtom(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='C ' 
        k=k+1       
        !Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
        ! Rotation 
        
        koordH(k,:)=hAtom1(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='H ' 
        k=k+1           
		!Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
        koordH(k,:)=hAtom2(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='H ' 
        k=k+1     
		!Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom3(1),hAtom3(2),hAtom3(3)
        koordH(k,:)=hAtom3(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='H ' 
        k=k+1 
	END IF
END DO

IF(.NOT. HEnde) THEN  ! Anfügen des Letzten H-Atoms zum Absättigen
    c1Atom=koord(MAXATM-7,:)   ! C1-Atom 6 nach der Hexyl-Kette + 1
    c2Atom=koord(MAXATM-8,:)   ! C2-Atom 6 nach der Hexyl-Kette + 1
    hAtom1(:)=c2Atom(:) 
    sAtom=koord(MAXATM-6,:)    ! S-Atom an Potsiton 6 nach der Hexyl-Kette
    CALL setDistance(c1Atom,hAtom1,1.095) !Abstand anpassen
    Call setAngle(c1Atom,c2Atom,sAtom,hAtom1,(360.0-119.9)) ! 119.9Winkel in Grad
    koordH(k,:)= hAtom1
    atomsorteH(k)='H ' 
    k=k+1
END IF

IF(.NOT. HStart) THEN  ! Anfügen des Ersten H-Atoms zum Absättigen
    c1Atom=koord(1,:)   ! C1-Atom 
    c2Atom=koord(2,:)   ! C2-Atom
    hAtom1(:)=c2Atom    ! Start-H-Position
    sAtom=koord(6,:)    ! S-Atom an Potsiton 6
    CALL setDistance(c1Atom,hAtom1,1.095) !Abstand anpassen
    Call setAngle(c1Atom,c2Atom,sAtom,hAtom1,(360.0-119.9)) ! 119.9 Winkel in Grad
    ! Tauschen H-Atom an den Arrayanfang
    atomsorteH(2:anzahlneu)=atomsorteH(1:anzahlneu-1)
    atomsorteH(1)='H '
    koordH(2:anzahlneu,:)=koordH(1:anzahlneu-1,:)
    koordH(1,:)=hAtom1
    k=k+1
END IF

DEALLocate(atomsortedummy)

IF(DEBUG) THEN
	WRITE(*,*) 'Koordinaten nach h_P3HT ' 
	DO i=1,size(atomsorteH),1
        	write(*,*)  atomsorteH(i),(koordH(i,j),j=1,3)
	END DO
END IF
END SUBROUTINE h_P3HT


!! creates P3MT / adds methy H atoms
SUBROUTINE h_P3MT(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel)
IMPLICIT NONE
Integer, INTENT (IN):: MAXATM
Real, Dimension(MAXATM,3), INTENT (IN) :: koord
Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
Character(5), ALLOCATABLE, Dimension(:), INTENT (OUT) :: atomsorteH
Character(5), ALLOCATABLE, Dimension(:)::atomsortedummy
Character(500), INTENT (IN) :: ziel
Integer ::i,j,k,ind,anzahl,anzahlneu,anzahlringe,ind5RingEnde,ierror1,cAtomRingInd,sAtomInd,zaehler
Character(3) :: dummy
Real, Dimension(3) ::X,cAtom,cAtomRing,hAtom1,hAtom2,hAtom3,vorigesC,naechstesC,kreuzprod,c1Atom,c2Atom,sAtom
Integer, Dimension(MAXATM) :: ind5Ring
LOGICAL::HStart,HEnde
Real::Rotwinkel

HStart=.false.
HEnde=.false.

! Bestimmung ob erstes und letztes H-Atom existiert
IF( TRIM(adjustl(atomsorte(1))) == 'HC' .AND. TRIM(adjustl(atomsorte(2))) == 'CA' ) HStart=.true.
IF( TRIM(adjustl(atomsorte(MAXATM-7))) == 'HC' ) HEnde=.true.

!Indexarray fur 5Ringe und Wasserstoffe erzeugen
anzahl=MAXATM
i=1
ind=1
ind5Ring=0

DO WHILE(ind <= anzahl)
	ind5Ring(i) = ind
	IF(index(atomsorte(ind),'S') >0) THEN !Hexylgruppe ueberspringen
		ind = ind +6 !jetzt: ind letztes Atom der Hexylgruppe
	END IF
	i=i+1
	ind=ind+1
END DO
ind5RingEnde = i-1 !letzter belegter Platz

IF( HStart .AND. HEnde ) THEN
    anzahlringe = (ind5RingEnde -2)/6              !anzahlringe = (ind5RingEnde -2)/6 !fuer ganze Kette 
    anzahlneu = anzahl - anzahlringe*5 + anzahlringe*3 !Gesamtanzahl Atome in neuer Datei, 5 C weg, 3 H dazu je Ring
else
    IF(HStart .OR. HEnde) THEN
        anzahlringe = (ind5RingEnde -1)/6  
        anzahlneu = anzahl - anzahlringe*5 + anzahlringe*3 + 1 !Gesamtanzahl Atome in neuer Datei, 5 C weg, 3 H dazu je Ring +1 fuer Anfangs/Abschluss H-Atom
    ELSE 
        anzahlringe = (ind5RingEnde)/6 
        anzahlneu = anzahl - anzahlringe*5 + anzahlringe*3 + 2 !Gesamtanzahl Atome in neuer Datei, 5 C weg, 3 H dazu je Ring +2 fuer Anfangs+Abschluss H-Atom
    END IF
END IF

!write(*,*) 'HStart=',HStart
!write(*,*) 'HEnde= ',HEnde
!WRITE(*,*) "anzahlringe",anzahlringe
!write(*,*) "anzahlneu",anzahlneu
!write(*,*) "ind5RingEnde",ind5RingEnde
ALLOCATE(koordH(anzahlneu,3))
ALLOCATE(atomsorteH(anzahlneu))
ALLocate(atomsortedummy(size(atomsorte)))

ind =1
i=1
DO WHILE(i <= MAXATM)
     dummy = adjustl(atomsorte(i))
     atomsorteH(ind) = dummy(1:1)//'  '
     atomsortedummy(i) = dummy(1:1)//'  '
     IF(index(atomsorte(i),'S') >0) THEN
        ind = ind+4                ! Ueberspringen um 4 Atome =1*C+3*H
        i=i+6
     END IF
    ind=ind+1
    i=i+1
END DO 

k=1       ! Index fuer das Speichern der Koordinaten in KoordH
!neue Koordinaten schreiben
DO i=1,ind5RingEnde,1
	!Write(23,'(2X,A2,1X,3F16.6)') atomsortedummy(ind5Ring(i)),koord(ind5Ring(i),1),koord(ind5Ring(i),2),koord(ind5Ring(i),3)  
    koordH(k,:)=koord(ind5Ring(i),:) !Uebergabe der Koordinaten an koordH
	k=k+1

    IF(index(atomsortedummy(ind5Ring(i)),'S') >0) THEN !H-Atome an Hexylgruppe anfangen
		DO ind=1,3,1 !Atompositionen speichern
			cAtomRing(ind) = koord(ind5Ring(i)-2,ind) !C-Atom an das Hexylgruppe gebunden ist
			cAtom(ind) = koord(ind5Ring(i)+1,ind)  !1. Atom Hexylgruppe
			hAtom1(ind) = koord(ind5Ring(i)+2,ind) !2. Atom der Hexylgruppe liefert ca. Position 1. H
		END DO
		!3 H-Atome an das C-Atom anfügen
		CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen

		DO ind =1,3,1
		hAtom2(ind) = hAtom1(ind) !schon mal die ungefaehre Position
		hAtom3(ind) = hAtom1(ind)
		END DO

		CALL rotate(cAtomRing,cAtom,hAtom2,120.) !in korrekte Position bringen
		CALL rotate(cAtomRing,cAtom,hAtom3,240.)

                !Start Rotation um Dih zum H-Atom anzupassen
                call set_dih_by_rotation_around_axis(hAtom1,cAtom,cAtomRing,koord(ind5Ring(i)-3,:),11.3,Rotwinkel)
                X=hAtom2(:)-cAtom(:)
                call Rotationsmatrix_alpha(X,Rotwinkel,cAtomRing-cAtom,.true.)
                hAtom2(:)=X(:)+cAtom(:)

                X=hAtom3(:)-cAtom(:)
                call Rotationsmatrix_alpha(X,Rotwinkel,cAtomRing-cAtom,.true.)
                hAtom3(:)=X(:)+cAtom(:)

                !Ende Extra Rotation

		!Methylgruppe ausgeben
		!Write(23,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
        koordH(k,:)=cAtom(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='C ' 
        k=k+1
		!Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
        koordH(k,:)=hAtom1(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='H ' 
        k=k+1       
		!Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
        koordH(k,:)=hAtom2(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='H ' 
        k=k+1         
		!Write(23,'(2X,A2,1X,3F16.6)') 'H ',hAtom3(1),hAtom3(2),hAtom3(3)
        koordH(k,:)=hAtom3(:) !Uebergabe der Koordinaten an koordH
        atomsorteH(k)='H ' 
        k=k+1 

    END IF
END DO

CLOSE(23)

IF(.NOT. HEnde) THEN  ! Anfügen des Letzten H-Atoms zum Absättigen
    c1Atom=koord(MAXATM-7,:)   ! C1-Atom 6 nach der Hexyl-Kette + 1
    c2Atom=koord(MAXATM-8,:)   ! C2-Atom 6 nach der Hexyl-Kette + 1
    hAtom1(:)=c2Atom(:) 
    sAtom=koord(MAXATM-6,:)    ! S-Atom an Potsiton 6 nach der Hexyl-Kette
    CALL setDistance(c1Atom,hAtom1,1.095) !Abstand anpassen
    Call setAngle(c1Atom,c2Atom,sAtom,hAtom1,(360.0-119.9)) ! 119.9Winkel in Grad
    koordH(k,:)= hAtom1
    atomsorteH(k)='H ' 
    k=k+1
END IF

IF(.NOT. HStart) THEN  ! Anfügen des Ersten H-Atoms zum Absättigen
    c1Atom=koord(1,:)   ! C1-Atom 
    c2Atom=koord(2,:)   ! C2-Atom
    hAtom1(:)=c2Atom    ! Start-H-Position
    sAtom=koord(6,:)    ! S-Atom an Potsiton 6
    CALL setDistance(c1Atom,hAtom1,1.095) !Abstand anpassen
    Call setAngle(c1Atom,c2Atom,sAtom,hAtom1,(360.0-119.9)) ! 119.9 Winkel in Grad
    ! Tauschen H-Atom an den Arrayanfang
    atomsorteH(2:anzahlneu)=atomsorteH(1:anzahlneu-1)
    atomsorteH(1)='H '
    koordH(2:anzahlneu,:)=koordH(1:anzahlneu-1,:)
    koordH(1,:)=hAtom1
    k=k+1
END IF

!WRITE(*,*) ' P3MT ENDE KOORDINATEN:'
!DO i=1,size(atomsorteH),1
!    write(*,*)  atomsorteH(i),(koordH(i,j),j=1,3)
!END DO

END SUBROUTINE h_P3MT

SUBROUTINE h_diPBI(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel) ! Schneidet die Ketten ab!
        !Subroutine to cut the chains and substitute by methyl groups (CH3) to a DiPBI-molecule 
        !obtained by a GROMACS-file.
        IMPLICIT NONE
        Integer, INTENT(IN) :: MAXATM
        Real, Dimension(MAXATM,3), INTENT(IN) :: koord
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
        Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
        Character(500) :: ziel
        Real, Dimension(3) ::cAtom,NAtom,X,hAtom1,hAtom2,hAtom3,kreuzprod
        INTEGER::i,j,p,anzahl_methylgruppen=4
        ! Mindex C1-N2-C5-O3-H1-H2-H3
        Integer, Dimension(4,7)::Mindex
        Real::Rotwinkel
        ALLOCATE(koordH(86,3))
        ALLOCATE(atomsorteH(86))

        Mindex=transpose(reshape( &
    &(/ 1, 2, 5, 3, 75, 76, 77, 37, 35, 32, 34, 78, 79, 80, 38, 39, 42, 40, 81, 82, 83, 74, 72, 69,71, 84, 85, 86/), (/ 7, 4/)))
        
        DO i=1,74
                koordH(i,:)=koord(i,:)           !Koordinaten an neuen Array übergeben
                atomsorteH(i)=atomsorte(i)
        END DO
        ! Addition von zusaetzlichen H atomen an C mit Abschneiden der Ketteh
        
   DO p=1,anzahl_methylgruppen
        
        cAtom=koord(Mindex(p,1),:)
        NAtom=koord(Mindex(p,2),:)
        CALL kreuzprodukt(CAtom - NAtom, cAtom-koord(Mindex(p,3),:),kreuzprod)
        hAtom1 = NAtom
        CALL setDistance(cAtom,hAtom1,1.095)
        CALL rotate(cAtom,(cAtom+kreuzprod),hAtom1,112.5)
        hAtom2 = hAtom1
        hAtom3 = hAtom1
        CALL rotate(NAtom,cAtom,hAtom2,120.) ! Rotationen von hAtom2
        CALL rotate(NAtom,cAtom,hAtom3,240.) ! Rotation von hAtom3
    
       !Start Rotation um Dih zum H-Atom anzupassen
        call set_dih_by_rotation_around_axis(hAtom1,cAtom,NAtom,koord(Mindex(p,4),:),7.8,Rotwinkel)
        X=hAtom2(:)-cAtom(:)
        call Rotationsmatrix_alpha(X,Rotwinkel,NAtom-cAtom,.true.)
        hAtom2(:)=X(:)+cAtom(:)
 
        X=hAtom3(:)-cAtom(:)
        call Rotationsmatrix_alpha(X,Rotwinkel,NAtom-cAtom,.true.)
        hAtom3(:)=X(:)+cAtom(:)
        !Ende Extra Rotation
    
        koordH(Mindex(p,5),:)=hAtom1(:)
        koordH(Mindex(p,6),:)=hAtom2(:)
        koordH(Mindex(p,7),:)=hAtom3(:)
        atomsorteH(Mindex(p,5))='H '
        atomsorteH(Mindex(p,6))='H '
        atomsorteH(Mindex(p,7))='H '

   END DO ! i-Schleife Methylgruppen
   !write(*,*) 'Ende h_diPBI new'
END SUBROUTINE h_diPBI



SUBROUTINE h_PPDI_propyl(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel) ! Schneidet die Ketten ab!
        !Subroutine to cut the chains and substitute by propyl groups (CH3-CH-CH3)-N to a PPDI-molecule 
        IMPLICIT NONE
        Integer, INTENT(IN) :: MAXATM
        Real, Dimension(MAXATM,3), INTENT(IN) :: koord
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
        Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
        Character(500) :: ziel
        Character(500) :: molA
        CHARACTER(1)   ::dummy
        Real, Dimension(3) ::cAtom,NAtom,X,hAtom1,hAtom2,hAtom3,kreuzprod,OAtom,hAtomA,hAtomB
        INTEGER::i,j,p,anzahl_methylgruppen=4
        ! Mindex C5-C6-N7-O8-H1-H2-H3
        Integer, Dimension(4,7)::Mindex
        Integer, Dimension(6,2)::Cindex
        Real::Rotwinkel,dist_sum_1, dist_sum_2
        ALLOCATE(koordH(96,3))
        ALLOCATE(atomsorteH(96))
        
        ! Die ersten 4 Zahlen geben jeweils den alten Index an (C5-C6-N7-O9) -- dann kommen die neuen Indizes (H77-H78-H79)
        ! O kann O9 oder O20 am N7 werden oder (O81 oder O74 am N81)
        ! 4 Methylgruppen werden gesetzt.
        Mindex=transpose(reshape( &
    &(/ 5,6,7,9,78,79,80, 94,6,7,9,84,85,86,  84,83,82,81,88,89,90, 89,83,82,81,94,95,96 /), (/ 7, 4/)))
        
    j=0
    DO i=7,82  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
            j=j+1
            koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
            dummy = TRIM(adjustl(atomsorte(i)))
            atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
            !WRITE(*,*) j,atomsorteH(j),koordH(j,:)
    END DO 
            
            
    !!! Uebergabe des C-Atome aus der Seitenkette an die Position im Array 1=alter Index ,2= neuer Index Anzahl der C-Atome durchlaufen als i=1,6 
    !! C5 aus koord wird zum ersten neuen C in koordH C77 daher 5,77 ... dann wird Platz fuer 3 H-Atome gelassen (neuer Index 78 79 80) und C6 wird 
    Cindex=transpose(reshape( (/ 5,77, 6,81, 94,83, 84,87, 83,91, 89,93 /), (/ 2, 6/)))
    Do i=1,6
        !WRITE(*,*) Cindex(i,1),Cindex(i,2)
        koordH(Cindex(i,2),:)=koord(Cindex(i,1),:)
        atomsorteH(Cindex(i,2))='C'
    END DO
    
    !H-Start
    cAtom=koordH(81,:)   !  C-Atom am Randstart
    NAtom=koordH(1,:)  
    hAtomA(:)=koordH(83,:) 
    ! Test welches O als Referenz fuer das H-ausgewaehlt werden soll, in dem der groessere Abstand der beiden C-Atome an denen CH3- gebaut werden soll bestimmt wird
    ! Test O9 und O20 in gro Nummerierung
    dist_sum_1=norm(koord(9,:)-koord(5,:))+norm(koord(9,:)-koord(94,:))
    dist_sum_2=norm(koord(20,:)-koord(5,:))+norm(koord(20,:)-koord(94,:))
    IF (dist_sum_1 > dist_sum_2) THEN
        OAtom=koord(9,:)
    ELSE 
        OAtom=koord(20,:)
    END IF
    CALL setDistance(cAtom,hAtomA,1.095) 
    call set_dih_by_rotation_around_axis(hAtomA,CAtom,NAtom,OAtom,2.9,Rotwinkel)
    koordH(82,:)= hAtomA
    atomsorteH(82)='H' 

    ! H-Ende
    cAtom=koordH(91,:)  ! C-Atom am Randende
    NAtom=koordH(76,:)  
    hAtomB(:)=koordH(93,:)
    ! Test welches O als Referenz fuer das H-ausgewaehlt werden soll, in dem der groessere Abstand der beiden C-Atome an denen CH3- gebaut werden soll bestimmt wird
    ! Test O81 und O74in gro Nummerierung
    dist_sum_1=norm(koord(81,:)-koord(84,:))+norm(koord(81,:)-koord(89,:))
    dist_sum_2=norm(koord(74,:)-koord(84,:))+norm(koord(74,:)-koord(89,:))
    IF (dist_sum_1 > dist_sum_2) THEN
        OAtom=koord(81,:)
    ELSE 
        OAtom=koord(74,:)
    END IF
    CALL setDistance(cAtom,hAtomB,1.095)
    call set_dih_by_rotation_around_axis(hAtomB,CAtom,NAtom,OAtom,2.9,Rotwinkel)
    koordH(92,:)= hAtomB
    atomsorteH(92)='H' 
    
   ! Schleife zum Anbau der 4 Methylendgruppen 
   DO p=1,anzahl_methylgruppen     !! Vorsicht die Benennung der Atome stimmt ab hier nicht mehr mit den Benennungen im PPDI ueberein, da die DIPBI Routine verwendet wurde 
        cAtom=koord(Mindex(p,1),:)
        NAtom=koord(Mindex(p,2),:)
        CALL kreuzprodukt(CAtom - NAtom, cAtom-koord(Mindex(p,3),:),kreuzprod)
        hAtom1 = NAtom
        CALL setDistance(cAtom,hAtom1,1.095)
        CALL rotate(cAtom,(cAtom+kreuzprod),hAtom1,112.5)
        hAtom2 = hAtom1
        hAtom3 = hAtom1
        CALL rotate(NAtom,cAtom,hAtom2,120.) ! Rotationen von hAtom2
        CALL rotate(NAtom,cAtom,hAtom3,240.) ! Rotation von hAtom3
    
       !Start Rotation um Dih zum H-Atom anzupassen
        if (p==1 .or. p==2) THEN
            call set_dih_by_rotation_around_axis(hAtom1,cAtom,NAtom,hAtomA,60.0,Rotwinkel)
        ELSE IF (p==3 .or. p==4) THEN
            call set_dih_by_rotation_around_axis(hAtom1,cAtom,NAtom,hAtomB,60.0,Rotwinkel)
        END IF
        !call set_dih_by_rotation_around_axis(hAtom1,cAtom,NAtom,koord(Mindex(p,4),:),7.8,Rotwinkel)
        X=hAtom2(:)-cAtom(:)
        call Rotationsmatrix_alpha(X,Rotwinkel,NAtom-cAtom,.true.)
        hAtom2(:)=X(:)+cAtom(:)
 
        X=hAtom3(:)-cAtom(:)
        call Rotationsmatrix_alpha(X,Rotwinkel,NAtom-cAtom,.true.)
        hAtom3(:)=X(:)+cAtom(:)
        !Ende Extra Rotation
    
        koordH(Mindex(p,5),:)=hAtom1(:)
        koordH(Mindex(p,6),:)=hAtom2(:)
        koordH(Mindex(p,7),:)=hAtom3(:)
        atomsorteH(Mindex(p,5))='H'
        atomsorteH(Mindex(p,6))='H'
        atomsorteH(Mindex(p,7))='H'

   END DO ! i-Schleife Methylgruppen

END SUBROUTINE h_PPDI_propyl





SUBROUTINE h_PPDI(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel) ! Schneidet die Ketten ab!
        !!Subroutine to cut the chains and substitute by methyl groups (CH3) to a PPI-molecule 
        !!obtained by a GROMACS-file.
        !! H-Atoms will be added with at the N-Atom with reference to the C-Atom on the oposite ring position (para)position.
        IMPLICIT NONE
        Integer, INTENT(IN) :: MAXATM
        Real, Dimension(MAXATM,3), INTENT(IN) :: koord
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
        Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
        Character(500) :: ziel
        Character(500) :: molA
        CHARACTER(1)::dummy
        Real, Dimension(3) ::c2Atom,c3Atom,NAtom,hAtom1
        INTEGER::i,j,p,anzahl_gruppen=2
        ALLOCATE(koordH(78,3))
        ALLOCATE(atomsorteH(78))
        WRITE(*,*) 'Start h_PPDI' 

        j=0
        DO i=7,82  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                j=j+1
                koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                dummy = TRIM(adjustl(atomsorte(i)))
                atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
                !WRITE(*,*) j,atomsorteH(j),koordH(j,:)
        END DO
        ! Addition von zusaetzlichen H atomen an C mit Abschneiden der Ketteh

        ! H-Start
        NAtom=koordH(1,:)   !  N-Atom am Randstart
        c2Atom=koordH(9,:)   ! Referenz C -Atom, welches im eingebauten Ring in der gegenueberliegenden Position (Parastellung) liegt
        hAtom1(:)=c2Atom(:) 
        C3Atom=koord(2,:)    ! S-Atom an Potsiton 6 nach der Hexyl-Kette
        CALL setDistance(NAtom,hAtom1,-1.095) !Abstand anpassen -- das negative Vorzeichen wandelt die Anbaurichtung um
        koordH(77,:)= hAtom1
        atomsorteH(77)='H' 
    
        ! H-Ende
        NAtom=koordH(76,:)   !  N-Atom am Randende
        c2Atom=koordH(63,:)  ! Referenz C -Atom, welches im eingebauten Ring in der gegenueberliegenden Position (Parastellung) liegt
        hAtom1(:)=c2Atom(:) 
        C3Atom=koord(75,:)    ! C-Atom seitlich vom anzudockenden N im Ring (in Richtung O) beliebig eines ist auswaehlbar
        CALL setDistance(NAtom,hAtom1,-1.095) !Abstand anpassen -- das negative Vorzeichen wandelt die Anbaurichtung um 
        koordH(78,:)= hAtom1
        atomsorteH(78)='H' 
 
END SUBROUTINE h_PPDI

SUBROUTINE h_PPDI_kette(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel) 
!! Subroutine to add H-Atoms to the chains so you end up with two CH3-CH2-CH2-CH3-CH-CH2-CH2-CH2-CH2-CH3 chains
!! First copies the Data for the N-Ringsystem-N, then adds the two CH-groups and the big loop adds the CH2-...-CH3 chains.
!! Addaption von h_DIPBI (uses definition with S (even if there is no S in the molecule)
!! Im Optput werden die Ketten nach dem N-Ringsystem-N geschrieben. -> Die Reihenfolge wird nicht wie im *.gro beibehalten!
IMPLICIT NONE
Integer, INTENT(IN) :: MAXATM
Real, Dimension(MAXATM,3), INTENT(IN) :: koord
Character(5), Dimension(MAXATM), INTENT(INOUT) :: atomsorte
REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
Character(500), INTENT(IN) :: ziel
CHARACTER(500) :: molA
CHARACTER(3) ::dummy
Integer :: i,j,ind,anzahl,anzahlneu,anzahlringe,ind5RingEnde,cAtomRingInd,sAtomInd,zaehler
 Real, Dimension(3) :: cAtom,cAtomRing,hAtom1,hAtom2,hAtom3,vorigesC,naechstesC,kreuzprod,hAtomold,hAtomA,hAtomB,OAtom,NAtom
 Integer, Dimension(20) :: ind5Ring ! 48?
 Integer, Dimension(20) ::ind_C_Ketten
 LOGICAL:: Cvorher, erste_vorher
 Integer:: ketten, erste_Kette, ketten_2
 Integer, Dimension(4):: ind_erste
 Real::Rotwinkel,dist_sum_1,dist_sum_2

!WRITE(*,*) 'h_PPDI_kette start'
!write(*,*) 'MAXATM', MAXATM
anzahl=MAXATM
anzahlneu = MAXATM + 46 !Gesamtanzahl Atome in neuer Datei

allocate(atomsorteH(anzahlneu))
allocate(koordH(anzahlneu,3))
! Übergabe der ersten Atome des DIPBI ohne Seitenketten - jeweils vom N 
j=0
DO i=7,82  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
    j=j+1
    koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
    dummy = TRIM(adjustl(atomsorte(i)))
    atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
    !WRITE(*,*) j,atomsorteH(j),koordH(j,:)
END DO 


! Indizes der 48 C Atome in den 4 Seitenketten
ind_C_Ketten(1)=1   ! C
ind=1
do i=5,1,-1
    ind_C_Ketten(ind)=i ! von 1 bis 5 ; Reihenfolge rueckwaerts
    ind=ind+1
END DO

! 2. Kettenanteil
do i=94,98,1
    ind_C_Ketten(ind)=i ! von 94-98
    ind=ind+1
END DO

do i=89,93,1
    ind_C_Ketten(ind)=i ! von 89-83
    ind=ind+1
END DO

do i=84,88,1
    ind_C_Ketten(ind)=i ! von 84-88
    ind=ind+1
END DO


do i=2,5,1
    atomsorte(i)="CH2"
END DO 
do i=84,97,1
    atomsorte(i)="CH2"
END DO 
atomsorte(1)="CH3"
atomsorte(98)="CH3"
atomsorte(93)="CH3"
atomsorte(88)="CH3"

    !H-Start
    cAtom=koord(6,:)   !  C-Atom am Randstart
    NAtom=koord(7,:)  
    hAtomA(:)=koord(5,:) 
    ! Test welches O als Referenz fuer das H-ausgewaehlt werden soll, in dem der groessere Abstand der beiden C-Atome an denen CH3- gebaut werden soll bestimmt wird
    ! Test O9 und O20 in gro Nummerierung
    dist_sum_1=norm(koord(9,:)-koord(5,:))+norm(koord(9,:)-koord(94,:))
    dist_sum_2=norm(koord(20,:)-koord(5,:))+norm(koord(20,:)-koord(94,:))
    IF (dist_sum_1 > dist_sum_2) THEN
        OAtom=koord(9,:)
    ELSE 
        OAtom=koord(20,:)
    END IF
    CALL setDistance(cAtom,hAtomA,1.095) 
    call set_dih_by_rotation_around_axis(hAtomA,CAtom,NAtom,OAtom,2.9,Rotwinkel)
    ! CH und H ueberfeveb
    koordH(77,:)= koord(6,:)
    atomsorteH(77)='C' 
    koordH(78,:)= hAtomA
    atomsorteH(78)='H' 

    ! H-Ende
    cAtom=koord(83,:)  ! C-Atom am Randende
    NAtom=koord(82,:)  
    hAtomB(:)=koord(84,:)
    ! Test welches O als Referenz fuer das H-ausgewaehlt werden soll, in dem der groessere Abstand der beiden C-Atome an denen CH3- gebaut werden soll bestimmt wird
    ! Test O81 und O74 in gro Nummerierung
    dist_sum_1=norm(koord(81,:)-koord(84,:))+norm(koord(81,:)-koord(89,:))
    dist_sum_2=norm(koord(74,:)-koord(84,:))+norm(koord(74,:)-koord(89,:))
    IF (dist_sum_1 > dist_sum_2) THEN
        OAtom=koord(81,:)
    ELSE 
        OAtom=koord(74,:)
    END IF
    CALL setDistance(cAtom,hAtomB,1.095)
    call set_dih_by_rotation_around_axis(hAtomB,CAtom,NAtom,OAtom,2.9,Rotwinkel)
    koordH(79,:)= koord(83,:)
    atomsorteH(79)='C'   
    koordH(80,:)= hAtomB
    atomsorteH(80)='H' 

    !DO j=1, size(atomsorteH),1
     !   WRITE(*,*) j,atomsorteH(j),koordH(j,:)
    !END DO
    
    !Do i=1,SiZE(ind_C_Ketten)
    !    WRITE(*,*) ind_C_Ketten(i),atomsorte(ind_C_Ketten(i)),koord(ind_C_Ketten(i),:) 
    !END DO

cvorher=.false.  !Logikvariable für das erste CH2 in jeder Kette
ketten=0 		 !Abzählen der 4 C12-Ketten
erste_vorher= .false. 	!Logikvariable für das zweite CH2 jeder Kette
ketten_2=1
j=81  ! Index zum hochzaehlen der Atome
!neue Koordinaten schreiben 
!Die entsprechenden H werden jeweils nach dem C-Atom angefuegt
DO i=1,size(ind_C_Ketten),1
    If (TRIM(atomsorte(ind_C_Ketten(i)))=='CH3' .OR. (TRIM(atomsorte(ind_C_Ketten(i)))=='CH2')) Then
        koordH(j,:)=koord(ind_C_Ketten(i),:)    ! Ausnahme der 4 schon geschriebenen ersten C-Atome der Kette
        atomsorteH(j)='C'        ! Zuerst jeweils das C in der Seitenkette an die Stelle j schreiben
        j=j+1
        !WRITE(*,*) j,ind_C_Ketten(i),atomsorte(ind_C_Ketten(i)),koord(ind_C_Ketten(i),:) 
    End If

	IF(atomsorte(ind_C_Ketten(i)) == 'CH2') THEN !H-Atome an Hexylgruppe anfügen
        !Falls es das erste CH2 der Kette ist:
    	IF(cvorher .eqv. .false.) Then  ! Beginn am Kettenanfang 
            ind_erste(ketten+1)=ind_C_Ketten(i)
            if (ketten==0) Then
              erste_Kette=ind_C_Ketten(i)
            End if

            !1. C-Atom mit H bestuecken
			DO ind=1,3,1 !Koordinaten merken
				cAtom(ind) = koord(ind_C_Ketten(i),ind)
				naechstesC(ind) = koord(ind_C_Ketten(i+1),ind)
                if ((ketten == 0) .or. (ketten == 1)) Then
					vorigesC(ind) = koord(6,ind)   ! Auswahl von C6 
                ELSE IF(ketten == 2 .or. ketten ==3 ) THEN
                    vorigesC(ind) = koord(83,ind)   ! Auswahl von C83
                End if

			END DO
            ketten=ketten+1
            !Write(*,*) 'Ketten', ketten
			if(ketten==4) then
             			erste_vorher = .true.
						cvorher = .true.
            End if
        
            CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
            hAtom1 = cAtom + kreuzprod ! Vektor senkrecht zur Ebene der 3 C-Atome
            hAtom2 = cAtom - kreuzprod
        
            CALL rotate(naechstesC,vorigesC,hAtom1,32.)  !Winkel anpassen
            CALL rotate(naechstesC,vorigesC,hAtom2,-32.) ! Rotation nach außen um die Verbinungdachse der Cs
        
            CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen
            CALL setDistance(cAtom,hAtom2,1.095) !Abstand anpassen
                    
            !Write(2,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
            koordH(j,:)=hAtom1(:)    
            atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
            j=j+1
            koordH(j,:)=hAtom2(:)    
            atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
            j=j+1
            cvorher=.true.   

        Else 
            if (erste_vorher .eqv. .true.) Then
                !2. C-Atom mit H bestuecken - hier muess wieder auf das 1. C zugegriffen werden 
                !Daher wird eine Adressierung mit ketten_2 eingefuehrt
                DO ind=1,3,1 !Koordinaten merken
                    cAtom(ind) = koord(ind_C_Ketten(i),ind)
                    vorigesC(ind) = koord(ind_erste(ketten_2+1),ind) !ketten_2=2 muss erst bei der zweiten Kette loslaufen
                    naechstesC(ind) = koord(ind_C_Ketten(i+1),ind)   ! daher ist ketten_2=1 zu Beginn 
                END DO
                          
                erste_vorher = .false.
                ketten_2=ketten_2+1

                
                CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
    
                hAtom1 = cAtom + kreuzprod
                hAtom2 = cAtom - kreuzprod
 
                CALL rotate(naechstesC,vorigesC,hAtom1,32.) !Winkel anpassen
                CALL rotate(naechstesC,vorigesC,hAtom2,-32.)
                
                CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen
                CALL setDistance(cAtom,hAtom2,1.095) !Abstand anpassen
                   
                !Write(2,'(2X,A2,1X,3F16.6)') ' C ',cAtom(1),cAtom(2),cAtom(3)
                !Write(2,'(2X,A2,1X,3F16.6)') ' H ',hAtom1(1),hAtom1(2),hAtom1(3)
                !Write(2,'(2X,A2,1X,3F16.6)') ' H ',hAtom2(1),hAtom2(2),hAtom2(3)
                    koordH(j,:)=hAtom1(:)    
                    atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                    j=j+1
                    koordH(j,:)=hAtom2(:)    
                    atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                    j=j+1       
            
            Else    ! Durchlaufen vom 3. C-Atom der Kette bis zum vorletzen C-Atom.
                    !1. C-Atom mit H bestuecken
                DO ind=1,3,1 !Koordinaten merken
                    cAtom(ind) = koord(ind_C_Ketten(i),ind)
                    vorigesC(ind) = koord(ind_C_Ketten(i-1),ind)
                    naechstesC(ind) = koord(ind_C_Ketten(i+1),ind)
                END DO
                        
                CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
                hAtom1 = cAtom + kreuzprod
                hAtom2 = cAtom - kreuzprod
                
                CALL rotate(naechstesC,vorigesC,hAtom1,32.) !Winkel anpassen
                CALL rotate(naechstesC,vorigesC,hAtom2,-32.)
        
                CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen
                CALL setDistance(cAtom,hAtom2,1.095) !Abstand anpassen
                hAtomOld=hAtom1
                
                !Write(2,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
                !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
                !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
                koordH(j,:)=hAtom1(:)    
                atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                j=j+1
                koordH(j,:)=hAtom2(:)    
                atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                j=j+1
                        
                
            END IF!Write(*,*) 'CH2 Ende'
        End if
	Else
        If(atomsorte(ind_C_Ketten(i))=='CH3') THEN 
            cAtomRingInd=  ind_C_Ketten(i - 2)
            sAtomInd= ind_C_Ketten(i)
            !letztes C-Atom mit 3 H-Atomen bestücken
            DO ind=1,3,1 !Koordinaten merken
                cAtom(ind) = koord(sAtomInd,ind)
                vorigesC(ind) = koord(cAtomRingInd+1,ind)
                naechstesC(ind) = koord(sAtomInd+1,ind) !vorvoriges C-atom
            END DO		
            
            CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
    
            hAtom1 = vorigesC
            CALL setDistance(cAtom,hAtom1,1.095)
            CALL rotate(cAtom,cAtom+kreuzprod,hAtom1,112.)
    
            hAtom2 = hAtom1  !schon mal die ungefähre Position
            hAtom3 = hAtom1
    
            CALL rotate(vorigesC,cAtom,hAtom2,120.) !in korrekte Position bringen
            CALL rotate(vorigesC,cAtom,hAtom3,240.)
            
           !Rotation um C-C Achse sodass der Diederwinkel H-C-C-H ~60° wird
           ! (Korrektur zur passenden Newman-Projektion)
           call set_dih_by_rotation_around_axis(hAtom1,cAtom,vorigesC,hAtomold,60.,Rotwinkel)
           hAtom2=hAtom2-cAtom(:)
           call Rotationsmatrix_alpha(hAtom2,Rotwinkel,vorigesC(:)-cAtom(:),.true.)
           hAtom2(:)=hAtom2(:)+cAtom(:)
          
           hAtom3=hAtom3-cAtom(:)
           call Rotationsmatrix_alpha(hAtom3,Rotwinkel,vorigesC(:)-cAtom(:),.true.)
           hAtom3(:)=hAtom3(:)+cAtom(:)

            !letztes C-Atom ausgeben
            !Write(2,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom3(1),hAtom3(2),hAtom3(3)
            koordH(j,:)=hAtom1(:)    
            atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
            j=j+1
            koordH(j,:)=hAtom2(:)    
            atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
            j=j+1
            koordH(j,:)=hAtom3(:)    
            atomsorteH(j)='H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
            j=j+1
            
            cvorher=.false. !Nächste Kette
            erste_vorher=.true. !Nächste Kette
        END IF
    End IF
END DO
 !CLOSE(2)

! molA='.'
!CALL make_coord_to_xyz_in_dir(atomsorteH,koordH,size(atomsorteH),ziel,molA)
!WRITE(*,*) 'h_PPDI_Kette Ende'

END SUBROUTINE h_PPDI_kette



SUBROUTINE h_PBDT_TS1(atomsorte,koord,atomsorte_H_all,koord_H_all,MAXATM,ziel) ! Schneidet die Ketten ab!
    !! Subroutine to cut the chains and substitute by methyl groups (CH3) to a PPI-molecule 
    !! obtained by a GROMACS-file.
    !! H-Atoms will be added with at the N-Atom with reference to the C-Atom on the oposite ring position (para)position.
    IMPLICIT NONE
    Integer, INTENT(IN) :: MAXATM
    Real, Dimension(MAXATM,3), INTENT(IN) :: koord
    REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) :: koord_H_all
    Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorte_H_all
    Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
    Character(500) :: ziel
    Character(500) :: molA
    CHARACTER(1)::dummy
    Real, Dimension(3) ::c1Atom,c2Atom,c3Atom,NAtom,hAtom1,hAtom2,SAtom,kreuzprod
    INTEGER::i,j,p,i_start=0,Count_atoms_all=1,anzahl_gruppen=2
    INTEGER::N_CarboxyS2,N_PBDT,N_HC,N_gesamt
    REAL, ALLOCATABLE, DIMENSION(:,:)       :: koordH,koordCarboxyS2,koordPBDT
    Character(5),ALLOCATABLE,  Dimension(:) :: atomsorteH,atomsorteCarboxyS2,atomsortePBDT
    LOGICAL:: poly1=.false.,segment8=.false.,HStart=.true.,DEBUG=.false.
    LOGICAL, ALLOCATABLE, DIMENSION(:)::PBDT_chain
    Integer, Dimension(15,2)::Aindex
    Integer, Dimension(33,2)::Change_index
    REAL::Rotwinkel
    HStart=.true.
    poly1=.false.
    segment8=.false.
    
    !! Arrays as temporarry data storage, to send coordinates to subroutines h_CarboxS2_propyl, h_PBDT_CH3
    ALLOCATE(koordCarboxyS2(15,3))
    ALLOCATE(atomsorteCarboxyS2(15))
    
    ALLOCATE(koordPBDT(32,3))
    ALLOCATE(atomsortePBDT(32))
    
    IF (DEBUG) THEN
        WRITE(*,*) 'Start h_PBDT_TS1' 
        DO j=1,MAXATM,1
            WRITE(*,*)  atomsorte(j),koord(j,:)
        END DO 
    END IF
    
    
    ! Zaehle die Anzahl der Elemente im Molekuel N_CarboxyS2 ueber Anzahl der F-Atome und N_PBDT=(N_HC-Gruppen/6) uber die Anzahl der HC-Gruppen
    N_CarboxyS2=0
    N_PBDT=0
    DO i=1,MAXATM
        if (TRIM(adjustl(atomsorte(i)))=='F') N_CarboxyS2=N_CarboxyS2+1
        if (TRIM(adjustl(atomsorte(i)))=='HC') N_PBDT=N_PBDT+1
    END DO ! Zahle Anzahl
    N_PBDT =N_PBDT/6
    IF (DEBUG) WRITE(*,*) ' Die Anzahl der Carboxy-Gruppen ist: ',N_CarboxyS2
    IF (DEBUG) WRITE(*,*) ' Die Anzahl der N_PBDT-Gruppen ist: ',N_PBDT   
    N_gesamt=38*N_PBDT+22*N_CarboxyS2+2
    IF (DEBUG) WRITE(*,*) ' Die Gesamtanzahl nach dem Anfuegen der H-Atome betraegt',N_gesamt
    ALLOCATE(PBDT_chain(N_PBDT+N_CarboxyS2))
    ALLOCATE(koord_H_all(N_gesamt,3))
    ALLOCATE(atomsorte_H_all(N_gesamt))
    
    IF( (N_PBDT ==0) .AND. (N_CarboxyS2==0) ) THEN
        WRITE(*,*) 'Fehler: PBDT-TS_oder_CaboxyS2_wurden_nicht_richtig_eingelesen_pruefe_gro-Inputdatei. ENDE!'
        CALL EXIT (1)
    END IF
    PBDT_chain(:)=.false.
    
    ! Festlegen der Reihenfolge
    j=0
    N_HC=0
    DO i=1,MAXATM,1
        if (TRIM(adjustl(atomsorte(i)))=='F') THEN
            j=j+1
            PBDT_chain(j)=.false.
        END IF
        if (TRIM(adjustl(atomsorte(i)))=='HC') THEN
            N_HC=N_HC+1
            iF (N_HC == 6) THEN
                j=j+1
                PBDT_chain(j)=.true.
                N_HC=0
            END IF
        END IF
    END DO ! Zahle Anzahl
    
    IF(DEBUG) WRITE(*,*) 'Reihenfolge T= PBDT, F=CarboxylS2'
    IF(DEBUG) THEN
	    DO i=1,(N_PBDT+N_CarboxyS2)
		IF( PBDT_chain(i)) THEN
		    IF(i==1) THEN
		        WRITE(*,'(A5)',ADVANCE='NO') 'PBDT-'
		    ELSE IF (i==size(PBDT_chain)) THEN
		        WRITE(*,'(A5)') 'PBDT'
		    ELSE
		        WRITE(*,'(A5)',ADVANCE='NO') 'PBDT-'
		    END IF
		ELSE
		    IF(i==1) THEN
		        WRITE(*,'(A10)',ADVANCE='NO') 'CarboxyS2-'
		    ELSE IF (i==size(PBDT_chain)) THEN
		        WRITE(*,'(A10)') '-CarboxyS2'
		    ELSE 
		        WRITE(*,'(A10)',ADVANCE='NO') 'CarboxyS2-'
		    END IF
		END IF ! PBDT_chain
	    END DO 
     END IF    

    Count_atoms_all=1 ! Platzlassen fuer H1
    i_start=0
    ! Schleife uber alle Teilabschnitte
    DO p=1,(N_PBDT+N_CarboxyS2) 
        IF (PBDT_chain(p)) THEN !PBDT-Abschnitt?
            IF(DEBUG) THEN
                WRITE(*,*) 'Start PBDT_Abschnitt'
                j=0
                DO i=i_start+1,i_start+32  !Bsp 17,48 Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                    j=j+1
                    WRITE(*,*) j,atomsorte(i),koord(i,:)
                END DO
            END IF
            !!! Anpassung der Atomreihung in der gro-Datei, da diese nicht immer gleich ist, muss dies nun hier korrigiert werden (vor dem Anfuegen der H-Atome). 
            IF ( TRIM(adjustl(atomsorte(i_start+9))) == 'CH3' .AND. TRIM(adjustl(atomsorte(i_start+2))) == 'S'  ) THEN ! Befindet sich am ersten Teil des 8-mers
                !! Teil 1 (keine Anderung der Anordnung
                Change_index=transpose(reshape( (/ 1,1, 2,2, 3,3, 4,4, 5,5, 6,6, 7,7, 8,8, 9,9, 10,10, 11,11, 12,12, &
                &13,13, 14,14, 15,15, 16,16, 17,17, 18,18, 19,19, 20,20, 21,21, 22,22, 23,23, 24,24, 25,25, 26,26, 27,27, 28,28,&
                &29,29, 30,30, 31,31, 32,32,   32,32 /), (/ 2, 33/)))   
                IF (DEBUG) WRITE(*,*) 'Fall 1'
            ELSE IF ( TRIM(adjustl(atomsorte(i_start+15))) == 'CH3' .AND. TRIM(adjustl(atomsorte(i_start+2))) == 'S' &
                    & .AND. TRIM(adjustl(atomsorte(i_start+24))) /= 'CH3' ) THEN !! Teil 2
                !! Change Teil 1 to Teil 2 (letzte indizes sind doppelt, da H-Abschluss moeglich ist)
                Change_index=transpose(reshape( (/ 1,1, 2,2, 3,3, 4,19, 5,20, 6,28, 7,25, 8,26, 9,27, 10,23, 11,24, 12,21, &
                &13,22, 14,18, 15,17, 16,7, 17,8, 18,16, 19,13, 20,14, 21,15, 22,11, 23,12, 24,9, 25,10, 26,4, 27,5, 28,6, &
                &29,32, 30,31, 31,29, 32,30,   32,30 /), (/ 2, 33/)))       
                IF (DEBUG)WRITE(*,*) 'Fall 2'                
            ELSE IF ( TRIM(adjustl(atomsorte(i_start+12))) == 'CH3' ) THEN !! Teil 3,4
                !! Change Teil 1 to Teil 3
                Change_index=transpose(reshape( (/ 1,1, 2,6, 3,5, 4,7, 5,8, 6,9, 7,10, 8,11, 9,12, 10,13, 11,14, 12,15, 13,16, &
                & 14,17, 15,18, 16,19, 17,20, 18,28, 19,25, 20,26, 21,27, 22,23, 23,24, 24,21, 25,22, 26,4, 27,2, 28,3, 29,29, &
                &30,30, 31,31, 32,32,   32,32 /), (/ 2, 33/)))   
                IF (DEBUG) WRITE(*,*) 'Fall 3'    
            ELSE IF ( TRIM(adjustl(atomsorte(i_start+15))) == 'CH3' .AND.  TRIM(adjustl(atomsorte(i_start+6))) == 'S' ) THEN !! Teil 5,6
                Change_index=transpose(reshape( (/ 1,1, 2,6, 3,5, 4,7, 5,8, 6,16, 7,13, 8,14, 9,15, 10,11, 11,12, 12,9, 13,10, &
                & 14,17, 15,18, 16,19, 17,20, 18,28, 19,25, 20,26, 21,27, 22,23, 23,24, 24,21, 25,22, 26,4, 27,2, 28,3, 29,29, &
                &30,30, 31,31, 32,32,   32,32 /), (/ 2, 33/)))   
               IF (DEBUG) WRITE(*,*) 'Fall 4'            
            ELSE IF ( TRIM(adjustl(atomsorte(i_start+24))) == 'CH3' )  THEN !! Teil 7
                Change_index=transpose(reshape( (/ 1,1, 2,6, 3,5, 4,7, 5,8, 6,16, 7,13, 8,14, 9,15, 10,11, 11,12, 12,9, 13,10, &
                & 14,17, 15,18, 16,19, 17,20, 18,21, 19,22, 20,23, 21,24, 22,25, 23,26, 24,27, 25,28, 26,4, 27,2, 28,3, 29,29, &
                &30,30, 31,31, 32,32,   32,32 /), (/ 2, 33/)))             
            
            ELSE IF ( TRIM(adjustl(atomsorte(i_start+10))) == 'CH3' ) THEN !! Teil 8
                !! Change Teil 1 to Teil 8
                Change_index=transpose(reshape( (/ 1,1, 2,24, 3,23, 4,22, 5,25, 6,33, 7,30, 8,31, 9,32, 10,28, 11,29, 12,26, &
                &13,27, 14,21, 15,15, 16,5, 17,6, 18,7, 19,8, 20,9, 21,10, 22,11, 23,12, 24,13, 25,14, 26,4, 27,2, 28,3, 29,16, &
                &30,17, 31,19, 32,20, 33,18  /), (/ 2, 33/)))  !! mit H-Atom-Ende an letzter Position (33 bzw 18) 
                IF (DEBUG) WRITE(*,*) 'Fall 8'   
                segment8=.true.              
            END IF
        
        
            koordPBDT(:,:)=0
            IF (DEBUG) WRITE(*,*) 'i_start:',i_start,i_start+1,i_start+32
            IF (DEBUG) WRITE(*,*) 'Start vor h_PBDT_CH3'
            j=0
            DO i=i_start+1,i_start+32  !Bsp 17,48 Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                j=j+1
                koordPBDT(Change_index(j,1),:)=koord(i_start+Change_index(j,2),:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                dummy = TRIM(adjustl(atomsorte(i_start+Change_index(j,2))))
                atomsortePBDT(Change_index(j,1)) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
                IF (DEBUG) THEN 
        WRITE(*,*) Change_index(j,1),Change_index(j,2),j,atomsorte(i),koord(i,:),' changed_to ',atomsortePBDT(j),koordPBDT(j,:)
                END IF !DEBUG 
            END DO

        
            i_start=i_start+32
            If (allocated(koordH))  deallocate (koordH)
            If (allocated(atomsorteH))  deallocate (atomsorteH)
            CALL h_PBDT_CH3(atomsortePBDT,koordPBDT,atomsorteH,koordH,size(atomsortePBDT),ziel,DEBUG) 
            
            IF (DEBUG) WRITE(*,*) 'ENDE PBDT_Abschnitt'
          
        ELSE ! CarboxylS2-Gruppe   
            ! Test ob man sich am Anfang des 8meres befindet, welches mit F-Atom anfaengt
            if (TRIM(adjustl(atomsorte(i_start+1)))=='F') THEN ! start ist i_start=0 also erstes Atom im Array
                poly1=.true. ! Erstes Polymer!
                ! Reihenfolge muss angepasst werden!
                j=0
                Aindex=transpose(reshape( (/ 1,4, 2,3, 3,5, 4,6, 5,7, 6,8, 7,9, 8,10, &
                                                & 9,11, 10,12, 11,13, 12,2, 13,1, 15,15, 16,14  /), (/ 2, 15/)))

                Do i=i_start+1,i_start+15
                    IF (DEBUG) THEN
                        WRITE(*,*) Aindex(i,1),Aindex(i,2)
                    END IF
                    koordCarboxyS2(Aindex(i,2),:)=koord(Aindex(i,1),:)
                    atomsorteCarboxyS2(Aindex(i,2))=atomsorte(Aindex(i,1))
                END DO
                i_start=i_start+16
                
                IF (DEBUG) THEN
                    WRITE(*,*) 'rein Carb gro-format:'
                    DO i =1,15
                        WRITE(*,*) atomsorteCarboxyS2(i),koordCarboxyS2(i,:)
                    END Do
                END IF 
                
                CALL h_CarboxS2_propyl(atomsorteCarboxyS2,koordCarboxyS2,atomsorteH,koordH,size(atomsorteCarboxyS2),ziel,DEBUG)            
                
                IF (DEBUG) WRITE(*,*) 'End Carboxy1'
                
                !molA='.'
                !CALL make_coord_to_xyz_in_dir(atomsorte_H_all(2:23),koord_H_all(2:23,:),22,ziel,molA)
                
                !!!!Ende Test fuer ersten Polymerabschnitt poly1
            ELSE   ! Hauptteil, nicht erstes Element ist CarboxylS2
                
                IF (DEBUG) THEN
                    WRITE(*,*) 'Start Carb'
                    DO i=i_start+1,i_start+15 
                        WRITE(*,*) atomsorte(i),koord(i,:)
                    END DO
                END IF ! DEBUG 
                
                ! Abb 1 auf 2 ! Abbildungsvorschriften der einzelnen Carboxyl Abschnittel von 1-8 auf die Reihung der Atome im zweiten Abschnitt!! Abbildungen sind in Gruppen zusammengeordnet
                if (TRIM(adjustl(atomsorte(i_start+1)))=='F') THEN ! start ist i_start=0 also erstes Atom im Array
                    Aindex=transpose(reshape( (/ 13, 	1, 	12, 	2, 	2, 	3, 	1, 	4, 	3, 	5, 	4, 	6, 	5, 	7, 	6, 	8, 	7, 	9, 	8, 	10,&
                    & 	9, 	11, 	10, 	12, 	11, 	13, 	16, 	14, 	15, 	15 /), (/ 2, 15/)))
                    IF (DEBUG) WRITE(*,*) 'FALL_A'
                ELSE IF (TRIM(adjustl(atomsorte(i_start+1)))=='CH3') THEN ! start ist i_start=0 also erstes Atom im Array, Fall, wenn die Atome schon vorsortiert wurden
                    Aindex=transpose(reshape( (/ 2,1, 3,2, 4,3, 5,4, 6,5, 7,6, 8,7, 9,8, 10,9, 11,10, 12,11, 13,12, &
                    & 14,13, 15,14, 16,15 /), (/ 2, 15/))) ! Ende 14,1 C an der H-Start
                    IF (DEBUG) WRITE(*,*) 'FALL_A sorted'
                ! Abb 2 auf 2
                ELSE IF (TRIM(adjustl(atomsorte(i_start+7)))=='O') THEN
                    Aindex=transpose(reshape( (/ 1, 	1, 	2, 	2, 	3, 	3, 	4, 	4, 	5, 	5, 	6, 	6, 	7, 	7, 	8, 	8, 	9, 	9, 	10, 	10,&
                    & 	11, 	11, 	12, 	12, 	13, 	13, 	14, 	14, 	15, 	15 /), (/ 2, 15/)))
                    IF (DEBUG) WRITE(*,*) 'FALL_B'
                ELSE IF (TRIM(adjustl(atomsorte(i_start+12)))=='O') THEN
                    ! Abb 3/4/7/8 (vier Abschnitte als Abbildung ) auf 2
                    Aindex=transpose(reshape( (/ 1, 	1, 	15, 	2, 	13, 	3, 	14, 	4, 	6, 	5, 	7, 	6, 	12, 	7, 	8, 	8, 	9, 	9,&
                    & 	10, 10, 11, 	11, 	5, 	12, 	4, 	13, 	3, 	14, 	2, 	15 /), (/ 2, 15/)))
                    IF (DEBUG) WRITE(*,*) 'FALL_C'
                ELSE IF (TRIM(adjustl(atomsorte(i_start+11)))=='O') THEN
                ! Abb 5/6 auf 2
                    Aindex=transpose(reshape( (/ 1, 	1, 	2, 	2, 	3, 	3, 	4, 	4, 	5, 	5, 	6, 	6, 	11, 	7, 	7, 	8, 	8, 	9, 	9,&
                    & 	10, 	10, 	11, 	12, 	12, 	13, 	13, 	14, 	14, 	15, 	15 /), (/ 2, 15/)))
                    IF (DEBUG) WRITE(*,*) 'FALL_D'
                END IF ! Ende Abbildung der indizes um auf indexiserung im zweiten Abschnitt abzubilden.
                
            
                j=0
                IF (DEBUG) WRITE(*,*) 'i_start:',i_start,i_start+1,i_start+15
                DO i=i_start+1,i_start+15  !49,63 Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                        j=j+1
                        koordCarboxyS2(Aindex(j,2),:)=koord(i_start+Aindex(j,1),:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                        dummy = TRIM(adjustl(atomsorte(i_start+Aindex(j,1))))
                        atomsorteCarboxyS2(Aindex(j,2)) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
                        !WRITE(*,*) j,atomsorte(i),atomsorteCarboxyS2(j),koordCarboxyS2(j,:)
                END DO
                
                IF (TRIM(adjustl(atomsorte(i_start+1)))=='CH3') THEN ! Erstes Atom im Array bei Vorsortierung
                    i_start=i_start+16
                ELSE
                    i_start=i_start+15
                END IF
                CALL h_CarboxS2_propyl(atomsorteCarboxyS2,koordCarboxyS2,atomsorteH,koordH,size(atomsorteCarboxyS2),ziel,DEBUG)
            END IF !Ende Test fuer ersten Polymerabschnitt poly1 / Ende Hauptteil CarboxylS2
        END IF ! PBDT_chain(p)? 
        
        DO i=1,size(atomsorteH) ! Ubergabe an gesamtes file ab Count_atoms_all=2
            Count_atoms_all=Count_atoms_all+1
            atomsorte_H_all(Count_atoms_all)=atomsorteH(i)
            koord_H_all(Count_atoms_all,:)=koordH(i,:)
        END DO
        
        IF(HStart) THEN ! Anbau des ersten H Atomes
            ! H-Start
            
            IF( PBDT_chain(1) )THEN
                c1Atom=koordPBDT(1,:)     !  C-Atom am Randende
                SAtom=koordCarboxyS2(29,:)
                hAtom1(:)=koordPBDT(3,:)  ! Referenz C -Atom, welches im eingebauten Ring in der gegenueberliegenden Position (Parastellung) liegt
            ELSE 
                c1Atom=koordCarboxyS2(1,:)
                SAtom=koordCarboxyS2(15,:)
                hAtom1(:)=koordCarboxyS2(14,:)
            END IF    
            CALL setDistance(c1Atom,hAtom1,-1.095) !Abstand anpassen -- das negative Vorzeichen wandelt die Anbaurichtung um 
            !CALL setAngle(c1Atom,SAtom,hAtom1,hAtom1,(3.9))
            Rotwinkel=angle(SAtom,c1Atom,hAtom1)
            
            !IF (DEBUG) !WRITE(*,*) 'Winkel vorher', Rotwinkel
            CALL kreuzprodukt(SAtom-c1Atom, c1Atom-hAtom1,kreuzprod)
            hAtom1=hAtom1-c1Atom(:)
            call Rotationsmatrix_alpha(hAtom1,(119.1-Rotwinkel),kreuzprod-c1Atom,.true.)
            hAtom1(:)=hAtom1(:)+c1Atom(:) 
            !IF (DEBUG) !WRITE(*,*) 'Winkel HStart, hinterher:',angle(SAtom,c1Atom,hAtom1)
            IF (DEBUG) WRITE(*,*) 'HStart','H ',hatom1
            koord_H_all(1,:)=hAtom1
            atomsorte_H_all(1)='H  ' 
            HStart=.false.
        END IF
    
    END DO ! p-Schleife ueber die einzelnen Unterabschnitte
   

! H-Ende setze das letze H-Atom     
    IF( PBDT_chain( size(PBDT_chain)) )THEN
        c1Atom=koordPBDT(30,:)     ! C-Atom am Randende
        SAtom=koordCarboxyS2(29,:)
        hAtom2(:)=koordPBDT(15,:)  ! Referenz C -Atom, welches im eingebauten Ring in der gegenueberliegenden Position (Parastellung) liegt
    ELSE 
        c1Atom=koordCarboxyS2(14,:)
        SAtom=koordCarboxyS2(15,:)
        hAtom2(:)=koordCarboxyS2(1,:)
    END IF    
    CALL setDistance(c1Atom,hAtom2,-1.095) !Abstand anpassen -- das negative Vorzeichen wandelt die Anbaurichtung um 
    
    Rotwinkel=angle(SAtom,c1Atom,hAtom2)
    !CALL setAngle(c1Atom,SAtom,hAtom2,hAtom2,-119.9)
    
    CALL kreuzprodukt(SAtom-c1Atom, c1Atom-hAtom2,kreuzprod)
    hAtom2=hAtom2-c1Atom(:)
    call Rotationsmatrix_alpha(hAtom2,(Rotwinkel-119.1),kreuzprod-c1Atom,.true.)
    hAtom2(:)=hAtom2(:)+c1Atom(:) 

    Count_atoms_all=Count_atoms_all+1
    koord_H_all(Count_atoms_all,:)=hAtom2
    atomsorte_H_all(Count_atoms_all)='H  ' 
    IF(DEBUG) write(*,*) 'H-Winkel: ',angle(SAtom,c1Atom,hAtom2)


    !molA='.'
    !CALL make_coord_to_xyz_in_dir(atomsorte_H_all,koord_H_all,size(atomsorte_H_all),ziel,molA)
    
    IF (DEBUG) WRITE(*,*) 'h_PBDT_TS1 Ende'
    !DO j=1,size(atomsorteH),1
     !   WRITE(*,*)  j,atomsorteH(j),koordH(j,:)
    !END DO 
   
END SUBROUTINE h_PBDT_TS1


SUBROUTINE h_PBDT_CH3(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel,DEBUG) ! Schneidet die Ketten ab!
        !! Subroutine to add H-Atoms to form two methyl-CH3 group  to in h_PBDT_CH3-Part PBDT_TS1-molecule 
        IMPLICIT NONE
        Integer, INTENT(IN) :: MAXATM
        Real, Dimension(MAXATM,3), INTENT(IN) :: koord
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
        Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
	LOGICAL,INTENT(IN)::DEBUG
        Character(500) :: ziel
        Character(500) :: molA
        CHARACTER(1)   ::dummy
        Real, Dimension(3) ::cAtom,NAtom,X,hAtom1,hAtom2,hAtom3,kreuzprod,OAtom,hAtomA,OEsterAtom
        INTEGER::i,j,p,anzahl_methylgruppen=2
        ! Mindex C10-C9-O9-HAtomA-H1-H2-H3
        Integer, Dimension(2,7)::Mindex
        Integer, Dimension(3,2)::Cindex
        Real::Rotwinkel,dist_sum_1, dist_sum_2
        ALLOCATE(koordH(38,3))
        ALLOCATE(atomsorteH(38))
        
    	IF (DEBUG) WRITE(*,*) 'start h_PBDT_CH3'
        ! Die ersten 4 Zahlen geben jeweils den alten Index an (C9-S8-C7-S6) -- dann kommen die neuen Indizes (H10-H11-H12)
        ! 2 Methylgruppen werden gesetzt.
        Mindex=transpose(reshape( &
    &(/ 9,8,7,6,10,11,12  ,21,20,19,18,25,26,27 /), (/ 7, 2/)))
    ! Uebergabe in neue Form
    j=0
    DO i=1,9  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
            j=j+1
            koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
            dummy = TRIM(adjustl(atomsorte(i)))
            atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
            !WRITE(*,*) j,atomsorteH(j),koordH(j,:)
    END DO 
    j=12 ! ueberspringe 3 H-Atome
    DO i=10,21  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
        j=j+1
        koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
        dummy = TRIM(adjustl(atomsorte(i)))
        atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
        !WRITE(*,*) i,atomsorteH(j),koordH(j,:)
    END DO        
    
    j=27 ! ueberspringe 3 H-Atome
    DO i=22,32  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
        j=j+1
        koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
        dummy = TRIM(adjustl(atomsorte(i)))
        atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
        !WRITE(*,*) i,atomsorteH(j),koordH(j,:)
    END DO   

    
   ! Schleife zum Anbau der 2 Methylendgruppen ( Hier gibt es natuerlich kein NAtom!!!!)
   DO p=1,anzahl_methylgruppen     !! Vorsicht die Benennung der Atome stimmt ab hier nicht mehr mit den Benennungen im h_PBDT_CH3 ueberein, da die DIPBI Routine verwendet wurde 
        cAtom=koord(Mindex(p,1),:)
        NAtom=koord(Mindex(p,2),:)
        CALL kreuzprodukt(CAtom - NAtom, cAtom-koord(Mindex(p,3),:),kreuzprod)
        hAtom1 = NAtom
        CALL setDistance(cAtom,hAtom1,1.095)
        CALL rotate(cAtom,(cAtom+kreuzprod),hAtom1,112.5)
        hAtom2 = hAtom1
        hAtom3 = hAtom1
        CALL rotate(NAtom,cAtom,hAtom2,120.) ! Rotationen von hAtom2
        CALL rotate(NAtom,cAtom,hAtom3,240.) ! Rotation von hAtom3
    
       !Start Rotation um Dih zum H-Atom anzupassen
        call set_dih_by_rotation_around_axis(hAtom1,cAtom,NAtom,koord(Mindex(p,4),:),60.0,Rotwinkel)

        !call set_dih_by_rotation_around_axis(hAtom1,cAtom,NAtom,koord(Mindex(p,4),:),7.8,Rotwinkel)
        X=hAtom2(:)-cAtom(:)
        call Rotationsmatrix_alpha(X,Rotwinkel,NAtom-cAtom,.true.)
        hAtom2(:)=X(:)+cAtom(:)
 
        X=hAtom3(:)-cAtom(:)
        call Rotationsmatrix_alpha(X,Rotwinkel,NAtom-cAtom,.true.)
        hAtom3(:)=X(:)+cAtom(:)
        !Ende Extra Rotation
    
        koordH(Mindex(p,5),:)=hAtom1(:)
        koordH(Mindex(p,6),:)=hAtom2(:)
        koordH(Mindex(p,7),:)=hAtom3(:)
        atomsorteH(Mindex(p,5))='H '
        atomsorteH(Mindex(p,6))='H '
        atomsorteH(Mindex(p,7))='H '

   END DO ! i-Schleife Methylgruppen
   
    !DO j=1,38
    !    WRITE(*,*) j,atomsorteH(j),koordH(j,:)
    !END DO  
    !   molA='.'
    !CALL make_coord_to_xyz_in_dir(atomsorteH,koordH,size(atomsorteH),ziel,molA)
   IF(DEBUG) WRITE(*,*) 'Ende h_PBDT_CH3'

END SUBROUTINE h_PBDT_CH3



SUBROUTINE h_CarboxS2_propyl(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel,debug) ! Schneidet die Ketten ab!
        !! Subroutine to add H-Atoms to form a propyl group (CH3-CH-CH3)-O-C=O --Ring to in CarboxS2-Part(with F-Atom) PBDT_TS1-molecule 
        IMPLICIT NONE
        Integer, INTENT(IN) :: MAXATM
        Real, Dimension(MAXATM,3), INTENT(IN) :: koord
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
        Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
	LOGICAL, INTENT(IN) :: debug
        Character(500) :: ziel
        Character(500) :: molA
        CHARACTER(1)   ::dummy
        Real, Dimension(3) ::cAtom,c2Atom,c3Atom,NAtom,X,hAtom1,hAtom2,hAtom3,kreuzprod,kreuzprodB,OAtom,hAtomA,hAtomB,OEsterAtom
        INTEGER::i,j,p,anzahl_methylgruppen=2
        ! Mindex C10-C9-O9-HAtomA-H1-H2-H3
        Integer, Dimension(2,7)::Mindex
        Integer, Dimension(3,2)::Cindex
        Real::Rotwinkel,dRotwinkel,winkel_A,winkel_B,winkel_opt,dist_sum_1, dist_sum_2
        INTEGER::winkel,i_start,i_ende
        ALLOCATE(koordH(22,3))
        ALLOCATE(atomsorteH(22))
    IF(debug) WRITE(*,*) 'start h_CarboxS2_propyl'
        ! Die ersten 4 Zahlen geben jeweils den alten Index an (C10-C9-O8-HAtomA) -- dann kommen die neuen Indizes (H12-H13-H14)
        ! 2 Methylgruppen werden gesetzt.
        Mindex=transpose(reshape( &
    &(/ 10,9,8,6,12,13,14  ,11,9,8,6,16,17,18 /), (/ 7, 2/)))
    !!! Uebergabe des C-Atome aus der Seitenkette an die Position im Array 1=alter Index ,2= neuer Index Anzahl der C-Atome durchlaufen als i=1,6 
    !!  C10 aus koord wird zum ersten neuen C in koordH C11 daher 10,11 ... dann wird Platz fuer 3 H-Atome gelassen (neuer Index 12 13 14) und altes C11 wird dann zu C15 mit (H16 H17 H18)
    Cindex=transpose(reshape( (/ 9,9, 10,11, 11,15 /), (/ 2, 3/)))
    Do i=1,3
        !WRITE(*,*) Cindex(i,1),Cindex(i,2)
        koordH(Cindex(i,2),:)=koord(Cindex(i,1),:)
        atomsorteH(Cindex(i,2))='C'
    END DO
    j=0
    DO i=1,9  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
            j=j+1
            koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
            dummy = TRIM(adjustl(atomsorte(i)))
            atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
            !WRITE(*,*) j,atomsorteH(j),koordH(j,:)
    END DO 
    j=18
    DO i=12,15  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
        j=j+1
        koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
        dummy = TRIM(adjustl(atomsorte(i)))
        atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
        !WRITE(*,*) j,atomsorteH(j),koordH(j,:)
    END DO        
            
    !H-Start
    cAtom=koord(9,:)   !  C-Atom am Randstart
    c2Atom=koord(10,:) 
    c3Atom=koord(11,:) 
    OEsterAtom=koord(8,:)  
    hAtomA(:)=c2Atom
    OAtom=koord(7,:)
    
    
    ! Test fuer positive oder negative Drehrichtung 
    ! Erst werden zwei testwinkel fuer hAtomA und hAtomB bestimmt um zu schauen, in welchen Bereich die Rotation durchgefuehrt werden muss
    Rotwinkel=120
    hAtomA(:)=c2Atom
    CALL setDistance(cAtom,hAtomA,1.095)
    hAtomA=hAtomA-cAtom(:)
    call Rotationsmatrix_alpha(hAtomA,Rotwinkel,OEsterAtom(:)-cAtom(:),.true.)
    hAtomA(:)=hAtomA(:)+cAtom(:)
    
    Rotwinkel=240
    hAtomB(:)=c2Atom
    CALL setDistance(cAtom,hAtomB,1.095) 
    hAtomB=hAtomB-cAtom(:)
    call Rotationsmatrix_alpha(hAtomB,Rotwinkel,OEsterAtom(:)-cAtom(:),.true.)
    hAtomB(:)=hAtomB(:)+cAtom(:)
    
    ! Dann werden fuer die neuen H-Atome die Abstaende zu C2 und C3 gemessen und der Bereich, fuer den die Abstaende groesser ausfallen wird gewaehlt
    !WRITE (*,*) 'Dist 1' ,(norm(c2Atom-hAtomA)+norm(c3Atom-hAtomA))
    !WRITE (*,*) 'Dist 2' ,(norm(c2Atom-hAtomB)+norm(c3Atom-hAtomB))
    
    
     IF( (norm(c2Atom-hAtomA)+norm(c3Atom-hAtomA)) > (norm(c2Atom-hAtomB)+norm(c3Atom-hAtomB)) ) THEN
        i_start=90
        i_ende=140
        !WRITE(*,*) 'FALL A'
     ELSE
        i_start=220
        i_ende=270
        !WRITE(*,*) 'FALL B'
    END IF

    ! Es soll nun durch schrittweisen Drehen, aus der Position con C2 um die C1-OEsterAtom-Achse, die Rotation gefunden werden, fuer
    ! die die Differenz (dRotwinkel) der beiden Winkel winkel_A(C2-C1-hAtomA) und winkel_B(hAtomA-C1-C3) minimal wird.
    winkel_opt=real(i_start)
    winkel_A=angle(c2Atom,cAtom,hAtomA)
    winkel_B=angle(hAtomA,cAtom,c3Atom)
    dRotwinkel=abs(winkel_A-winkel_B)
    DO winkel=i_start,i_ende,2
        Rotwinkel=REal(winkel)
        hAtomA(:)=c2Atom !Reset
        ! Finde den passenden Winkel so, dass die Abweiung dRotwinkel minimal wird
        hAtomA=hAtomA-cAtom(:)
        call Rotationsmatrix_alpha(hAtomA,Rotwinkel,OEsterAtom(:)-cAtom(:),.true.)
        hAtomA(:)=hAtomA(:)+cAtom(:)
        
        winkel_A=angle(c2Atom,cAtom,hAtomA)
        winkel_B=angle(hAtomA,cAtom,c3Atom)
        IF ( abs(winkel_A-winkel_B) < dRotwinkel ) THEN ! Winkel sollen gleich gross sein
            ! winkel_opt enthaelt zum Schluss die optimale Rotation
            winkel_opt=Rotwinkel
            dRotwinkel=abs(winkel_A-winkel_B)
        END IF
        
        !WRITE(*,*) 'Winkel',Rotwinkel,winkel_A,winkel_B,winkel_opt
    END DO 
    
    ! Rotation mit optimalem Winkel: winkel_opt
    hAtomA(:)=c2Atom
    CALL setDistance(cAtom,hAtomA,1.095) 
    hAtomA=hAtomA-cAtom(:)
    call Rotationsmatrix_alpha(hAtomA,winkel_opt,OEsterAtom(:)-cAtom(:),.true.)
    hAtomA(:)=hAtomA(:)+cAtom(:)
    koordH(10,:)=hAtomA 
    atomsorteH(10)='H  ' 
    !winkel_A=angle(c2Atom,cAtom,hAtomA)
    !winkel_B=angle(hAtomA,cAtom,c3Atom)
    !WRITE(*,*) 'Final Winkel',Rotwinkel,winkel_A,winkel_B,winkel_opt
    
   !!! Schleife zum Anbau der 2 Methylendgruppen ( Hier gibt es natuerlich kein NAtom!!!!)
   DO p=1,anzahl_methylgruppen     !! Vorsicht die Benennung der Atome stimmt ab hier nicht mehr mit den Benennungen im CarboxyS2 ueberein, da die DIPBI Routine verwendet wurde 
        cAtom=koord(Mindex(p,1),:)
        NAtom=koord(Mindex(p,2),:)
        CALL kreuzprodukt(CAtom - NAtom, cAtom-koord(Mindex(p,3),:),kreuzprod)
        hAtom1 = NAtom
        CALL setDistance(cAtom,hAtom1,1.095)
        CALL rotate(cAtom,(cAtom+kreuzprod),hAtom1,112.5)
        hAtom2 = hAtom1
        hAtom3 = hAtom1
        CALL rotate(NAtom,cAtom,hAtom2,120.) ! Rotationen von hAtom2
        CALL rotate(NAtom,cAtom,hAtom3,240.) ! Rotation von hAtom3
    
       !Start Rotation um Dih zum H-Atom anzupassen
        call set_dih_by_rotation_around_axis(hAtom1,cAtom,NAtom,hAtomA,60.0,Rotwinkel)

        !call set_dih_by_rotation_around_axis(hAtom1,cAtom,NAtom,koord(Mindex(p,4),:),7.8,Rotwinkel)
        X=hAtom2(:)-cAtom(:)
        call Rotationsmatrix_alpha(X,Rotwinkel,NAtom-cAtom,.true.)
        hAtom2(:)=X(:)+cAtom(:)
 
        X=hAtom3(:)-cAtom(:)
        call Rotationsmatrix_alpha(X,Rotwinkel,NAtom-cAtom,.true.)
        hAtom3(:)=X(:)+cAtom(:)
        !Ende Extra Rotation
    
        koordH(Mindex(p,5),:)=hAtom1(:)
        koordH(Mindex(p,6),:)=hAtom2(:)
        koordH(Mindex(p,7),:)=hAtom3(:)
        atomsorteH(Mindex(p,5))='H '
        atomsorteH(Mindex(p,6))='H '
        atomsorteH(Mindex(p,7))='H '

   END DO ! i-Schleife Methylgruppen

END SUBROUTINE h_CarboxS2_propyl






SUBROUTINE h_diPBI_kette(atomsorte,koord,atomsorteH,koordH,MAXATM,ziel) 
!! Adds H Atoms to the total side chain in DIPBI
IMPLICIT NONE
Integer, INTENT(IN) :: MAXATM
Real, Dimension(MAXATM,3), INTENT(IN) :: koord
Character(5), Dimension(MAXATM), INTENT(INOUT) :: atomsorte
REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
Character(500), INTENT(IN) :: ziel
Integer :: i,j,ind,anzahl,anzahlneu,anzahlringe,ind5RingEnde,cAtomRingInd,sAtomInd,zaehler
 Real, Dimension(3) :: cAtom,cAtomRing,hAtom1,hAtom2,hAtom3,vorigesC,naechstesC,kreuzprod,hAtomold
 Integer, Dimension(48) :: ind5Ring,ind_C_Ketten
 LOGICAL:: Cvorher, erste_vorher
 Integer:: ketten, erste_Kette, ketten_2
 Integer, Dimension(4):: ind_erste
 Real::Rotwinkel

!write(*,*) 'MAXATM', MAXATM
anzahl=MAXATM
anzahlneu = MAXATM + 100 !Gesamtanzahl Atome in neuer Datei
!write(*,*) 'Anzahl neu H-DIPBI_KETTE',anzahlneu ! 218 ?
allocate(atomsorteH(anzahlneu))
allocate(koordH(anzahlneu,3))
! Übergabe der ersten 74 Atome des DIPBI ohne Seitenketten - jeweils mit N und erstem C-Atom
do j=1,74,1
    koordH(j,:)=koord(j,:)    
    atomsorteH(j)=atomsorte(j)   
END DO ! Wert von j nach der Schleife ist j=75

! Indizes der 48 C Atome in den 4 Seitenketten
ind_C_Ketten(1)=1   ! C
do i=0,10,1
    ind_C_Ketten(i+2)=75+i ! von 75-85
END DO
ind_C_Ketten(13)=37     ! C
do i=0,10,1
ind_C_Ketten(i+14)=86+i ! von 86-96
END DO
ind_C_Ketten(25)=38   ! C
do i=0,10,1
ind_C_Ketten(i+26)=97+i ! von 97-107
END DO
ind_C_Ketten(37)=74   ! C
do i=0,10,1
ind_C_Ketten(i+38)=108+i ! von 108-118
END DO


atomsorte(1)="CH2"
atomsorte(37)="CH2"
atomsorte(38)="CH2"
atomsorte(74)="CH2"
do i=75,118,1
    atomsorte(i)="CH2"
END DO 
atomsorte(85)="CH3"
atomsorte(96)="CH3"
atomsorte(107)="CH3"
atomsorte(118)="CH3"

cvorher=.false.  !Logikvariable für das erste CH2 in jeder Kette
ketten=0 		 !Abzählen der 4 C12-Ketten
erste_vorher= .false. 	!Logikvariable für das zweite CH2 jeder Kette
ketten_2=1

!neue Koordinaten schreiben 
!Die entsprechenden H werden jeweils nach dem C-Atom angefuegt
DO i=1,size(ind_C_Ketten),1
        If (TRIM(atomsorte(ind_C_Ketten(i)))=='CH3' .OR. (TRIM(atomsorte(ind_C_Ketten(i)))=='CH2')) Then
                IF( ind_C_Ketten(i)/=1 .AND. ind_C_Ketten(i)/=37 .AND. ind_C_Ketten(i)/=38 .AND. ind_C_Ketten(i)/=74 ) THEN
                koordH(j,:)=koord(ind_C_Ketten(i),:)    ! Ausnahme der 4 schon geschriebenen ersten C-Atome der Kette
                atomsorteH(j)=' C'        ! Zuerst jeweils das C in der Seitenkette an die Stelle j schreiben
                j=j+1
                END IF
        End If

	IF(atomsorte(ind_C_Ketten(i)) == 'CH2') THEN !H-Atome an Hexylgruppe anfügen
        !Falls es das erste CH2 der Kette ist:
    	IF(cvorher .eqv. .false.) Then  ! Beginn am Kettenanfang 
			sAtomInd=ind_C_Ketten(i)
            ind_erste(ketten+1)=ind_C_Ketten(i)
            if (ketten==0) Then
              erste_Kette=ind_C_Ketten(i)
            End if

            !1. C-Atom mit H bestuecken
			DO ind=1,3,1 !Koordinaten merken
				cAtom(ind) = koord(sAtomInd,ind)
				naechstesC(ind) = koord(ind_C_Ketten(i+1),ind)
                if ((ketten == 0) .or. (ketten == 2)) Then
					vorigesC(ind) = koord(sAtomInd +1,ind)   ! Auswahl von C1 und C38 
                    Else
                    vorigesC(ind) = koord(sAtomInd -2,ind)    ! Auswahl von C37 und C74
                End if

			END DO
            ketten=ketten+1
            !Write(*,*) 'Ketten', ketten
			if(ketten==4) then
             			erste_vorher = .true.
						cvorher = .true.
            End if
        
            CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
            hAtom1 = cAtom + kreuzprod ! Vektor senkrecht zur Ebene der 3 C-Atome
            hAtom2 = cAtom - kreuzprod
        
            CALL rotate(naechstesC,vorigesC,hAtom1,32.)  !Winkel anpassen
            CALL rotate(naechstesC,vorigesC,hAtom2,-32.) ! Rotation nach außen um die Verbinungdachse der Cs
        
            CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen
            CALL setDistance(cAtom,hAtom2,1.095) !Abstand anpassen
        
            !Write(2,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
            koordH(j,:)=hAtom1(:)    
            atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
            j=j+1
            koordH(j,:)=hAtom2(:)    
            atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
            j=j+1
            cvorher=.true.
            
        Else 

            if (erste_vorher .eqv. .true.) Then
                !2. C-Atom mit H bestuecken - hier muess wieder auf das 1. C zugegriffen werden 
                !Daher wird eine Adressierung mit ketten_2 eingefuehrt
                DO ind=1,3,1 !Koordinaten merken
                    cAtom(ind) = koord(ind_C_Ketten(i),ind)
                    vorigesC(ind) = koord(ind_erste(ketten_2+1),ind) !ketten_2=2 muss erst bei der zweiten Kette loslaufen
                    naechstesC(ind) = koord(ind_C_Ketten(i+1),ind)   ! daher ist ketten_2=1 zu Beginn 
                END DO
                          
                erste_vorher = .false.
                ketten_2=ketten_2+1

                
                CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
    
                hAtom1 = cAtom + kreuzprod
                hAtom2 = cAtom - kreuzprod
 
                CALL rotate(naechstesC,vorigesC,hAtom1,32.) !Winkel anpassen
                CALL rotate(naechstesC,vorigesC,hAtom2,-32.)
                
            
                CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen
                CALL setDistance(cAtom,hAtom2,1.095) !Abstand anpassen
    
                !Write(2,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
                !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
                !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
                    koordH(j,:)=hAtom1(:)    
                    atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                    j=j+1
                    koordH(j,:)=hAtom2(:)    
                    atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                    j=j+1       
            
            Else    ! Durchlaufen vom 3. C-Atom der Kette bis zum vorletzen C-Atom.
                    !1. C-Atom mit H bestuecken
                DO ind=1,3,1 !Koordinaten merken
                    cAtom(ind) = koord(ind_C_Ketten(i),ind)
                    vorigesC(ind) = koord(ind_C_Ketten(i-1),ind)
                    naechstesC(ind) = koord(ind_C_Ketten(i+1),ind)
                END DO
                        
                CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
                hAtom1 = cAtom + kreuzprod
                hAtom2 = cAtom - kreuzprod
                
                CALL rotate(naechstesC,vorigesC,hAtom1,32.) !Winkel anpassen
                CALL rotate(naechstesC,vorigesC,hAtom2,-32.)
        
                CALL setDistance(cAtom,hAtom1,1.095) !Abstand anpassen
                CALL setDistance(cAtom,hAtom2,1.095) !Abstand anpassen
                hAtomOld=hAtom1
                
                !Write(2,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
                !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
                !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
                        koordH(j,:)=hAtom1(:)    
                        atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                        j=j+1
                        koordH(j,:)=hAtom2(:)    
                        atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                        j=j+1
                        
                
            END IF!Write(*,*) 'CH2 Ende'
        End if
	Else
        If(atomsorte(ind_C_Ketten(i))=='CH3') THEN 
            cAtomRingInd=  ind_C_Ketten(i - 2)
            sAtomInd= ind_C_Ketten(i)
            !letztes C-Atom mit 3 H-Atomen bestücken
            DO ind=1,3,1 !Koordinaten merken
                cAtom(ind) = koord(sAtomInd,ind)
                vorigesC(ind) = koord(cAtomRingInd+1,ind)
                naechstesC(ind) = koord(sAtomInd+1,ind) !vorvoriges C-atom
            END DO		
            
            CALL kreuzprodukt(cAtom - vorigesC, cAtom -naechstesC,kreuzprod)
    
            hAtom1 = vorigesC
            CALL setDistance(cAtom,hAtom1,1.095)
            CALL rotate(cAtom,cAtom+kreuzprod,hAtom1,112.)
    
            hAtom2 = hAtom1  !schon mal die ungefähre Position
            hAtom3 = hAtom1
    
            CALL rotate(vorigesC,cAtom,hAtom2,120.) !in korrekte Position bringen
            CALL rotate(vorigesC,cAtom,hAtom3,240.)
            
           !Rotation um C-C Achse sodass der Diederwinkel H-C-C-H ~60° wird
           ! (Korrektur zur passenden Newman-Projektion)
           call set_dih_by_rotation_around_axis(hAtom1,cAtom,vorigesC,hAtomold,60.,Rotwinkel)

           hAtom2=hAtom2-cAtom(:)
           call Rotationsmatrix_alpha(hAtom2,Rotwinkel,vorigesC(:)-cAtom(:),.true.)
           hAtom2(:)=hAtom2(:)+cAtom(:)
          
           hAtom3=hAtom3-cAtom(:)
           call Rotationsmatrix_alpha(hAtom3,Rotwinkel,vorigesC(:)-cAtom(:),.true.)
           hAtom3(:)=hAtom3(:)+cAtom(:)

            !letztes C-Atom ausgeben
            !Write(2,'(2X,A2,1X,3F16.6)') 'C ',cAtom(1),cAtom(2),cAtom(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom1(1),hAtom1(2),hAtom1(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom2(1),hAtom2(2),hAtom2(3)
            !Write(2,'(2X,A2,1X,3F16.6)') 'H ',hAtom3(1),hAtom3(2),hAtom3(3)
                    koordH(j,:)=hAtom1(:)    
                    atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                    j=j+1
                    koordH(j,:)=hAtom2(:)    
                    atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                    j=j+1
                    koordH(j,:)=hAtom3(:)    
                    atomsorteH(j)=' H'        ! Zuerst jeweils H in der Seitenkette an die Stelle j schreiben
                    j=j+1
            
            cvorher=.false. !Nächste Kette
            erste_vorher=.true. !Nächste Kette
        END IF
    End IF
END DO
 !CLOSE(2)

END SUBROUTINE h_diPBI_kette



SUBROUTINE h_cut_chain_HDI(atomsorte,koord,MAXATM,atomsorteH,koordH,ziel) ! Schneidet die Ketten ab!
        !! Subroutine to cut the chains and substitute by methyl groups (CH3) to a HDI-molecule 
        !! obtained by a GROMACS-file.
        IMPLICIT NONE
        Integer, INTENT(IN) :: MAXATM
        Real, Dimension(MAXATM,3), INTENT(IN) :: koord
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
        Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
        Character(500) :: ziel
        Character(500) :: molA
        CHARACTER(1)::dummy
        Real, Dimension(3) ::CAtom,hAtom1
        INTEGER::i,j,p
        LOGICAL::DEBUG=.false.
 

        IF( MAXATM==170) THEN ! HDI
                WRITE(*,*) 'h_cut_chain_HDI' 
                ALLOCATE(koordH(56,3))
                ALLOCATE(atomsorteH(56))
                j=0
                DO i=1,56  !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                        j=j+1
                        koordH(j,:)=koord(i,:)           !Koordinaten an neuen Array übergeben nicht die Seitenketten! 
                        dummy = TRIM(adjustl(atomsorte(i)))
                        atomsorteH(j) = dummy(1:1)//'  ' ! Schneidet alles rechts ab
                END DO

                ! C28 durch H ersetzen 
                CAtom=koordH(25,:)   !  C-Atom am Randende
                hAtom1(:)=koordH(28,:) 
                CALL setDistance(CAtom,hAtom1,1.095) !Abstand anpassen -- das negative Vorzeichen wandelt die Anbaurichtung um 
                koordH(28,:)= hAtom1
                atomsorteH(28)='H ' 
                
                ! C31 durch H ersetzen 
                CAtom=koordH(29,:)   !  C-Atom am Randende
                hAtom1(:)=koordH(31,:) 
                CALL setDistance(CAtom,hAtom1,1.095) !Abstand anpassen -- das negative Vorzeichen wandelt die Anbaurichtung um 
                koordH(31,:)= hAtom1
                atomsorteH(31)='H ' 
                
                IF (DEBUG) THEN
                    DO i=1,56
                        WRITE(*,*) atomsorteH(i),koordH(i,:)
                    END DO 
                    WRITE(*,*) 'Ende h_cut_chain_HDI'
                END IF
        ELSE ! Not the correct size 
                ! Uebergabe ohne Bearbeitung der H-Atome
                ALLOCATE(koordH(MAXATM,3))
                ALLOCATE(atomsorteH(MAXATM))
                DO i=1,MAXATM
                    koordH(i,:)=koord(i,:)
                    atomsorteH(i)=atomsorte(i)
                END DO    
                write(*,*) ' Keine Bearbeitung der H-Atome in h_cut_chain_HDI'
        END IF ! h_cut_chain_HDI
END SUBROUTINE h_cut_chain_HDI

SUBROUTINE alq3_rename_elements(atomsorte,koord,MAXATM,atomsorteH,koordH,ziel) 
        !! Subroutine alq3_rename_elements, H-Atoms or bonds are not affected!
        !! obtained by a GROMACS-file.
        IMPLICIT NONE
        Integer, INTENT(IN) :: MAXATM
        Real, Dimension(MAXATM,3), INTENT(IN) :: koord
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
        REAL, ALLOCATABLE, DIMENSION(:,:), INTENT (OUT) ::koordH
        Character(5),ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteH
        Character(500) :: ziel
        Character(500) :: molA
        INTEGER::i
        LOGICAL::DEBUG=.false.
	IF(MAXATM==52)THEN ! check for Alq3 number of elements
		IF (DEBUG) WRITE(*,*) 'Start alq3_rename_elements'
		ALLOCATE(koordH(52,3))
        	ALLOCATE(atomsorteH(52))
		DO i=1,MAXATM
			koordH(i,:)=koord(i,:)
			IF (index(atomsorte(i),'H') .ne. 0) THEN
				 atomsorteH(i)='H'
			ELSE IF (index(atomsorte(i),'C') .ne. 0) THEN
				 atomsorteH(i)='C'
			ELSE IF (index(atomsorte(i),'N') .ne. 0) THEN
				 atomsorteH(i)='N'
			ELSE IF (index(atomsorte(i),'O') .ne. 0) THEN
				 atomsorteH(i)='O'
			ELSE IF (index(atomsorte(i),'Al') .ne. 0) THEN
				!WRITE(*,*) 'Warning: use H Atom instead of Al in Alq3!'
				atomsorteH(i)='Al'   ! 'Al'  ! use H as a dummy element for Alq3 calculations
			ELSE
				WRITE(*,*) 'Error: Element was not expected in Alq3 molecule: CALL EXIT (1)'
				CALL EXIT (1)		
			END IF ! check elements
		END DO ! loop element array
		
		IF (DEBUG) THEN
                    DO i=1,MAXATM
                        WRITE(*,*) atomsorteH(i),koordH(i,:)
                    END DO 
                    WRITE(*,*) 'Ende alq3_rename_elements'
                END IF
	ELSE
		WRITE(*,*) 'Error: molecule is not Alq3'
	END IF
END SUBROUTINE alq3_rename_elements


Real Function norm(atom1)
IMPLICIT NONE
Real, Dimension(3), Intent(IN) :: atom1
norm = SQRT(atom1(1)**2 + atom1(2)**2 + atom1(3)**2)
END Function norm


SUBROUTINE SortArray(Array1,n_max)
!! Subroutine um eine Real Array aufsteigend zu ordnen
!! Array1 liest den Array ein
!! n_max gibt die Anzahl der sortierten Elemente an
IMPLICIT NONE
INTEGER::i,j,i_min,n_max
Real, Dimension(:), INTENT (INOUT)::Array1 
REAL::dummy

! Selection Sort
!n_max=size(Array1)
do j=1,n_max-1 ! Durchgehen der Atome in Molekuel     
    i_min=j
        DO i=j+1,n_max ! Durchgehen von B          
           IF( Array1(i) < Array1(i_min) ) THEN
               i_min=i
           END IF
        END DO    
        IF( i_min /= j ) THEN ! tauschen
               dummy=Array1(j)
               Array1(j)=Array1(i_min)
               Array1(i_min)=dummy
        END IF           
END DO

END SUBROUTINE SortArray

REAL FUNCTION angle(pointA,pointB,pointC)
!! calculates the angle beween 3 Points in A,B,C ( A-to-B-to-C)
IMPLICIT NONE
Real, Dimension(3), INTENT (IN)::pointA,pointB,pointC
REAL, DIMENSION(3)::V1,V2
REAL, PARAMETER::PI=4.0*atan(1.0)
V1=pointA-pointB
V2=pointC-pointB
angle=acos( DOT_PRODUCT(V1,V2)/(norm(V1)*norm(V2)) )*180.0/PI
END FUNCTION angle

Subroutine setDistance(festesAtom,varAtom,newDistance) !! varAtom verschieben Richtung beibehalten
IMPLICIT NONE
Real, Dimension(3) :: festesAtom,varAtom
Real :: newDistance,oldDistance
Integer :: j

oldDistance = getDistance(festesAtom,varAtom)
DO j=1,3,1
	varAtom(j) = newDistance/oldDistance *(varAtom(j)-festesAtom(j)) + festesAtom(j)
END DO

RETURN
END Subroutine setDistance

Real Function getDistance(atom1,atom2)
IMPLICIT NONE
Real, Dimension(3), Intent(IN) :: atom1,atom2
getDistance = SQRT((atom1(1) - atom2(1))**2 + (atom1(2) - atom2(2))**2 + (atom1(3) - atom2(3))**2)
END Function getDistance

Subroutine setAngle(Atom1,Atom2,Atom3,SetAtom4,winkel)
!! Atom1: Zentrum der Drehung!!!
!! Atom2 und Atom3 sind Punkte, die die Ebene n1 mit Atom1 aufspannen.
!! winkel: Rotationswinkel um Normalenvektor der Ebene n1
!! SetAtom4 wird gedreht
!! n1= Normalenvekotor der Ebene 
!! RotAtom=dummyvektor
IMPLICIT NONE
Real, Dimension(3), INTENT(IN)  ::Atom1,Atom2,Atom3
Real, Dimension(3), INTENT(OUT) ::SetAtom4
Real, INTENT(IN) :: winkel
Real, Dimension(3) :: n1,RotAtom

CALL kreuzprodukt(Atom2(:)-Atom1(:),Atom3-Atom1(:),n1)

RotAtom=SetAtom4(:)-Atom1(:)
CALL Rotationsmatrix_alpha(RotAtom,winkel,n1,.true.)
SetAtom4(:)=RotAtom(:)+Atom1(:)

END Subroutine 

Subroutine rotate(achse1,achse2,varAtom,winkel) !rotiert varAtom um winkel Grad um Achse
IMPLICIT NONE
Real, Dimension(3) :: achse1,achse2,varAtom,achsvektor,rotiervektor,kreuzprod
Real :: winkel,bogen,skalar,achsLaenge,cosinus,sinus
Integer :: i

bogen=winkel/180. *ACOS(-1.)  !Winkel in bogenmass umrechnen
cosinus = COS(bogen)
sinus = SIN(bogen)

!achsvektor, rotiervektor und skalarprodukt berechnen
skalar =0.
achsLaenge = getDistance(achse1,achse2)
DO i=1,3,1
	achsvektor(i) = (achse2(i) - achse1(i))/achsLaenge !normierter AchsVektor
	rotiervektor(i) = varAtom(i) - achse2(i)
	skalar = skalar + achsvektor(i)*rotiervektor(i)
END DO

!kreuzprodukt berechnen
 CALL kreuzprodukt(achsvektor,rotiervektor,kreuzprod)

!neue Koordinaten
DO i=1,3,1
	varAtom(i) = (1-cosinus)*achsvektor(i)*skalar + cosinus*rotiervektor(i) +sinus*kreuzprod(i) + achse2(i)
END DO


END Subroutine rotate

SUBROUTINE Rotationsmatrix_alpha(x,winkel,Achsenvektor,beliebigeachse) 
        !!Rotation von Vektor x um winkel um eine Ursprungsgerade um die Achse mit dem Einheistvekotor n=(n_1,n_2,n_3)
        !!winkel in Grad
        !! .true. fuer Rotation um beliebige Achse n
        !! .false. fuer Ursprungsrotation um Richtung n 
        IMPLICIT NONE
        Real, INTENT(IN)::winkel
        Real, Dimension(3), INTENT(INOUT) ::x
        Real, Dimension(3), INTENT(IN) ::Achsenvektor
        Real, Dimension(3,3) ::Rotmat
        REAL, PARAMETER::Pi=4.0*atan(1.0)
        Real, Dimension(3)::n,cross1,cross2
        REAL::alpha
        LOGICAL, INTENT(IN)::beliebigeachse  ! .true. sonst Ursprungsrotation
        alpha=Pi/180.0*winkel
        n(:)=Achsenvektor(:)
        n(:)=n(:)/norm(n(:))  !Erzeuge Einheitsvekor

        IF ( beliebigeachse .EQV. .false.) THEN ! Berechnung der Rotation um Ursprungsgerade       
                RotMat(1,1)=n(1)**2*(1-cos(alpha))+cos(alpha)
                RotMat(1,2)=n(1)*n(2)*(1-cos(alpha))-n(3)*sin(alpha)
                RotMat(1,3)=n(1)*n(3)*(1-cos(alpha))+n(2)*sin(alpha)

                RotMat(2,1)=n(1)*n(2)*(1-cos(alpha))+n(3)*sin(alpha)
                RotMat(2,2)=n(2)**2*(1-cos(alpha))+cos(alpha)
                RotMat(2,3)=n(2)*n(3)*(1-cos(alpha))-n(1)*sin(alpha)


                RotMat(3,1)=n(1)*n(3)*(1-cos(alpha))-n(2)*sin(alpha)
                RotMat(3,2)=n(2)*n(3)*(1-cos(alpha))+n(1)*sin(alpha)
                RotMat(3,3)=n(3)**2*(1-cos(alpha))+cos(alpha)
        
                !write(*,*) 'Determinante der Rotationsmatrix',Det_3x3(RotMat)        
                !write(*,*) x
                !DO i=1,3
                !        write(*,*) (RotMat(i,j) j=1,3)
                !END DO
                IF ( (1.0-Det_3x3(Rotmat)) > 1.0E-4 ) THEN
                    write(*,*) 'Fehler bei der Berechnung der Rotationsmatrix: Det(RotMat) /= 1.0 '
                    Write(*,*) 'Beende Programm'
                    CALL EXIT (1)
                END IF        
                x=MATMUL(RotMat,x)
                RETURN
        END IF
        
        ! Drehung um eine beliebige Achse n
        ! R_{\hat{n}}(\alpha)\vec{x}=\hat{n}(\hat{n}\cdot\vec{x})
        ! +\cos\left(\alpha\right)(\hat{n}\times\vec{x})\times\hat{n}
        ! +\sin\left(\alpha\right)(\hat{n}\times\vec{x})       
        IF ( beliebigeachse .EQV. .true. ) THEN  ! Drehung um Beliebige Achse
                call kreuzprodukt(n,x,cross1)
                call kreuzprodukt(cross1,n,cross2)
                x(:)=n(:)*DOT_PRODUCT(n,x)+cos(alpha)*cross2(:)+sin(alpha)*cross1(:)
                RETURN
        END IF

END SUBROUTINE Rotationsmatrix_alpha

REAL FUNCTION Det_3x3(Mat)
!! Berechnen Determinante der 3x3 Matrix mit der Regel von Sarrus
       IMPLICIT NONE
       REAL, DIMENSION(3,3), INTENT(IN)::Mat
 Det_3x3=MAT(1,1)*MAT(2,2)*MAT(3,3)-MAT(3,1)*MAT(2,2)*MAT(1,3)+&
         MAT(1,2)*MAT(2,3)*MAT(3,1)-MAT(3,2)*MAT(2,3)*MAT(1,1)+&
         MAT(1,3)*MAT(2,1)*MAT(3,2)-MAT(3,3)*MAT(2,1)*MAT(1,2)
END FUNCTION 

SUBROUTINE set_dih_by_rotation_around_axis(atom1,atom2,atom3,atom4,dih_ziel,Rotwinkel)
 !! Rotiert ein Atom1 (H) um die Achse Atom2-Atom3 bis der diederwinkel 
 !! zu Atom 4 dem gewuenschten Diederwinkel:dih_ziel (Grand) A1-A2-A3-A4 entspricht.
 IMPLICIT NONE
 REAL, Dimension(3), INTENT(INOUT):: atom1
 REAL, Dimension(3), INTENT(IN):: atom2,atom3,atom4
 Real, INTENT(IN) ::dih_ziel
 Real, INTENT(OUT)::Rotwinkel
 Real::dih_best,dih,alpha
 Real, Dimension(3)::X,C
 Integer :: i

! Bestimmung des besten Rotwinkels

 call calc_dihedral(atom4,atom3,atom2,atom1,dih_best)
 !write(*,*) "dih start",dih_best
 Rotwinkel=0.0
 C(:)=atom2(:)
 do i=1,720,1    ! Durchlaufen der Rotationen in 1-Schritten um besten Rotationswinkel zu finden.
        alpha=0.5*i
        X=atom1(:)-C(:)
        call Rotationsmatrix_alpha(X,alpha,atom3(:)-C(:),.true.)
        call calc_dihedral(atom4,atom3,atom2,X(:)+C(:),dih)
        IF ( abs(abs(dih)-abs(dih_ziel)) < abs(abs(dih_best)-abs(dih_ziel)) ) THEN ! Abweichung von dih_ziel soll gering sein. 
                dih_best=dih
                Rotwinkel=alpha    
        END IF
end do

! write(*,*) 'dih=',dih,'dih_ziel=',dih_ziel,'Rotationswinkel=',Rotwinkel
! write(*,*) 'dih_best=',dih_best
 !Rotation von atom1 um Rotwinkel um Achse (Atom2-Atom3)

 X=atom1(:)-C(:)
 call Rotationsmatrix_alpha(X,Rotwinkel,atom3(:)-C(:),.true.)
 atom1(:)=X(:)+C(:)

  RETURN
END SUBROUTINE set_dih_by_rotation_around_axis 


SUBROUTINE calc_dihedral(A1,A2,A3,A4,diederwinkel)
!! Berechnet den Diederwinkel zwischen 4 Atomen
!! A1,A2,A3,A4 Eingangsvektoren
!! B1,B2,B3    Verbindungsvektoren
!! N1,N2       Ebenen-Normalenvektoren
!! M           Orthogonalvektor, fuer lok. cart. Koordinatensystem
IMPLICIT NONE
REAL, DIMENSION (3), INTENT(IN) ::A1,A2,A3,A4
REAL, INTENT (OUT)  :: diederwinkel
REAL, DIMENSION (3) :: B1,B2,B3,N1,N2,M
REAL::PI=4.0*atan(1.0)
! Verbindungsvekotren festlegen
B1(:)=A2(:)-A1(:)
B2(:)=A3(:)-A2(:)
B3(:)=A4(:)-A3(:)
! Vektoren Normieren
B1=B1/norm(B1)
B2=B2/norm(B2)
B3=B3/norm(B3)
! Ebenen Vektoren
call kreuzprodukt(B1,B2,N1)
call kreuzprodukt(B2,B3,N2)
call kreuzprodukt(B2,N1,M)

diederwinkel=180.0/PI*atan2(DOT_PRODUCT(M,N2),DOT_PRODUCT(N1,N2))

END SUBROUTINE calc_dihedral        


Subroutine kreuzprodukt(vektor1,vektor2,kreuzprod) !berechnet vektor1 x vektor2 in kreuzprod
IMPLICIT NONE
Real, Dimension(3) :: vektor1, vektor2,kreuzprod

kreuzprod(1) = vektor1(2)*vektor2(3)-vektor1(3)*vektor2(2)
kreuzprod(2) = vektor1(3)*vektor2(1)-vektor1(1)*vektor2(3)
kreuzprod(3) = vektor1(1)*vektor2(2)-vektor1(2)*vektor2(1)

END Subroutine kreuzprodukt


SUBROUTINE how_many_h_atoms_to_mol(R_u_I,atomsorte,residue_names,NAtoms,N_Resids,GesamtanzahlH,KETTE)
    !! SUBROUTINE to determine the number of all atoms, if H-Atoms are included
    !! OUTPUT:  GesamtanzahlH
    IMPLICIT NONE
    INTEGER ::NAtoms,N_Resids
    INTEGER, Dimension(N_Resids+1,2), INTENT(IN)    :: R_u_I
    Character(5), Dimension(NAtoms), INTENT(IN)     :: atomsorte  
    Character(5),Dimension(N_Resids), INTENT(IN)    :: residue_names
    LOGICAL, INTENT(IN)    :: KETTE
    INTEGER, INTENT(OUT)   :: GesamtanzahlH
    Character(5) :: residue_name
    INTEGER::i,j,k,NAtomsH
    LOGICAL::DEBUG=.true.
    LOGICAL::  DIPBI,DIPBI_KETTE,P3HT,P3MT,P3HT_KETTE,PPDI,PPDI_KETTE,PBDT_TS1,HDI_KETTE,PPDI_Propyl
    !!! fuer h_P3HT
    INTEGER::ind,ind5RingEnde,anzahl,anzahlringe,anzahlneu,N_CarboxyS2,N_PBDT
    LOGICAL::HStart=.false.,HEnde=.false.
    INTEGER, ALLOCATABLE, DIMENSION(:) ::ind5Ring
	
	
    GesamtanzahlH=0
    NAtomsH=0
    ! Schleife, um alles Resids durchzulaufen und um daraus die Anzahl der Atome mit Wasserstoff abgesaettigt zu bestimmen. 
    DO i=1,N_Resids
		 	residue_name=residue_names(i)
		 	IF ('THP'==TRIM(residue_name(1:3))) THEN     !!1) P3HT / P3MT Check residue_name   
			    IF(KETTE) THEN   !!!!! P3HT  !!!!!!
					HStart=.false.
					HEnde=.false.

					! Bestimmung ob erstes und letztes H-Atom existiert
					IF(TRIM(adjustl(atomsorte(R_u_I(i,2))))=='HC' .AND. &
						&TRIM(adjustl(atomsorte(R_u_I(i,2)+1)))=='CA') HStart=.true.
					IF(TRIM(adjustl(atomsorte(R_u_I(i+1,2)-8))) == 'HC' ) HEnde=.true.

					!Indexarray fur 5Ringe und Wasserstoffe erzeugen
					anzahl=R_u_I(i+1,2)-R_u_I(i,2)
					ALLOCATE(ind5Ring(anzahl))
					k =1
					ind=1
					ind5Ring=0
					DO WHILE(ind <= anzahl)
						ind5Ring(k) = ind
						IF(index(atomsorte(R_u_I(i,2)+ind),'S') >0) THEN !Hexylgruppe ueberspringen
							ind = ind +6 !jetzt: ind letztes Atom der Hexylgruppe
						END IF
						k=k+1
						ind=ind+1
					END DO
					ind5RingEnde = k-1 !letzter belegter Platz

					!Write(*,*) 'HStart',HStart,'HEnde',HEnde

					IF( HStart .AND. HEnde ) THEN
					    anzahlringe = (ind5RingEnde -2)/6                   !anzahlringe = (ind5RingEnde -2)/6 !fuer ganze Kette 
					    anzahlneu = anzahl + anzahlringe*13                 !Gesamtanzahl Atome in neuer Datei
					ELSE
					    IF(HStart .OR. HEnde) THEN
						anzahlringe = (ind5RingEnde -1)/6  
						anzahlneu = anzahl + anzahlringe*13 + 1 !Gesamtanzahl Atome in neuer Datei +1 fuer Anfangs/Abschluss H-Atom
					    ELSE 
						anzahlringe = (ind5RingEnde)/6 
						anzahlneu = anzahl + anzahlringe*13 + 2 !Gesamtanzahl Atome in neuer Datei +2 fuer Anfangs+Abschluss H-Atom
					    END IF
					END IF
					DEALLOCATE(ind5Ring)
					NAtomsH=anzahlneu
			    ELSE !! h_P3MT  !!! P3MT
					!!! P3MT Determine number of atoms in P3MT
					HStart=.false.
					HEnde=.false.

					! Bestimmung ob erstes und letztes H-Atom existiert
	      				IF( TRIM(adjustl(atomsorte(R_u_I(i,2)))) == 'HC' .AND. &
						&TRIM(adjustl(atomsorte(R_u_I(i,2)+1))) == 'CA' ) HStart=.true.
					IF( TRIM(adjustl(atomsorte(R_u_I(i+1,2)-8))) == 'HC' ) HEnde=.true.

					!Indexarray fur 5Ringe und Wasserstoffe erzeugen
					anzahl=R_u_I(i+1,2)-R_u_I(i,2)
					ALLOCATE(ind5Ring(anzahl))
					k=1
					ind=1
					ind5Ring=0
					DO WHILE(ind <= anzahl)
						ind5Ring(k) = ind
						IF(index(atomsorte(R_u_I(i,2)+ind),'S') >0) THEN !Hexylgruppe ueberspringen
							ind = ind +6 !jetzt: ind letztes Atom der Hexylgruppe
						END IF
						k=k+1
						ind=ind+1
					END DO
					ind5RingEnde = k-1 !letzter belegter Platz
					IF( HStart .AND. HEnde ) THEN
					    anzahlringe = (ind5RingEnde -2)/6              !anzahlringe = (ind5RingEnde -2)/6 !fuer ganze Kette 
					    anzahlneu = anzahl - anzahlringe*5 + anzahlringe*3 !Gesamtanzahl Atome in neuer Datei, 5 C weg, 3 H dazu je Ring
					ELSE
					    IF(HStart .OR. HEnde) THEN
						anzahlringe = (ind5RingEnde -1)/6  
						anzahlneu = anzahl - anzahlringe*5 + anzahlringe*3 + 1 !Gesamtanzahl Atome in neuer Datei, 5 C weg, 3 H dazu je Ring +1 fuer Anfangs/Abschluss H-Atom
					    ELSE 
						anzahlringe = (ind5RingEnde)/6 
						anzahlneu = anzahl - anzahlringe*5 + anzahlringe*3 + 2 !Gesamtanzahl Atome in neuer Datei, 5 C weg, 3 H dazu je Ring +2 fuer Anfangs+Abschluss H-Atom
					    END IF
					END IF
					DEALLOCATE(ind5Ring)
					NAtomsH=anzahlneu
					IF(DEBUG) THEN
						WRITE(*,*) 'P3MT NAtomsH 1:',NAtomsH,' anzahlringe:',anzahlringe,' ind',ind,HStart,HEnde
					END IF
			    END IF    ! h_P3MT oder h_P3HT mit Kette


			ELSE IF( TRIM(residue_name)=='DIPBI' ) THEN  ! DIPBI
			    IF(KETTE) THEN !h_diPBI_kette 
				NAtomsH=218
			    ELSE  ! Methyl-Gruppe
				NAtomsH=86
			    END IF  
			ELSE IF ('HDI'==TRIM(residue_name(1:3)) ) THEN
			    IF(KETTE)THEN
				anzahl=R_u_I(i+1,2)-R_u_I(i,2) 
				IF( anzahl .LE. 57) THEN
					NAtomsH=56
				ELSE
					NAtomsH=170
				END IF   
			    ELSE ! h_cut_chain_HDI => NatomsH=56
				NAtomsH=56
			    END IF
			ELSE IF ('PPDI'==TRIM(residue_name(1:4)) ) THEN
			    IF(KETTE) THEN
				NAtomsH=144
			    ELSE IF(PPDI_Propyl) THEN
				NAtomsH=96
			    ELSE ! H-Atoms fuer h_PPDI 
				NAtomsH=78
			    END IF !PPDI
			ELSE IF ( '8poly'==TRIM(residue_name(1:5)) ) THEN
			    ! Zaehle die Anzahl der Elemente im Molekuel N_CarboxyS2 ueber Anzahl der F-Atome und N_PBDT=(N_HC-Gruppen/6) uber die Anzahl der HC-Gruppen
			    N_CarboxyS2=0
			    N_PBDT=0
			    DO j=R_u_I(i,2),R_u_I(i+1,2)-1  ! TO DO INDEX FALSCH!!!
				if (TRIM(adjustl(atomsorte(j)))=='F')  N_CarboxyS2=N_CarboxyS2+1
				if (TRIM(adjustl(atomsorte(j)))=='HC') N_PBDT=N_PBDT+1
			    END DO ! Zahle Anzahl
			    N_PBDT =N_PBDT/6
			    IF (DEBUG) WRITE(*,*) ' Die Anzahl der Carboxy-Gruppen ist: ',N_CarboxyS2
			    IF (DEBUG) WRITE(*,*) ' Die Anzahl der N_PBDT-Gruppen ist: ',N_PBDT   
			    NAtomsH=38*N_PBDT+22*N_CarboxyS2+2
			    IF (DEBUG) WRITE(*,*) ' Die Gesamtanzahl nach dem Anfuegen der H-Atome betraegt',NAtomsH
			ELSE
			    NAtomsH=R_u_I(i+1,2)-R_u_I(i,2)
			END IF 	 !! Check residue_name   	


		GesamtanzahlH=GesamtanzahlH+NAtomsH
    END DO !NAtoms
    IF(DEBUG) THEN
	WRITE(*,*)  'how_many_h_atoms_to_mol GesamtanzahlH: ', GesamtanzahlH
    END IF
END SUBROUTINE how_many_h_atoms_to_mol


INTEGER FUNCTION Calc_N_electrons_in(atomsorte,N_Atoms,Kette) ! 
!! Berechnet der Anzahl der Elektronen        
!! Wenn Kette=false, dann wird die Ladung fuer CH1,CH2,CH3 in gro format =0 gesetzte
        IMPLICIT NONE
        Integer, INTENT(IN) :: N_Atoms
        Character(5), Dimension(N_Atoms), INTENT(IN) :: atomsorte
        Real, Dimension(3)::schwerpunkt
        integer::i
        REAL::m_i,m_ges,z_ges,z_i
        LOGICAL, INTENT(IN)::Kette 

       ! Berechnung der Anzahl der Elektronen im neutralen Molekül z_ges 
        z_ges=0
        DO i=1,N_Atoms
           z_i=0
           ! .gro Gruppen
           if (TRIM(adjustl(atomsorte(i)))=='HC')  z_i=1    ! fuer H
           if (KETTE) THEN 
                if (TRIM(adjustl(atomsorte(i)))=='CH1') z_i=12   ! fuer 1H und 1C 
                if (TRIM(adjustl(atomsorte(i)))=='CH2') z_i=14   ! fuer 2H und 1C  
                if (TRIM(adjustl(atomsorte(i)))=='CH3') z_i=15   ! fuer 3H und 1C  
           else 
                if (TRIM(adjustl(atomsorte(i)))=='CH1') z_i=0   ! fuer 1H und 1C 
                if (TRIM(adjustl(atomsorte(i)))=='CH2') z_i=0   ! fuer 2H und 1C  
                if (TRIM(adjustl(atomsorte(i)))=='CH3') z_i=0   ! fuer 3H und 1C  
           end if
           if (TRIM(adjustl(atomsorte(i)))=='CA')  z_i=12   ! fuer CA = C  
           if (TRIM(adjustl(atomsorte(i)))=='NR')  z_i=7   ! fuer NR = N 
           if (TRIM(adjustl(atomsorte(i)))=='NT')  z_i=7   ! fuer NT = N 
           ! Normale Atome
           if (TRIM(adjustl(atomsorte(i)))=='H')  z_i=1!Wasserstoff     H       1
           if (TRIM(adjustl(atomsorte(i)))=='He') z_i=2!Helium          He      2
           if (TRIM(adjustl(atomsorte(i)))=='Li') z_i=3!Lithium         Li      3
           if (TRIM(adjustl(atomsorte(i)))=='Be') z_i=4!Beryllium       Be      4
           if (TRIM(adjustl(atomsorte(i)))=='B')  z_i=5!Bor     B       5
           if (TRIM(adjustl(atomsorte(i)))=='C')  z_i=6!Kohlenstoff     C       6
           if (TRIM(adjustl(atomsorte(i)))=='N')  z_i=7!        Stickstoff      N       7
           if (TRIM(adjustl(atomsorte(i)))=='O')  z_i=8!        Sauerstoff      O       8
           if (TRIM(adjustl(atomsorte(i)))=='F')  z_i=9!        Fluor   F       9
           if (TRIM(adjustl(atomsorte(i)))=='Ne') z_i=10!        Neon    Ne      10
           if (TRIM(adjustl(atomsorte(i)))=='Na') z_i=11!        Natrium         Na      11
           if (TRIM(adjustl(atomsorte(i)))=='Mg') z_i=12!Magnesium       Mg      12
           if (TRIM(adjustl(atomsorte(i)))=='Al') z_i=13!        Aluminum        Al      13
           if (TRIM(adjustl(atomsorte(i)))=='Si') z_i=14!        Silizium        Si      14
           if (TRIM(adjustl(atomsorte(i)))=='P')  z_i=15!        Phosphor        P       15
           if (TRIM(adjustl(atomsorte(i)))=='S')  z_i=16!Schwefel        S       16
           if (TRIM(adjustl(atomsorte(i)))=='Cl') z_i=17!Chlor   Cl      17
           if (TRIM(adjustl(atomsorte(i)))=='K')  z_i=19!        Kalium  K       19
           if (TRIM(adjustl(atomsorte(i)))=='Ar') z_i=17!Argon   Ar      18
           if (TRIM(adjustl(atomsorte(i)))=='Ca') z_i=20!Kalzium         Ca      20
           if (TRIM(adjustl(atomsorte(i)))=='Sc') z_i=21!        Scandium        Sc      21
           if (TRIM(adjustl(atomsorte(i)))=='Ti') z_i=22!Titan   Ti      22
           if (TRIM(adjustl(atomsorte(i)))=='V')  z_i=23 !        Vanadium        V       23
           if (TRIM(adjustl(atomsorte(i)))=='Cr') z_i=24!        Chrom   Cr      24
           if (TRIM(adjustl(atomsorte(i)))=='Mn') z_i=25!Mangan  Mn      25
           if (TRIM(adjustl(atomsorte(i)))=='Fe') z_i=26!Eisen   Fe      26
           if (TRIM(adjustl(atomsorte(i)))=='Ni') z_i=28!        Nickel  Ni      28
           if (TRIM(adjustl(atomsorte(i)))=='Co') z_i=27!        Kobalt  Co      27
           if (TRIM(adjustl(atomsorte(i)))=='Cu') z_i=29!Kupfer  Cu      29
           if (TRIM(adjustl(atomsorte(i)))=='Zn') z_i=30!Zink    Zn      30
           if (TRIM(adjustl(atomsorte(i)))=='Ga') z_i=31!Gallium         Ga      31
           if (TRIM(adjustl(atomsorte(i)))=='Ge') z_i=32!Germanium       Ge      32
           if (TRIM(adjustl(atomsorte(i)))=='As') z_i=33!        Arsen   As      33
           if (TRIM(adjustl(atomsorte(i)))=='Se') z_i=34!Selen   Se      34
           if (TRIM(adjustl(atomsorte(i)))=='Br') z_i=35!Brom    Br      35
           if (TRIM(adjustl(atomsorte(i)))=='Kr') z_i=36!Krypton         Kr      36
           if (TRIM(adjustl(atomsorte(i)))=='Rb') z_i=37!        Rubidium        Rb      37
           if (TRIM(adjustl(atomsorte(i)))=='Sr') z_i=38!Strontium       Sr      38
           if (TRIM(adjustl(atomsorte(i)))=='Y')  z_i=39!        Yttrium         Y       39
           if (TRIM(adjustl(atomsorte(i)))=='Zr') z_i=40!Zirkonium       Zr      40
           if (TRIM(adjustl(atomsorte(i)))=='Nb') z_i=41!        Nobium  Nb      41
           if (TRIM(adjustl(atomsorte(i)))=='Mo') z_i=42!Molybdän        Mo      42
           if (TRIM(adjustl(atomsorte(i)))=='Tc') z_i=43!Technetium      Tc      43
           if (TRIM(adjustl(atomsorte(i)))=='Ru') z_i=44!Ruthenium       Ru      44
           if (TRIM(adjustl(atomsorte(i)))=='Rh') z_i=45!        Rhodium         Rh      45
           if (TRIM(adjustl(atomsorte(i)))=='Pd') z_i=46!Palladium       Pd      46
           if (TRIM(adjustl(atomsorte(i)))=='Ag') z_i=47!        Silber  Ag      47
           if (TRIM(adjustl(atomsorte(i)))=='Cd') z_i=48!        Kadmium         Cd      48
           if (TRIM(adjustl(atomsorte(i)))=='In') z_i=49!        Indium  In      49
           if (TRIM(adjustl(atomsorte(i)))=='Sn') z_i=50!Zinn    Sn      50
           if (TRIM(adjustl(atomsorte(i)))=='Sb') z_i=51!Antimon         Sb      51
           if (TRIM(adjustl(atomsorte(i)))=='I')  z_i=53!        Jod     I       53
           if (TRIM(adjustl(atomsorte(i)))=='Te') z_i=52!Tellur  Te      52
           if (TRIM(adjustl(atomsorte(i)))=='Xe') z_i=54!        Xenon   Xe      54
           if (TRIM(adjustl(atomsorte(i)))=='Re') z_i=75!        Rhenium         Re      75
           if (TRIM(adjustl(atomsorte(i)))=='Os') z_i=76!Osmium  Os      76
           if (TRIM(adjustl(atomsorte(i)))=='Ir') z_i=77!        Iridium         Ir      77
           if (TRIM(adjustl(atomsorte(i)))=='Pt') z_i=78!        Platin  Pt      78
           if (TRIM(adjustl(atomsorte(i)))=='Au') z_i=79!        Gold    Au      79
           if (TRIM(adjustl(atomsorte(i)))=='Hg') z_i=80!Quecksilber     Hg      80
           if (TRIM(adjustl(atomsorte(i)))=='Tl') z_i=81!        Thallium        Tl      81
           if (TRIM(adjustl(atomsorte(i)))=='Pb') z_i=82!Blei    Pb      82
           if (TRIM(adjustl(atomsorte(i)))=='Bi') z_i=83!        Wismut  Bi      83
           if (TRIM(adjustl(atomsorte(i)))=='Po') z_i=84    !Polonium        Po      84
           if (TRIM(adjustl(atomsorte(i)))=='At') z_i=85     !Astat   At      85
           if (TRIM(adjustl(atomsorte(i)))=='Rn') z_i=86      !    Radon
           
            IF (z_i==0 .and. KETTE) THEN
                WRITE(*,*) 'Es ist ein Fehler bei der Berechnung der Elektronenanzahl z_ges aufgetreten!'
            END IF
            z_ges=z_ges+z_i
        END DO
	Calc_N_electrons_in=z_ges
END FUNCTION Calc_N_electrons_in

REAL FUNCTION element_to_mass(atomsorte,Kette)
   !! gets the Element symbol and returns the atomic mass in a.u.
   !! Default Kette= .false.
   IMPLICIT NONE
   Character(LEN=*) , INTENT(IN):: atomsorte
   LOGICAL, INTENT(IN) :: Kette
   REAL::mass
   mass=0.0
   if (Kette) then
      if (TRIM(adjustl(atomsorte))=='CH1') THEN
                     element_to_mass=13.0186 !  fuer 1H und 1C 
                     return 
      ELSE if (TRIM(adjustl(atomsorte))=='CH2') THEN
                     element_to_mass=14.0265 ! fuer 2H und 1C  
                      return 
      ELSE if (TRIM(adjustl(atomsorte))=='CH3') THEN
                     element_to_mass=15.0344 ! fuer 3H und 1C  
                      Return 
      END IF
   END IF ! Kette = .true.
   
   SELECT CASE (TRIM(adjustl(atomsorte)))
         CASE('HC')  
            mass=1.0079  ! fuer H
         CASE('CH1') 
                     mass=0.0 !  fuer 1H und 1C
         CASE('CH2') 
                     mass=0.0 ! fuer 2H und 1C
         CASE('CH3') 
                     mass=0.0 ! fuer 3H und 1C
           CASE('CA')  
            mass=12.0107 ! fuer CA = C 
           CASE('NR')  
            mass=14.0067 ! fuer NR = N 
           CASE('NT')  
            mass=14.0067 ! fuer NT = N 
           !richtige Elemente
           CASE('H')  
            mass=1.0079   !Wasserstoff     H       1
           CASE('He')  
            mass=4.0026  !Helium          He      2
           CASE('Li')  
            mass=6.941   !Lithium         Li      3
           CASE('Be')  
            mass=9.0122  !Beryllium       Be      4
           CASE('B')  
            mass=10.811   !Bor     B       5
           CASE('C')  
            mass=12.0107  !Kohlenstoff     C       6
           CASE('N')  
            mass=14.0067  !        Stickstoff      N       7
           CASE('O')  
            mass=15.9994  !        Sauerstoff      O       8
           CASE('F')  
            mass=18.9984  !        Fluor   F       9
           CASE('Ne')  
            mass=20.1797 !        Neon    Ne      10
           CASE('Na')  
            mass=22.9897 !        Natrium         Na      11
           CASE('Mg')  
            mass=24.305  !Magnesium       Mg      12
           CASE('Al')  
            mass=26.9815 !        Aluminum        Al      13
           CASE('Si')  
            mass=28.0855 !        Silizium        Si      14
           CASE('P')  
            mass=30.9738  !        Phosphor        P       15
           CASE('S')  
            mass=32.065   !Schwefel        S       16
           CASE('Cl')  
            mass=35.453  !Chlor   Cl      17
           CASE('K')  
            mass=39.0983  !        Kalium  K       19
           CASE('Ar')  
            mass=39.948  !Argon   Ar      18
           CASE('Ca')  
            mass=40.078  !Kalzium         Ca      20
           CASE('Sc')  
            mass=44.9559 !        Scandium        Sc      21
           CASE('Ti')  
            mass=47.867  !Titan   Ti      22
           CASE('V')  
            mass=50.9415  !        Vanadium        V       23
           CASE('Cr')  
            mass=51.9961 !        Chrom   Cr      24
           CASE('Mn')  
            mass=54.938  !Mangan  Mn      25
           CASE('Fe')  
            mass=55.845  !Eisen   Fe      26
           CASE('Ni')  
            mass=58.6934 !        Nickel  Ni      28
           CASE('Co')  
            mass=58.9332  !        Kobalt  Co      27
           CASE('Cu')  
            mass=63.546   !Kupfer  Cu      29
           CASE('Zn')  
            mass=65.39    !Zink    Zn      30
           CASE('Ga')  
            mass=69.723   !Gallium         Ga      31
           CASE('Ge')  
            mass=72.64    !Germanium       Ge      32
           CASE('As')  
            mass=74.9216  !        Arsen   As      33
           CASE('Se')  
            mass=78.96    !Selen   Se      34
           CASE('Br')  
            mass=79.904   !Brom    Br      35
           CASE('Kr')  
            mass=83.8     !Krypton         Kr      36
           CASE('Rb')  
            mass=85.4678  !        Rubidium        Rb      37
           CASE('Sr')  
            mass=87.62    !Strontium       Sr      38
           CASE('Y')  
            mass=88.9059   !        Yttrium         Y       39
           CASE('Zr')  
            mass=91.224   !Zirkonium       Zr      40
           CASE('Nb')  
            mass=92.9064  !        Nobium  Nb      41
           CASE('Mo')  
            mass=95.94   !Molybdän        Mo      42
           CASE('Tc')  
            mass=98      !Technetium      Tc      43
           CASE('Ru')  
            mass=101.07  !Ruthenium       Ru      44
           CASE('Rh')  
            mass=102.9055!        Rhodium         Rh      45
           CASE('Pd')  
            mass=106.42  !Palladium       Pd      46
           CASE('Ag')  
            mass=107.8682!        Silber  Ag      47
           CASE('Cd')  
            mass=112.411 !        Kadmium         Cd      48
           CASE('In')  
            mass=114.818 !        Indium  In      49
           CASE('Sn')  
            mass=118.71  !Zinn    Sn      50
           CASE('Sb')  
            mass=121.76  !Antimon         Sb      51
           CASE('I')  
            mass=126.9045!        Jod     I       53
           CASE('Te')  
            mass=127.6   !Tellur  Te      52
           CASE('Xe')  
            mass=131.293 !        Xenon   Xe      54
           CASE('Re')  
            mass=186.207 !        Rhenium         Re      75
           CASE('Os')  
            mass=190.23  !Osmium  Os      76
           CASE('Ir')  
            mass=192.217 !        Iridium         Ir      77
           CASE('Pt')  
            mass=195.078 !        Platin  Pt      78
           CASE('Au')  
            mass=196.9665!        Gold    Au      79
           CASE('Hg')  
            mass=200.59  !Quecksilber     Hg      80
           CASE('Tl')  
            mass=204.3833!        Thallium        Tl      81
           CASE('Pb')  
            mass=207.2   !Blei    Pb      82
           CASE('Bi')  
            mass=208.9804!        Wismut  Bi      83
           CASE('Po')  
            mass=209     !Polonium        Po      84
           CASE('At')  
            mass=210 !Astat   At      85
           CASE('Rn')  
            mass=222 !    Radon
         CASE DEFAULT
            mass=0.0
            WRITE(*,*) 'WARNING! Element symbol was not recognized. Mass is set to zero for: '//TRIM(adjustl(atomsorte))
         END SELECT    
      element_to_mass=mass
END FUNCTION element_to_mass



SUBROUTINE Massenschwerpunkt(atomsorte,koord,MAXATM,schwerpunkt,Kette)! 
! subroutine berechnet den Massenschwerpunkt mit Ausgabe in schwerpunkt        
        IMPLICIT NONE
        Integer, INTENT(IN) :: MAXATM
        Real, Dimension(MAXATM,3), INTENT(IN) :: koord
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
        logical, INTENT(IN)::Kette
        Real, Dimension(3)::schwerpunkt
        integer::i
        REAL::m_i,m_ges
        !write(*,*) 'Start der Berechnung des Schwerpunktes'        
        schwerpunkt(:)=0
        m_ges=0
        m_i=0
        do i=1,MAXATM
           m_i=0
           ! GRO Ausgaben
           if (TRIM(adjustl(atomsorte(i)))=='HC')  m_i=1.0079  ! fuer H
           if (Kette) then
                   if (TRIM(adjustl(atomsorte(i)))=='CH1') m_i=13.0186 !  fuer 1H und 1C 
                   if (TRIM(adjustl(atomsorte(i)))=='CH2') m_i=14.0265 ! fuer 2H und 1C  
                   if (TRIM(adjustl(atomsorte(i)))=='CH3') m_i=15.0344 ! fuer 3H und 1C  
           else
                   if (TRIM(adjustl(atomsorte(i)))=='CH1') m_i=0.0 !  fuer 1H und 1C
                   if (TRIM(adjustl(atomsorte(i)))=='CH2') m_i=0.0 ! fuer 2H und 1C
                   if (TRIM(adjustl(atomsorte(i)))=='CH3') m_i=0.0 ! fuer 3H und 1C
           end if
           if (TRIM(adjustl(atomsorte(i)))=='CA')  m_i=12.0107 ! fuer CA = C 
           if (TRIM(adjustl(atomsorte(i)))=='NR')  m_i=14.0067 ! fuer NR = N 
           if (TRIM(adjustl(atomsorte(i)))=='NT')  m_i=14.0067 ! fuer NT = N 
           if (TRIM(adjustl(atomsorte(i)))=='OE')  m_i=15.9994 ! fuer OE = O 
           !richtige Elemente
           if (TRIM(adjustl(atomsorte(i)))=='H')  m_i=1.0079   !Wasserstoff     H       1
           if (TRIM(adjustl(atomsorte(i)))=='He')  m_i=4.0026  !Helium          He      2
           if (TRIM(adjustl(atomsorte(i)))=='Li')  m_i=6.941   !Lithium         Li      3
           if (TRIM(adjustl(atomsorte(i)))=='Be')  m_i=9.0122  !Beryllium       Be      4
           if (TRIM(adjustl(atomsorte(i)))=='B')  m_i=10.811   !Bor     B       5
           if (TRIM(adjustl(atomsorte(i)))=='C')  m_i=12.0107  !Kohlenstoff     C       6
           if (TRIM(adjustl(atomsorte(i)))=='N')  m_i=14.0067  !        Stickstoff      N       7
           if (TRIM(adjustl(atomsorte(i)))=='O')  m_i=15.9994  !        Sauerstoff      O       8
           if (TRIM(adjustl(atomsorte(i)))=='F')  m_i=18.9984  !        Fluor   F       9
           if (TRIM(adjustl(atomsorte(i)))=='Ne')  m_i=20.1797 !        Neon    Ne      10
           if (TRIM(adjustl(atomsorte(i)))=='Na')  m_i=22.9897 !        Natrium         Na      11
           if (TRIM(adjustl(atomsorte(i)))=='Mg')  m_i=24.305  !Magnesium       Mg      12
           if (TRIM(adjustl(atomsorte(i)))=='Al')  m_i=26.9815 !        Aluminum        Al      13
           if (TRIM(adjustl(atomsorte(i)))=='Si')  m_i=28.0855 !        Silizium        Si      14
           if (TRIM(adjustl(atomsorte(i)))=='P')  m_i=30.9738  !        Phosphor        P       15
           if (TRIM(adjustl(atomsorte(i)))=='S')  m_i=32.065   !Schwefel        S       16
           if (TRIM(adjustl(atomsorte(i)))=='Cl')  m_i=35.453  !Chlor   Cl      17
           if (TRIM(adjustl(atomsorte(i)))=='K')  m_i=39.0983  !        Kalium  K       19
           if (TRIM(adjustl(atomsorte(i)))=='Ar')  m_i=39.948  !Argon   Ar      18
           if (TRIM(adjustl(atomsorte(i)))=='Ca')  m_i=40.078  !Kalzium         Ca      20
           if (TRIM(adjustl(atomsorte(i)))=='Sc')  m_i=44.9559 !        Scandium        Sc      21
           if (TRIM(adjustl(atomsorte(i)))=='Ti')  m_i=47.867  !Titan   Ti      22
           if (TRIM(adjustl(atomsorte(i)))=='V')  m_i=50.9415  !        Vanadium        V       23
           if (TRIM(adjustl(atomsorte(i)))=='Cr')  m_i=51.9961 !        Chrom   Cr      24
           if (TRIM(adjustl(atomsorte(i)))=='Mn')  m_i=54.938  !Mangan  Mn      25
           if (TRIM(adjustl(atomsorte(i)))=='Fe')  m_i=55.845  !Eisen   Fe      26
           if (TRIM(adjustl(atomsorte(i)))=='Ni')  m_i=58.6934 !        Nickel  Ni      28
           if (TRIM(adjustl(atomsorte(i)))=='Co')  m_i=58.9332  !        Kobalt  Co      27
           if (TRIM(adjustl(atomsorte(i)))=='Cu')  m_i=63.546   !Kupfer  Cu      29
           if (TRIM(adjustl(atomsorte(i)))=='Zn')  m_i=65.39    !Zink    Zn      30
           if (TRIM(adjustl(atomsorte(i)))=='Ga')  m_i=69.723   !Gallium         Ga      31
           if (TRIM(adjustl(atomsorte(i)))=='Ge')  m_i=72.64    !Germanium       Ge      32
           if (TRIM(adjustl(atomsorte(i)))=='As')  m_i=74.9216  !        Arsen   As      33
           if (TRIM(adjustl(atomsorte(i)))=='Se')  m_i=78.96    !Selen   Se      34
           if (TRIM(adjustl(atomsorte(i)))=='Br')  m_i=79.904   !Brom    Br      35
           if (TRIM(adjustl(atomsorte(i)))=='Kr')  m_i=83.8     !Krypton         Kr      36
           if (TRIM(adjustl(atomsorte(i)))=='Rb')  m_i=85.4678  !        Rubidium        Rb      37
           if (TRIM(adjustl(atomsorte(i)))=='Sr')  m_i=87.62    !Strontium       Sr      38
           if (TRIM(adjustl(atomsorte(i)))=='Y')  m_i=88.9059   !        Yttrium         Y       39
           if (TRIM(adjustl(atomsorte(i)))=='Zr')  m_i=91.224   !Zirkonium       Zr      40
           if (TRIM(adjustl(atomsorte(i)))=='Nb')  m_i=92.9064  !        Nobium  Nb      41
           if (TRIM(adjustl(atomsorte(i)))=='Mo')  m_i=95.94   !Molybdän        Mo      42
           if (TRIM(adjustl(atomsorte(i)))=='Tc')  m_i=98      !Technetium      Tc      43
           if (TRIM(adjustl(atomsorte(i)))=='Ru')  m_i=101.07  !Ruthenium       Ru      44
           if (TRIM(adjustl(atomsorte(i)))=='Rh')  m_i=102.9055!        Rhodium         Rh      45
           if (TRIM(adjustl(atomsorte(i)))=='Pd')  m_i=106.42  !Palladium       Pd      46
           if (TRIM(adjustl(atomsorte(i)))=='Ag')  m_i=107.8682!        Silber  Ag      47
           if (TRIM(adjustl(atomsorte(i)))=='Cd')  m_i=112.411 !        Kadmium         Cd      48
           if (TRIM(adjustl(atomsorte(i)))=='In')  m_i=114.818 !        Indium  In      49
           if (TRIM(adjustl(atomsorte(i)))=='Sn')  m_i=118.71  !Zinn    Sn      50
           if (TRIM(adjustl(atomsorte(i)))=='Sb')  m_i=121.76  !Antimon         Sb      51
           if (TRIM(adjustl(atomsorte(i)))=='I')  m_i=126.9045!        Jod     I       53
           if (TRIM(adjustl(atomsorte(i)))=='Te')  m_i=127.6   !Tellur  Te      52
           if (TRIM(adjustl(atomsorte(i)))=='Xe')  m_i=131.293 !        Xenon   Xe      54
           if (TRIM(adjustl(atomsorte(i)))=='Re')  m_i=186.207 !        Rhenium         Re      75
           if (TRIM(adjustl(atomsorte(i)))=='Os')  m_i=190.23  !Osmium  Os      76
           if (TRIM(adjustl(atomsorte(i)))=='Ir')  m_i=192.217 !        Iridium         Ir      77
           if (TRIM(adjustl(atomsorte(i)))=='Pt')  m_i=195.078 !        Platin  Pt      78
           if (TRIM(adjustl(atomsorte(i)))=='Au')  m_i=196.9665!        Gold    Au      79
           if (TRIM(adjustl(atomsorte(i)))=='Hg')  m_i=200.59  !Quecksilber     Hg      80
           if (TRIM(adjustl(atomsorte(i)))=='Tl')  m_i=204.3833!        Thallium        Tl      81
           if (TRIM(adjustl(atomsorte(i)))=='Pb')  m_i=207.2   !Blei    Pb      82
           if (TRIM(adjustl(atomsorte(i)))=='Bi')  m_i=208.9804!        Wismut  Bi      83
           if (TRIM(adjustl(atomsorte(i)))=='Po')  m_i=209     !Polonium        Po      84
           if (TRIM(adjustl(atomsorte(i)))=='At')  m_i=210 !Astat   At      85
           if (TRIM(adjustl(atomsorte(i)))=='Rn')  m_i=222 !    Radon
           
           schwerpunkt(:)=schwerpunkt(:)+m_i*koord(i,:)     
           m_ges=m_ges+m_i     
           !write(*,*) m_i,schwerpunkt
        end do 
        schwerpunkt=schwerpunkt/real(m_ges)
        ! Berechnung der Anzahl der Elektronen im neutralen Molekül z_ges 
END SUBROUTINE Massenschwerpunkt 
       
       
SUBROUTINE N_electrons_in_molecule(atomsorte,MAXATM,z_ges) ! 
! subroutine berechnet die Anzahl der Elektronen z_ges im gesamten molekuel      
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: z_ges
        INTEGER, INTENT(IN) :: MAXATM
        Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte       
        integer::i,z_i
        
        z_ges=0
        do i=1,MAXATM
           z_i=0
           ! .gro Gruppen
           if (TRIM(adjustl(atomsorte(i)))=='HC')  z_i=1    ! fuer H
           if (TRIM(adjustl(atomsorte(i)))=='CH1') z_i=12   ! fuer 1H und 1C 
           if (TRIM(adjustl(atomsorte(i)))=='CH2') z_i=14   ! fuer 2H und 1C  
           if (TRIM(adjustl(atomsorte(i)))=='CH3') z_i=15   ! fuer 3H und 1C  
           if (TRIM(adjustl(atomsorte(i)))=='CA')  z_i=12   ! fuer CA = C  
           ! Normale Atome
           if (TRIM(adjustl(atomsorte(i)))=='H')  z_i=1!Wasserstoff     H       1
           if (TRIM(adjustl(atomsorte(i)))=='He') z_i=2!Helium          He      2
           if (TRIM(adjustl(atomsorte(i)))=='Li') z_i=3!Lithium         Li      3
           if (TRIM(adjustl(atomsorte(i)))=='Be') z_i=4!Beryllium       Be      4
           if (TRIM(adjustl(atomsorte(i)))=='B')  z_i=5!Bor     B       5
           if (TRIM(adjustl(atomsorte(i)))=='C')  z_i=6!Kohlenstoff     C       6
           if (TRIM(adjustl(atomsorte(i)))=='N')  z_i=7!        Stickstoff      N       7
           if (TRIM(adjustl(atomsorte(i)))=='O')  z_i=8!        Sauerstoff      O       8
           if (TRIM(adjustl(atomsorte(i)))=='F')  z_i=9!        Fluor   F       9
           if (TRIM(adjustl(atomsorte(i)))=='Ne') z_i=10!        Neon    Ne      10
           if (TRIM(adjustl(atomsorte(i)))=='Na') z_i=11!        Natrium         Na      11
           if (TRIM(adjustl(atomsorte(i)))=='Mg') z_i=12!Magnesium       Mg      12
           if (TRIM(adjustl(atomsorte(i)))=='Al') z_i=13!        Aluminum        Al      13
           if (TRIM(adjustl(atomsorte(i)))=='Si') z_i=14!        Silizium        Si      14
           if (TRIM(adjustl(atomsorte(i)))=='P')  z_i=15!        Phosphor        P       15
           if (TRIM(adjustl(atomsorte(i)))=='S')  z_i=16!Schwefel        S       16
           if (TRIM(adjustl(atomsorte(i)))=='Cl') z_i=17!Chlor   Cl      17
           if (TRIM(adjustl(atomsorte(i)))=='K')  z_i=19!        Kalium  K       19
           if (TRIM(adjustl(atomsorte(i)))=='Ar') z_i=17!Argon   Ar      18
           if (TRIM(adjustl(atomsorte(i)))=='Ca') z_i=20!Kalzium         Ca      20
           if (TRIM(adjustl(atomsorte(i)))=='Sc') z_i=21!        Scandium        Sc      21
           if (TRIM(adjustl(atomsorte(i)))=='Ti') z_i=22!Titan   Ti      22
           if (TRIM(adjustl(atomsorte(i)))=='V')  z_i=23 !        Vanadium        V       23
           if (TRIM(adjustl(atomsorte(i)))=='Cr') z_i=24!        Chrom   Cr      24
           if (TRIM(adjustl(atomsorte(i)))=='Mn') z_i=25!Mangan  Mn      25
           if (TRIM(adjustl(atomsorte(i)))=='Fe') z_i=26!Eisen   Fe      26
           if (TRIM(adjustl(atomsorte(i)))=='Ni') z_i=28!        Nickel  Ni      28
           if (TRIM(adjustl(atomsorte(i)))=='Co') z_i=27!        Kobalt  Co      27
           if (TRIM(adjustl(atomsorte(i)))=='Cu') z_i=29!Kupfer  Cu      29
           if (TRIM(adjustl(atomsorte(i)))=='Zn') z_i=30!Zink    Zn      30
           if (TRIM(adjustl(atomsorte(i)))=='Ga') z_i=31!Gallium         Ga      31
           if (TRIM(adjustl(atomsorte(i)))=='Ge') z_i=32!Germanium       Ge      32
           if (TRIM(adjustl(atomsorte(i)))=='As') z_i=33!        Arsen   As      33
           if (TRIM(adjustl(atomsorte(i)))=='Se') z_i=34!Selen   Se      34
           if (TRIM(adjustl(atomsorte(i)))=='Br') z_i=35!Brom    Br      35
           if (TRIM(adjustl(atomsorte(i)))=='Kr') z_i=36!Krypton         Kr      36
           if (TRIM(adjustl(atomsorte(i)))=='Rb') z_i=37!        Rubidium        Rb      37
           if (TRIM(adjustl(atomsorte(i)))=='Sr') z_i=38!Strontium       Sr      38
           if (TRIM(adjustl(atomsorte(i)))=='Y')  z_i=39!        Yttrium         Y       39
           if (TRIM(adjustl(atomsorte(i)))=='Zr') z_i=40!Zirkonium       Zr      40
           if (TRIM(adjustl(atomsorte(i)))=='Nb') z_i=41!        Nobium  Nb      41
           if (TRIM(adjustl(atomsorte(i)))=='Mo') z_i=42!Molybdän        Mo      42
           if (TRIM(adjustl(atomsorte(i)))=='Tc') z_i=43!Technetium      Tc      43
           if (TRIM(adjustl(atomsorte(i)))=='Ru') z_i=44!Ruthenium       Ru      44
           if (TRIM(adjustl(atomsorte(i)))=='Rh') z_i=45!        Rhodium         Rh      45
           if (TRIM(adjustl(atomsorte(i)))=='Pd') z_i=46!Palladium       Pd      46
           if (TRIM(adjustl(atomsorte(i)))=='Ag') z_i=47!        Silber  Ag      47
           if (TRIM(adjustl(atomsorte(i)))=='Cd') z_i=48!        Kadmium         Cd      48
           if (TRIM(adjustl(atomsorte(i)))=='In') z_i=49!        Indium  In      49
           if (TRIM(adjustl(atomsorte(i)))=='Sn') z_i=50!Zinn    Sn      50
           if (TRIM(adjustl(atomsorte(i)))=='Sb') z_i=51!Antimon         Sb      51
           if (TRIM(adjustl(atomsorte(i)))=='I')  z_i=53!        Jod     I       53
           if (TRIM(adjustl(atomsorte(i)))=='Te') z_i=52!Tellur  Te      52
           if (TRIM(adjustl(atomsorte(i)))=='Xe') z_i=54!        Xenon   Xe      54
           if (TRIM(adjustl(atomsorte(i)))=='Re') z_i=75!        Rhenium         Re      75
           if (TRIM(adjustl(atomsorte(i)))=='Os') z_i=76!Osmium  Os      76
           if (TRIM(adjustl(atomsorte(i)))=='Ir') z_i=77!        Iridium         Ir      77
           if (TRIM(adjustl(atomsorte(i)))=='Pt') z_i=78!        Platin  Pt      78
           if (TRIM(adjustl(atomsorte(i)))=='Au') z_i=79!        Gold    Au      79
           if (TRIM(adjustl(atomsorte(i)))=='Hg') z_i=80!Quecksilber     Hg      80
           if (TRIM(adjustl(atomsorte(i)))=='Tl') z_i=81!        Thallium        Tl      81
           if (TRIM(adjustl(atomsorte(i)))=='Pb') z_i=82!Blei    Pb      82
           if (TRIM(adjustl(atomsorte(i)))=='Bi') z_i=83!        Wismut  Bi      83
           if (TRIM(adjustl(atomsorte(i)))=='Po') z_i=84    !Polonium        Po      84
           if (TRIM(adjustl(atomsorte(i)))=='At') z_i=85     !Astat   At      85
           if (TRIM(adjustl(atomsorte(i)))=='Rn') z_i=86      !    Radon
           
           !if (z_i==0) THEN
            !write(*,*) 'Es ist ein Fehler bei der Berechnung der Elektronenanzahl z_ges aufgetreten!'
           !END IF
        z_ges=z_ges+z_i
        END DO
        
        !WRITE(*,*) 'Die Anzahl der Elekronen im neutralen Molekuel ist',z_ges
        !write(*,*) 'Der Schwerpunkt liegt bei ', schwerpunkt
END SUBROUTINE N_electrons_in_molecule







SUBROUTINE make_coord_to_xyz(atomsorte,koord,MAXATM,ziel)
!! Subroutine um aus den Inputkoordinaten eine xyz-Datei zu erstellen.
IMPLICIT NONE
Integer, INTENT (IN):: MAXATM
Real, Dimension(MAXATM,3), INTENT (IN) :: koord
Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
Character(500), INTENT (IN) :: ziel
Integer :: i,ierror

OPEN(UNIT=77,FILE=TRIM(ziel)//'.xyz',STATUS='REPLACE',IOSTAT=ierror)
Write(77,'(I8,/,"  genarated by ")') MAXATM
DO i= 1, MAXATM, 1
    write(77,'(2X,A5,1X,3F16.6)') atomsorte(i),koord(i,1),koord(i,2),koord(i,3)
END DO
ClOSE(77)
WRITE(*,*) 'Koordinaten in ',TRIM(ziel),'.xyz geschrieben'

END SUBROUTINE make_coord_to_xyz

SUBROUTINE make_coord_to_xyz_in_dir(atomsorte,koord,MAXATM,ziel,foldername)
!! Subroutine um aus den Inputkoordinaten eine xyz-Datei in einem Unterordner "foldername" zu erstellen.
    IMPLICIT NONE
    Integer, INTENT (IN):: MAXATM
    Real, Dimension(MAXATM,3), INTENT (IN) :: koord
    Character(5), Dimension(MAXATM), INTENT(IN) :: atomsorte
    Character(500), INTENT (IN) :: ziel
    Character(100) ::foldername
     Character(500) ::bashline
    Integer :: i,ierror
    LOGICAL::folder_vorhanden
    
    inquire(file=TRIM(foldername),exist=folder_vorhanden) ! Erstelle Ordner, wenn er noch nicht existiert
    if (.NOT. folder_vorhanden) THEN
            write(bashline,*) 'mkdir '//TRIM(foldername)
            CALL execute_command_line(TRIM(bashline)) 
            WRITE(*,*) TRIM(bashline)
    END IF
    
    OPEN(UNIT=77,FILE=TRIM(foldername)//'/'//TRIM(ziel)//'.xyz',STATUS='REPLACE',IOSTAT=ierror)
    Write(77,'(I8,/,"  genarated by ")') MAXATM
    DO i= 1, MAXATM, 1
        write(77,'(2X,A3,1X,3F16.6)') atomsorte(i),koord(i,1),koord(i,2),koord(i,3)
    END DO
    ClOSE(77)
    WRITE(*,*) 'Koordinaten in ',TRIM(foldername)//'/'//TRIM(ziel)//'.xyz geschrieben'

END SUBROUTINE make_coord_to_xyz_in_dir


SUBROUTINE distance_PBC_correction(dist,x_box,y_box,z_box)
	!! PBC Correction of distance vector, make sure the input vector and the box vector is in the same length units.
	IMPLICIT NONE
	REAL, INTENT(IN)::x_box,y_box,z_box
        REAL, DIMENSION (3),INTENT(INOUT)::dist
	REAL, DIMENSION (3)::box_size
	INTEGER::i
    	box_size(1)=x_box
    	box_size(2)=y_box
    	box_size(3)=z_box

    	DO i=1,3
        	IF( dist(i) > 0.5*box_size(i) ) THEN
        	        dist(i) = dist(i)-box_size(i)
        	ELSE IF ( dist(i) < -0.5*box_size(i) ) THEN
        	        dist(i) = dist(i)+box_size(i)
        	ELSE
        	        CYCLE
        	END IF
   	END DO

END SUBROUTINE distance_PBC_correction


FUNCTION countsubstring(s1, s2) result(c)
	!!  counts number of substrings s2 in a string s1
        !!  n = countsubstring("the three truths", "th")
 	!!  => n=3
	  character(*), intent(in) :: s1, s2
	  integer :: c, p, posn
	 
	  c = 0
	  if(len(s2) == 0) return
	  p = 1
	  do 
	    posn = index(s1(p:), s2)
	    if(posn == 0) return
	    c = c + 1
	    p = p + posn + len(s2)
	  end do
END FUNCTION countsubstring




!!!!!!!!!!!!! START SUBROUTINES FOR PPDI / PBDT-TS1!!!!!!!!!!
INTEGER FUNCTION PBDT_TS1_sequence_to_index(sequence)
	!! index to access the reference lambda value
	CHARACTER(LEN=50), INTENT(IN)::sequence
	LOGICAL::DEBUG=.false.
	INTEGER::i
	
	IF( TRIM(adjustl(sequence)) == 'A') THEN
		i=2
	ELSE IF ( TRIM(adjustl(sequence)) == 'B')THEN
		i=3
	ELSE IF ( (TRIM(adjustl(sequence)) == 'AB') .OR. (TRIM(adjustl(sequence)) == 'BA') )THEN
		i=4
	ELSE IF ( TRIM(adjustl(sequence)) == 'ABA')THEN
		i=5
	ELSE IF ( TRIM(adjustl(sequence)) == 'BAB')THEN
		i=6		
	ELSE IF ( (TRIM(adjustl(sequence)) == 'ABAB') .OR. (TRIM(adjustl(sequence)) == 'BABA') )THEN
		i=7	
	ELSE IF ( TRIM(adjustl(sequence)) == 'ABABA')THEN
		i=8		
	ELSE IF ( TRIM(adjustl(sequence)) == 'BABAB')THEN
		i=9		
	ELSE IF ( (TRIM(adjustl(sequence)) == 'ABABAB') .OR. (TRIM(adjustl(sequence)) == 'BABABA') )THEN
		i=10	
	ELSE IF ( TRIM(adjustl(sequence)) == 'ABABABA')THEN
		i=11		
	ELSE IF ( TRIM(adjustl(sequence)) == 'BABABAB')THEN
		i=12			
	ELSE IF ( (TRIM(adjustl(sequence)) == 'ABABABAB') .OR. (TRIM(adjustl(sequence)) == 'BABABABA') )THEN
		i=13		
	ELSE IF ( TRIM(adjustl(sequence)) == 'ABABABABA')THEN
		i=14		
	ELSE IF ( TRIM(adjustl(sequence)) == 'BABABABAB')THEN
		i=14
	ELSE IF ( (TRIM(adjustl(sequence)) == 'ABABABABAB') .OR. (TRIM(adjustl(sequence)) == 'BABABABABA') )THEN
		i=15
	ELSE IF ( TRIM(adjustl(sequence)) == 'ABABABABABA')THEN
		i=15		
	ELSE IF ( TRIM(adjustl(sequence)) == 'BABABABABAB')THEN
		i=15
	ELSE IF ( (TRIM(adjustl(sequence)) == 'ABABABABABAB') .OR. (TRIM(adjustl(sequence)) == 'BABABABABABA') )THEN ! 6mer
		i=16
	ELSE IF ( TRIM(adjustl(sequence)) == 'ABABABABABABA')THEN ! 6mer
		i=16		
	ELSE IF ( TRIM(adjustl(sequence)) == 'BABABABABABAB')THEN ! 7mer
		i=17	
	ELSE IF ( (TRIM(adjustl(sequence)) == 'ABABABABABABAB') .OR. (TRIM(adjustl(sequence)) == 'BABABABABABABA') )THEN ! 7mer
		i=17
	ELSE IF ( TRIM(adjustl(sequence)) == 'ABABABABABABABA')THEN ! 7mer
		i=17		
	ELSE IF ( TRIM(adjustl(sequence)) == 'BABABABABABABAB')THEN ! 8mer
		i=17
	ELSE IF ( (TRIM(adjustl(sequence)) == 'ABABABABABABABAB') .OR. (TRIM(adjustl(sequence)) == 'BABABABABABABABA') )THEN ! 8mer
		i=17
	ELSE IF ( TRIM(adjustl(sequence)) == 'ABABABABABABABABA')THEN								  
		i=17		
	ELSE IF ( TRIM(adjustl(sequence)) == 'BABABABABABABABAB')THEN
		i=17					
	ELSE IF ( TRIM(adjustl(sequence)) == 'PPDI' .OR. TRIM(adjustl(sequence)) == 'PPDI_') THEN
		i=1
	ELSE
		WRITE(*,*) 'Error: string input to PBDT_TS1_sequence_to_index was not determined: '//TRIM(sequence)
		CALL EXIT (1)
	END IF
	IF (DEBUG) WRITE(*,*) 'PBDT_TS1_sequence_to_index sequence:   ',TRIM(sequence),' index:',i	
	PBDT_TS1_sequence_to_index=i
	
END FUNCTION PBDT_TS1_sequence_to_index



INTEGER FUNCTION seq_to_k_ind(sequence,fin)
	!! Subroutine to determine the k index, when given the polymer sequence of PDBT-TS1 /PPDI 
	!! input  :: sequence = 'ABABAB', fin='M' (middle segment)
	!! output :: k = 7
	CHARACTER(LEN=20), INTENT(IN) ::  sequence
	CHARACTER(LEN=1), INTENT(IN)  :: fin
	INTEGER::k,N_A,N_B
	k=0
	N_A=countsubstring(sequence, 'A')
	N_B=countsubstring(sequence, 'B')
	
	IF( INDEX(sequence,'PPDI') /= 0) THEN
		k=1
	ELSE IF (fin == 'M') THEN
		IF( sequence(1:1) == 'A') THEN
			k=N_A+N_B+1
		ELSE 
			k=N_A+N_B+16
		END IF
	ELSE IF (fin == 'S') THEN
		IF( sequence(1:1) == 'A') THEN
			k=N_A+N_B+31
		ELSE 
			k=N_A+N_B+46
		END IF
	ELSE IF (fin == 'E') THEN
		IF( sequence(1:1) == 'A') THEN
			k=N_A+N_B+61
		ELSE 
			k=N_A+N_B+76
		END IF
	ELSE IF ( (N_A==8) .AND. (N_B==8))THEN
		k=92
	END IF
	seq_to_k_ind=k
END FUNCTION seq_to_k_ind



!!!!!!!!!!!!! END SUBROUTINES FOR PPDI / PBDT-TS1!!!!!!!!!!




!!!!!!!!!!!!!!! START make_min_d_COM_vector_list !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine which  calculates the center of masses for subhopping sites on DIPBI, P3HT, PPDI, PBDT_TS1
!!! It also sets a new data structure for system_d usind a module
SUBROUTINE make_min_d_COM_vector_list(system,N_Resids,R_u_I,MAXATM,residue_number,residue_names,atom,koordA,&
				  &neighbourlist_filename,gro_inputfile,N_DIPBI,N_P3HT,N_PPDI,N_PBDT_TS1,N_Species,N_Neighbours,neighbour_list_Res,&
				  &x_box,y_box,z_box)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: MAXATM,N_Resids
	INTEGER, DIMENSION(N_Resids+1,2), INTENT(IN)             :: R_u_I
	REAL,        DIMENSION(MAXATM,3), INTENT(IN)             :: koordA
	CHARACTER(LEN=5), DIMENSION(MAXATM), INTENT(IN)          :: residue_names,atom
	INTEGER         , DIMENSION(MAXATM), INTENT(IN)          :: residue_number
	INTEGER, INTENT(IN) ::N_DIPBI,N_P3HT,N_PPDI,N_PBDT_TS1
	CHARACTER(LEN=*), INTENT(IN)                             ::gro_inputfile,neighbourlist_filename
	INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT)        ::neighbour_list_Res
	REAL , INTENT(IN)   ::x_box,y_box,z_box
	INTEGER, INTENT(OUT)  ::N_Neighbours
	!!!!!!!!!!!! contains all system information after initialization
	type (system_t) , INTENT(INOUT):: system
	!! local array with minimum distance list
	REAL, ALLOCATABLE, DIMENSION (:,:,:) :: min_d_COM_arr
	!!! read d_COM of neighbours
	REAL, ALLOCATABLE, DIMENSION(:,:)::neighbour_list_COM
	!! DIPBI
	CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:)  :: DIPBI_atom
	REAL, ALLOCATABLE, DIMENSION(:,:)            :: DIPBI_koord,DIPBI_COM
	!! P3HT
	REAL, DIMENSION(6,3) ::Ringe_atome
	CHARACTER(5), DIMENSION(6) :: Ringe_atomsorte
	REAL, ALLOCATABLE, DIMENSION(:,:)::P3HT_COM
	!! PPDI
	CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:)      :: PPDI_atom
	REAL, ALLOCATABLE, DIMENSION(:,:)    		 :: PPDI_koord,PPDI_COM	
	!! PBDT_TS1
	CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:)      :: PBDT_TS1_atom,CarboxyS2_atom,PBDT_atom
	REAL, ALLOCATABLE, DIMENSION(:,:)    		 :: PBDT_TS1_koord,CarboxyS2_koord,PBDT_koord,PBDT_TS1_COM		
	INTEGER, ALLOCATABLE, DIMENSION(:)	:: ind_CarboxyS2,ind_PBDT
	!! Other species, all atoms are inclueded in COM determination
	INTEGER::N_Species,m_Spec
	CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:)  :: Spec_atom
	REAL, ALLOCATABLE, DIMENSION(:,:)    		 :: Spec_koord,Spec_COM	
	!! variables to determine nearest neighbour vector between (sub)-hopping sites.
	REAL, DIMENSION (3) ::dist_vec, dist_vec_min
	REAL::norm_dist_vec, norm_dist_vec_min
	INTEGER::molAid,molBid,ResidA,ResidB,i_offset
		
	CHARACTER(LEN=500)  ::COM_filename
	REAL, DIMENSION(3)  ::COM,COM_A,COM_B
	INTEGER::N_Schwefel_all,m_CarboxyS2,m_PBDT,n_CarboxyS2,n_PBDT
	INTEGER::i,j,k,k_max,l,l_max,m,ierror,i_resid,i_COM,ind,N_mini_sites
	LOGICAL::gro_format,DEBUG=.false.,print_COM_molecules_to_xyz=.false.

	
IF(DEBUG)	WRITE(*,*) 'start:   make_min_d_COM_vector_list: N_Neighbours: ',N_Neighbours
	! Einlesen der Daten aus der Nachbarschaftsliste
IF( .NOT. allocated(neighbour_list_Res)) THEN
		ALLOCATE(neighbour_list_COM(N_Resids,N_Neighbours+1))
		ALLOCATE(neighbour_list_Res(N_Resids,N_Neighbours+1))
		CALL Read_neighbourlist(neighbourlist_filename,N_Resids,N_Neighbours,neighbour_list_COM,neighbour_list_Res)
END IF
N_Schwefel_all=0
DO j=1,MAXATM
               IF (TRIM(adjustl(atom(j)))=='S') 	N_Schwefel_all=N_Schwefel_all+1
               IF (TRIM(adjustl(atom(j)))=='S-S_R') N_Schwefel_all=N_Schwefel_all+1
END DO

IF (N_P3HT /= 0)  THEN
	N_mini_sites=N_Schwefel_all
ELSE IF (N_PBDT_TS1 /= 0) THEN
	N_mini_sites=464*(8+8)
END IF

WRITE(*,*) N_mini_sites

IF( (N_mini_sites /= 416*32) .AND. (N_P3HT /= 0) ) THEN
	WRITE(*,*) '# S atoms: ',N_Schwefel_all
	WRITE(*,*) 'Only for 416 P3HT in the box!'
	CALL EXIT(1)
ELSE IF( (N_mini_sites /= 464*(8+8)) .AND. (N_PBDT_TS1 /= 0) ) THEN
	WRITE(*,*) '# S atoms: ',N_Schwefel_all
	WRITE(*,*) 'Only for 464 PBDT_TS1 in the box!'
	CALL EXIT(1)
ELSE
	WRITE(*,*) ' N_DIPBI: ',N_DIPBI,' N_P3HT: ',N_P3HT,' N_PPDI: ',N_PPDI,' N_PBDT_TS1: ',N_PBDT_TS1,' N_Species:',N_Species
END IF



i_resid=0
!! Get COM for centers nur fur DIPBI
IF( N_DIPBI /= 0 ) THEN
		!! Center of Masses in DIPBI subunits
		ALLOCATE(DIPBI_COM(N_DIPBI,3))
		ALLOCATE(DIPBI_koord(74,3))
		ALLOCATE(DIPBI_atom(74))
                
		i_COM=0   ! index fuer die Anzahl der Schwerpunkt von DIPBI
		DO i=1,MAXATM
			IF(  ((i == 1) .OR. (residue_number(i) /= residue_number(i-1))) .AND. & !New resid 
			&	  ((TRIM(adjustl(residue_names(i)))=='DIPBI') .OR. (TRIM(adjustl(residue_names(i)))=='DPBIK') &
			& .OR. (TRIM(adjustl(residue_names(i)))=='dipbi')) )THEN ! New DIPBI			
					i_COM=i_COM+1
					DO k=1,74
									DIPBI_atom(k)=atom(i+k-1)
                        	        DIPBI_koord(k,1:3)=koordA(i+k-1,1:3)
                    END DO !k
					CALL Massenschwerpunkt(DIPBI_atom,DIPBI_koord,74,COM,.false.)
                    DIPBI_COM(i_COM,:)=COM(:)
                    
                    !!!! ADD data to system
					!!!                 system,residA,NAtomsA,resnameA,molnameA,koordA,elementsA
					i_resid=i_resid+1
					call make_new_molecule_for_system(system,R_u_I(i_resid,1),R_u_I(i_resid+1,2)-R_u_I(i_resid,2),residue_names(R_u_I(i_resid,2))&
				&,residue_names(R_u_I(i_resid,2)),koordA(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1,1:3),atom(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1))
					call set_COM_in_Molecule(system%molecules(i_resid),COM(1),COM(2),COM(3))
					call set_COM_in_Hoppingsite(system%molecules(i_resid)%hoppingsite,COM(1),COM(2),COM(3))
					
			END IF ! new DIPBI ?
		END DO ! loop all atoms
				
		IF (i_COM /= N_DIPBI) THEN
			WRITE(*,*) 'Error: Missmatch in number of evaluated DIPBI: ',i_COM,' N_DIPBI:',N_DIPBI
			CALL EXIT(1)
		END IF		
END IF !N_DIPBI



!! Get COM for centers nur fur P3HT
IF( N_P3HT /= 0 ) THEN
				!! Center of Masses in thiophene subunits
				IF(DEBUG) WRITE(*,*) 'START P3HT'
                ALLOCATE(P3HT_COM(N_mini_sites,3))
                ind=0 		! index fur die Atome in den Ringen von Mol A
				i_COM=0   ! index fuer die Anzahl der Schwerpunkt von A
                DO i=1,MAXATM

						IF ( (TRIM(adjustl(residue_names(i)(1:2)))=='PH') .OR. (TRIM(adjustl(residue_names(i)(1:2)))=='PM' ) .OR. &
						&    (TRIM(adjustl(residue_names(i)(1:2)))=='P3') .OR. (TRIM(adjustl(residue_names(i)(1:3)))=='THP') ) THEN  !!! FIND P3HT/P3MT/THP
							IF( (i == 1) .OR. (residue_number(i) /= residue_number(i-1)) ) THEN ! New segment AND P3HT segment
									!WRITE(*,*) 'TMP:',residue_number(i),residue_names(i)
									!!!! ADD data to system
									!!!                 system,residA,NAtomsA,resnameA,molnameA,koordA,elementsA
									i_resid=i_resid+1
									call make_new_molecule_for_system(system,R_u_I(i_resid,1),R_u_I(i_resid+1,2)-R_u_I(i_resid,2),&
									&residue_names(R_u_I(i_resid,2)),residue_names(R_u_I(i_resid,2)),&
									&koordA(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1,1:3),atom(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1))							 
							END IF !!! add new segment
                
							!!! Evaluate sub-unit Center of Masses and add them to sub-Hopping-sites
							IF( (TRIM(adjustl(atom(i)))=='S') .OR. (TRIM(adjustl(atom(i)))=='S-S_R') )THEN
									DO j=(i-5),i,1
											ind=ind+1
											Ringe_atomsorte(ind)=atom(j)
											Ringe_atome(ind,:)=koordA(j,:)
											!WRITE(*,*) ind,'Ringsystem: ',Ringe_atomsorte(ind),Ringe_atome(ind,:)
									END DO
									! Abspeichern des COM des Ringes
									IF (ind == 6) THEN
											! Korrektur, wenn das Ende der THP Kette erreicht ist, da hier die Reihung geaendert ist,
											! ersetze End HC durch CA in der Liste 
											IF( (TRIM(adjustl(atom(i-1)))=='H') .OR. (TRIM(adjustl(atom(i-1)))=='HC') .OR.&
													& (TRIM(adjustl(atom(i-1)))=='H-H_') ) THEN
													Ringe_atomsorte(5)=atom(j-7)
													Ringe_atome(5,:)=koordA(j-7,:)
													!WRITE(*,*) 'Ersetzt CA:', Ringe_atomsorte(5),Ringe_atome(5,:)
											END IF

											ind=0
											i_COM=i_COM+1
											CALL Massenschwerpunkt(Ringe_atomsorte,Ringe_atome,6,COM,.false.)
											P3HT_COM(i_COM,:)=COM(:)
	  
											!!! Add subhoppingsite to molecule
											call make_new_subhoppingsite_for_hoppingsite(system%molecules(i_resid)%hoppingsite,COM)
					  
									END IF ! calc COM_1
							END IF !!! END  Evaluate sub-unit Center of Masses and add them to sub-Hopping-sites     
                        
							IF ( i == R_u_I(i_resid+1,2)-1) THEN !!! Last atom in segment
													!!! Select COM from subsegment COMs
													!!! set middle segment as COM for all with selected id = N_subhoppingsites/2+1
													ind=int(int(system%molecules(i_resid)%hoppingsite%N_subhoppingsites)/2+1)
													COM=system%molecules(i_resid)%hoppingsite%subhoppingsite(ind)%koord(1:3)  
													call set_COM_in_Molecule(system%molecules(i_resid),COM(1),COM(2),COM(3))
													call set_COM_in_Hoppingsite(system%molecules(i_resid)%hoppingsite,COM(1),COM(2),COM(3))      
													ind=0     
													!WRITE(*,*) ' This COM ',i_resid,' was added to the hopping site !', COM               
							END IF      !!! Last atom in segment        
                                                
                        ENDIF !P3HT/P3MT/THP found
                END DO ! i=,1,MAXATM
                         
                
 				IF (i_COM /= N_mini_sites) THEN
					WRITE(*,*) 'Error: Missmatch in number of evaluated P3HT: ',i_COM,' N_mini_sites: ', N_mini_sites
					CALL EXIT(1)
				END IF
          
END IF


!! Get COM for centers nur fur PPDI
IF( N_PPDI /= 0 ) THEN
				!! Center of Masses in PPDI subunits, PPDI ohne Seitenketten!
                ALLOCATE(PPDI_COM(N_PPDI,3))
                ALLOCATE(PPDI_koord(76,3))
                ALLOCATE(PPDI_atom(76))
                
				i_COM=0   ! index fuer die Anzahl der Schwerpunkt von DIPBI
                DO i=1,MAXATM
					IF( (TRIM(adjustl(residue_names(i)(1:4)))=='PPDI') .AND. ( (i == 1) .OR. (residue_number(i) /= residue_number(i-1)) ) )THEN ! New PPDI
						i_COM=i_COM+1
						 
						IF( TRIM(adjustl(atom(i)(1:1))) == 'N' ) THEN !!Koordinaten an neuen Array übergeben nicht die Seitenketten! N=1 wurde vorher schon abgeschnitten 
							ind=-1 ! Starte mit N atom an position 1
						ELSE ! Seitenkette ueberspringen starte mit atom 7
							ind=6
						END IF
						DO k=1,76
									PPDI_atom(k)=atom(i+k+ind)
                                    PPDI_koord(k,1:3)=koordA(i+k+ind,1:3)
                        END DO
						CALL Massenschwerpunkt(PPDI_atom,PPDI_koord,76,COM,.false.)
                        PPDI_COM(i_COM,:)=COM(:)	
                        

						!!!! ADD data to system
						!!!                 system,residA,NAtomsA,resnameA,molnameA,koordA,elementsA
						i_resid=i_resid+1
						call make_new_molecule_for_system(system,R_u_I(i_resid,1),R_u_I(i_resid+1,2)-R_u_I(i_resid,2),&
						&residue_names(R_u_I(i_resid,2)),residue_names(R_u_I(i_resid,2)),&
						&koordA(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1,1:3),atom(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1))
					
						call set_COM_in_Molecule(system%molecules(i_resid),COM(1),COM(2),COM(3))
						call set_COM_in_Hoppingsite(system%molecules(i_resid)%hoppingsite,COM(1),COM(2),COM(3))     
											                   	
					END IF ! new PPDI ?
				END DO ! loop all atoms
				
				IF (i_COM /= N_PPDI) THEN
					WRITE(*,*) 'Error: Missmatch in number of evaluated PPDI'
					CALL EXIT(1)
				END IF		
END IF !N_PPDI



!! Get COM for centers just for PBDT_TS1
!! No H-Atoms adde before.
IF( N_PBDT_TS1 /= 0 ) THEN 
		!! Center of Masses in thiophene subunits
                ALLOCATE(PBDT_TS1_COM(N_mini_sites,3))
                m_CarboxyS2=9
				ALLOCATE(CarboxyS2_koord(m_CarboxyS2,3))
				ALLOCATE(CarboxyS2_atom(m_CarboxyS2))
				ALLOCATE(ind_CarboxyS2(m_CarboxyS2))

                m_PBDT=16
				ALLOCATE(PBDT_koord(m_PBDT,3))
				ALLOCATE(PBDT_atom(m_PBDT))
				ALLOCATE(ind_PBDT(m_PBDT))
 			
				gro_format=.false.
				DO j=1,size(atom)/5	
					IF (TRIM(adjustl(atom(j)))=='CH3') gro_format=.true.
					IF (TRIM(adjustl(atom(j)))=='CH2') gro_format=.true.
					IF (TRIM(adjustl(atom(j)))=='CH1') gro_format=.true.
					IF (TRIM(adjustl(atom(j)))=='HC')  gro_format=.true.
					IF (TRIM(adjustl(atom(j)))=='NR')  gro_format=.true.
					IF (TRIM(adjustl(atom(j)))=='OE')  gro_format=.true.
					IF(gro_format) EXIT
				END DO ! find gro format ?
 			
				IF ( gro_format) THEN
					! gro format
					ind_CarboxyS2= (/ 1,2,3,4,5,12,13,14,15 /)
					ind_PBDT= (/ 1,2,3,4,5, 14,15,16,17, 26,27,28,29,30,31,32 /)

				ELSE ! H atoms added 
					ind_CarboxyS2= (/ 1,2,3,4,5,19,20,21,22 /)
					ind_PBDT= (/ 1,2,3,4,5, 17,18,19,20, 32,33,34, 35,36,37,38/)
				END IF
				
                ind=0 		! index fur die Atome in den Ringen von Mol A
				i_COM=0     ! index fuer die Anzahl der Schwerpunkt von A
				DO l=1,N_Resids
					DO i=R_u_I(l,2),R_u_I(l+1,2)-1  ! atom indizes in single resids
							IF( (.NOT. gro_format) .AND. &
							&( (TRIM(adjustl(residue_names(i)(1:1)))=='A') .OR. (TRIM(adjustl(residue_names(i)(1:1)))=='B') )) THEN
										! NEUES SEGMENT POSITION 
										IF(i==R_u_I(l,2)) THEN
												!!!! ADD data to system
												!!!                 system,residA,NAtomsA,resnameA,molnameA,koordA,elementsA
												i_resid=i_resid+1
												call make_new_molecule_for_system(system,R_u_I(i_resid,1),R_u_I(i_resid+1,2)-R_u_I(i_resid,2),&
													&residue_names(R_u_I(i_resid,2)),residue_names(R_u_I(i_resid,2)),&
													&koordA(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1,1:3),atom(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1))  										
												
												!!! START and determine subsegment COMs
												READ(residue_names(i)(2:2),'(I10)') n_CarboxyS2   
												READ(residue_names(i)(4:4),'(I10)') n_PBDT
												ind=0 
												DO WHILE ( n_CarboxyS2 + n_PBDT > 0)
												
														IF( ind /= 0 .OR. (TRIM(adjustl(residue_names(i)(1:1)))=='A') ) THEN 
															! segment A
															DO k=1,m_CarboxyS2
																!! select the index j = i_startindex_resid   +  ind_start_index_sub_segment  + ind_atom_in_CarboxyS2 inside the ring system selected via  ind_CarboxyS2(k) , where m_CarboxyS2 is the number of atom in the selecte ring-subsystem 
																j = i+ ind + ind_CarboxyS2(k)
																CarboxyS2_koord(k,1:3)=koordA(j,1:3) 
																CarboxyS2_atom(k)=atom(j)
														
															END DO !k
															ind=ind+ind_CarboxyS2(m_CarboxyS2)
															
															i_COM=i_COM+1
															CALL Massenschwerpunkt(CarboxyS2_atom,CarboxyS2_koord,m_CarboxyS2,COM,.false.)
															PBDT_TS1_COM(i_COM,:)=COM(:)
															IF(DEBUG) WRITE(*,*) 'COM_A',COM
															!!! Add subhoppingsite to molecule
															call make_new_subhoppingsite_for_hoppingsite(system%molecules(i_resid)%hoppingsite,COM)
															
															n_CarboxyS2=n_CarboxyS2-1
														END IF ! find segment A
														
														IF(n_CarboxyS2==0 .AND. n_PBDT==0 ) EXIT
														IF(DEBUG) WRITE(*,*) 'START PBDT'
														! segment B
														DO k=1,m_PBDT
																j= i + ind + ind_PBDT(k)
																PBDT_koord(k,1:3)=koordA(j,1:3)
																PBDT_atom(k)=atom(j)
														END DO !k
														ind=ind+ind_PBDT(m_PBDT)
																							
														i_COM=i_COM+1
														CALL Massenschwerpunkt(PBDT_atom,PBDT_koord,m_PBDT,COM,.false.)
														PBDT_TS1_COM(i_COM,:)=COM(:)
														IF(DEBUG) WRITE(*,*) 'COM_B',COM
														!!! Add subhoppingsite to molecule
														call make_new_subhoppingsite_for_hoppingsite(system%molecules(i_resid)%hoppingsite,COM)														
														
														n_PBDT=n_PBDT-1
													
												END DO ! WHILE loop all mini segments in Resid
												
												!!! Select COM from subsegment COMs
												!!! set middle segment as COM for all with selected id = N_subhoppingsites/2+1
												ind=int(int(system%molecules(i_resid)%hoppingsite%N_subhoppingsites)/2+1)
												COM=system%molecules(i_resid)%hoppingsite%subhoppingsite(ind)%koord(1:3)   
												call set_COM_in_Molecule(system%molecules(i_resid),COM(1),COM(2),COM(3))
												call set_COM_in_Hoppingsite(system%molecules(i_resid)%hoppingsite,COM(1),COM(2),COM(3))   
																																			
												IF(DEBUG) THEN
													WRITE(*,*) 'N_subhoppingsites: ',system%molecules(i_resid)%hoppingsite%N_subhoppingsites						
													call print_molecule_xyz(system%molecules(i_resid))
													WRITE(*,*) 'Final Selected Hopping site: Ag',system%molecules(i_resid)%hoppingsite%COM(1:3)
												END IF !DEBUG
															
									 END IF  !i==R_u_I(l,2) segmentstart 	
							
							ELSE IF( (TRIM(adjustl(residue_names(i)(1:5)))=='8poly') ) THEN  !!!! DETERMINATION FOR gromacs data input format.
							
										IF(    ((TRIM(adjustl(atom(i)))=='CH3') .AND. (TRIM(adjustl(atom(i-1)))=='CH3'))  &
											& .AND. ( (i == 1) .OR. (residue_number(i) /= residue_number(i-1)) ) .OR.	  &
											&  ((TRIM(adjustl(atom(i)))=='CH3') .AND. (TRIM(adjustl(atom(i-1)(1:4))) =='PPDI')) )THEN ! New PDBT_TS1 chains start
												
												!!!! ADD data to system
												!!!                 system,residA,NAtomsA,resnameA,molnameA,koordA,elementsA
												i_resid=i_resid+1
												call make_new_molecule_for_system(system,R_u_I(i_resid,1),R_u_I(i_resid+1,2)-R_u_I(i_resid,2),&
													&residue_names(R_u_I(i_resid,2)),residue_names(R_u_I(i_resid,2)),&
													&koordA(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1,1:3),atom(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1))  		
																							
												ind=0    
												DO m=1,8    ! segments
														! segment A
														DO k=1,m_CarboxyS2
															j= i + ind + ind_CarboxyS2(k)
															CarboxyS2_koord(k,1:3)=koordA(j,1:3) 
															CarboxyS2_atom(k)=atom(j)
														END DO 
														ind=ind+ind_CarboxyS2(m_CarboxyS2)
														
														i_COM=i_COM+1
														CALL Massenschwerpunkt(CarboxyS2_atom,CarboxyS2_koord,m_CarboxyS2,COM,.false.)
														PBDT_TS1_COM(i_COM,:)=COM(:)
														!WRITE(*,*) 'COM_A',COM
														!!! Add subhoppingsite to molecule
														call make_new_subhoppingsite_for_hoppingsite(system%molecules(i_resid)%hoppingsite,COM)	
														IF ( j >= R_u_I(i_resid+1,2) ) EXIT													
														! segment B
														DO k=1,m_PBDT
																j= i + ind + ind_PBDT(k)
																PBDT_koord(k,1:3)=koordA(j,1:3)
																PBDT_atom(k)=atom(j)
														END DO !k
														ind=ind+ind_PBDT(m_PBDT)
																							
														i_COM=i_COM+1
														CALL Massenschwerpunkt(PBDT_atom,PBDT_koord,m_PBDT,COM,.false.)
														PBDT_TS1_COM(i_COM,:)=COM(:)
														!WRITE(*,*) 'COM_B',COM
														!!! Add subhoppingsite to molecule
														call make_new_subhoppingsite_for_hoppingsite(system%molecules(i_resid)%hoppingsite,COM)		
														
														IF ( j >= R_u_I(i_resid+1,2) ) EXIT													
											END DO ! m=1,8 segments
											
												!!! Select COM from subsegment COMs
												!!! set middle segment as COM for all with selected id = N_subhoppingsites/2+1
												ind=int(int(system%molecules(i_resid)%hoppingsite%N_subhoppingsites)/2+1)
												COM=system%molecules(i_resid)%hoppingsite%subhoppingsite(ind)%koord(1:3)   
												call set_COM_in_Molecule(system%molecules(i_resid),COM(1),COM(2),COM(3))
												call set_COM_in_Hoppingsite(system%molecules(i_resid)%hoppingsite,COM(1),COM(2),COM(3))   												
										END IF ! New PDBT_TS1 chains start
									CALL EXIT (1)
									END IF ! 8poly
									
							END DO ! i=,1,MAXATM in den Resid abschnitten
							
                   END DO ! l=1,N_Resids       
                
 				IF (i_COM /= N_mini_sites) THEN
					WRITE(*,*) 'Error: Missmatch in number of evaluated PBDT_TS1: ',i_COM,' N_mini_sites: ', N_mini_sites
					CALL EXIT(1)
				END IF         
END IF ! N_PBDT_TS1





IF((N_DIPBI == 0) .AND. (N_PPDI==0) .AND. (N_P3HT==0) .AND. (N_PBDT_TS1==0)) THEN
	N_Species=(size(R_u_I(:,:))-2)/2
	IF( N_Species /= 0 ) THEN 
					ALLOCATE(Spec_COM(N_Species,3))
					i_COM=0     ! index fuer die Anzahl der Schwerpunkt von A			    
					DO i=1,N_Species
						m_Spec=R_u_I(i+1,2)-R_u_I(i,2)
						ALLOCATE(Spec_koord(m_Spec,3))
						ALLOCATE(Spec_atom(m_Spec))
						ind=0
						DO k=R_u_I(i,2),R_u_I(i+1,2)-1
							ind=ind+1
							Spec_atom(ind)=atom(k)
							Spec_koord(ind,1:3)=koordA(k,1:3)
						END DO ! index 
						i_COM=i_COM+1
						CALL Massenschwerpunkt(Spec_atom,Spec_koord,m_Spec,COM,.false.)
						Spec_COM(i_COM,:)=COM(:)
						
						!!!! ADD data to system
						!!!                 system,residA,NAtomsA,resnameA,molnameA,koordA,elementsA
						i_resid=i_resid+1
						call make_new_molecule_for_system(system,R_u_I(i_resid,1),R_u_I(i_resid+1,2)-R_u_I(i_resid,2),&
							&residue_names(R_u_I(i_resid,2)),residue_names(R_u_I(i_resid,2)),&
							&koordA(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1,1:3),atom(R_u_I(i_resid,2):R_u_I(i_resid+1,2)-1))
							
						call set_COM_in_Molecule(system%molecules(i_resid),COM(1),COM(2),COM(3))
						call set_COM_in_Hoppingsite(system%molecules(i_resid)%hoppingsite,COM(1),COM(2),COM(3))    						
											
						DEALLOCATE(Spec_koord)
						DEALLOCATE(Spec_atom)
					END DO !  i=1,N_Species
	END IF ! N_Species
END IF ! No other Molecules


IF( print_COM_molecules_to_xyz) THEN
				COM_filename='my_COM_'//TRIM(gro_inputfile(:(INDEX(gro_inputfile, '.gro')-1)))//'.xyz'
				OPEN(UNIT=80,FILE=TRIM(COM_filename),STATUS='REPLACE',IOSTAT=ierror)
				
				!!! COM to *.xyz
				WRITE(80,*) N_DIPBI+N_mini_sites+N_PPDI
				WRITE(80,*) ' '
				DO i=1,N_DIPBI
					WRITE(80,*)  'Se ',DIPBI_COM(i,:)
				END DO 	
				
                IF( N_P3HT /= 0 ) THEN
					DO i=1,N_mini_sites
							WRITE(80,*)  'F ',P3HT_COM(i,:)
					END DO 
				END IF
                DO i=1,N_PPDI
					WRITE(80,*)  'I ',PPDI_COM(i,:)
				END DO 
				IF(N_PBDT_TS1 /= 0 ) THEN
					DO i=1,N_mini_sites
						WRITE(80,*)  'Br ',PBDT_TS1_COM(i,:)
					END DO 
				END IF
END IF

!! Get the distance Vector between two centre of masse, in the reduced description. (e.g. Centre of ringes in polymer chain)
!! i =1,N_Resids,   j=2=Resid1,Nearest neighbour_1,Nearest neighbour_2,Nearest neighbour_3...N_Neighbours+1,      k= koordinates 1= dist_x,2=dist_y,3=dist_z  = COM_Resid2 - COM_Resid1 ;
!! 	min_d_COM_arr(i,1,4) = ResidA  
!!  min_d_COM_arr(i,1,5) = ResidB 
!!  min_d_COM_arr(i,1,6) = norm_dist_vec_min 

WRITE(*,*) 'Start calculation for minimal |d_COM|=COM_B-COM_A nearest neighbours '
ALLOCATE(min_d_COM_arr(N_Resids,N_Neighbours+1,3))

!init_hoppingpair_in_system(this,system)
i_offset=1-neighbour_list_Res(1,1)
min_d_COM_arr=0.0
DO i=1,N_Resids
	ResidA=neighbour_list_Res(i,1) 
	!min_d_COM_arr(i,1,4)=ResidA
	molAid=ResidA+i_offset
	DO j=2,N_Neighbours+1 ! nearest neighbours
		ResidB=neighbour_list_Res(i,j) !NN
		min_d_COM_arr(i,1,5)=ResidB
		molBid=ResidB+i_offset
		! dummy
		dist_vec_min=(/1000.0, 1000.0 ,1000./)
		norm_dist_vec_min=50000
		IF(system%molecules(molAid)%hoppingsite%hasSubhoppingsite)THEN
			k_max=system%molecules(molAid)%hoppingsite%N_subhoppingsites    !k = N_subhoppingsites in A
		ELSE
			k_max=1
		END IF
		
		IF(system%molecules(molBid)%hoppingsite%hasSubhoppingsite)THEN	
			l_max=system%molecules(molBid)%hoppingsite%N_subhoppingsites
		ELSE
			l_max=1
		END IF 
			
		
		DO k=1,k_max	
				!!! get COM for A 
				IF(k_max==1) THEN
					IF(system%molecules(molAid)%hoppingsite%hasCOM) THEN !hoppingsite%hasCOM
						COM_A=system%molecules(molAid)%hoppingsite%COM(1:3)
					ELSE !hoppingsite%hasCOM
						IF(system%molecules(molAid)%hasCOM) THEN
							COM_A=system%molecules(molAid)%COM(1:3)
						ELSE
							!! calc COM if not available for all coordinates
							CALL Massenschwerpunkt(atom(R_u_I(molAid,2):R_u_I(molAid+1,2)-1),&
									&koordA(R_u_I(molAid,2):R_u_I(molAid+1,2)-1,1:3),R_u_I(molAid+1,2)-R_u_I(molAid,2),COM,.false.)
								
							call set_COM_in_Molecule(system%molecules(molAid),COM(1),COM(2),COM(3))
							call set_COM_in_Hoppingsite(system%molecules(molAid)%hoppingsite,COM(1),COM(2),COM(3)) 
							COM_A=COM   								
						ENDIF !molecule%hasCOM
					ENDIF !hoppingsite%hasCOM
				ELSE !k_max=1
					COM_A=system%molecules(molAid)%hoppingsite%subhoppingsite(k)%koord(1:3)			
				END IF !k_max=1
				IF(DEBUG) WRITE(*,*) 'COM_A: ',COM_A						

				!!! get COM for B 	
				DO l=1,l_max
				!!! get COM for B 
						IF(l_max==1) THEN
								IF(system%molecules(molBid)%hoppingsite%hasCOM) THEN !hoppingsite%hasCOM
										COM_B=system%molecules(molBid)%hoppingsite%COM(1:3)
								ELSE !hoppingsite%hasCOM
										IF(system%molecules(molBid)%hasCOM) THEN
											COM_B=system%molecules(molBid)%COM(1:3)
										ELSE
											!! calc COM if not available for all coordinates
											call Massenschwerpunkt(atom(R_u_I(molBid,2):R_u_I(molBid+1,2)-1),&
													&koordA(R_u_I(molBid,2):R_u_I(molBid+1,2)-1,1:3),R_u_I(molBid+1,2)-R_u_I(molBid,2),COM,.false.)
												
											call set_COM_in_Molecule(system%molecules(molBid),COM(1),COM(2),COM(3))
											call set_COM_in_Hoppingsite(system%molecules(molBid)%hoppingsite,COM(1),COM(2),COM(3)) 
											COM_B=COM   								
										ENDIF !molecule%hasCOM
								ENDIF !hoppingsite%hasCOM
						ELSE !l_max=1
									COM_B=system%molecules(molBid)%hoppingsite%subhoppingsite(l)%koord(1:3)			
						END IF !l_max=1
						
						!!! Finally calc the distance vector
						dist_vec(1:3)=COM_B(1:3)-COM_A(1:3)
						!!! correct periodic boundary conditions
						CALL distance_PBC_correction(dist_vec,x_box,y_box,z_box)
						norm_dist_vec=norm(dist_vec)
						
						!!! check if it is the smallest
						IF( norm_dist_vec  < norm_dist_vec_min ) THEN
							dist_vec_min(1:3)=dist_vec(1:3)
							norm_dist_vec_min=norm_dist_vec
						END IF
						IF(DEBUG) WRITE(*,*) 'COM_B: ',COM_B,'|d|:',norm_dist_vec,' dist_vec_min:',dist_vec_min(1:3)
				END DO ! l = N_subhoppingsites in B
			END DO !k = N_subhoppingsites in A

		min_d_COM_arr(i,j,1:3)=dist_vec_min(1:3)
		min_d_COM_arr(i,j,4)=ResidA
		min_d_COM_arr(i,j,5)=ResidB
		min_d_COM_arr(i,j,6)=norm_dist_vec_min
		!WRITE(*,*) min_d_COM_arr(i,j,:)
		IF( norm_dist_vec_min==50000) THEN
			WRITE(*,*) 'ResidA:',ResidA,' ResidB:',ResidB,' |d|=',norm_dist_vec_min
			WRITE(*,*) 'd: ',dist_vec_min(1:3)
			WRITE(*,*) 'Error: no norm_dist_vec encountered, revise input or implementation!'
			CALL EXIT(1)
		END IF
	END DO !j Nearest neighbours
END DO !i =1,N_Resids


WRITE(*,*) 'tmp ENDE N_PBDT_TS1:',i_resid ,system%molecules(i_resid)%hoppingsite%COM(1)
CALL EXIT (2)  !MY EXIT

END SUBROUTINE make_min_d_COM_vector_list


SUBROUTINE Read_neighbourlist(filename,N_Residues,N_Neighbours,COM_Daten,Res_Daten)
        IMPLICIT NONE
        CHARACTER(500), INTENT(IN) :: filename
        INTEGER, INTENT(IN)::N_Residues,N_Neighbours
        REAL, DIMENSION(N_Residues,N_Neighbours+1), INTENT(OUT)::COM_Daten
        INTEGER, DIMENSION(N_Residues,N_Neighbours+1), INTENT(OUT)::Res_Daten
        INTEGER::i,j,NR,dummy,ierror
		LOGICAL::Datei_vorhanden
		
		
        !residue_offset=33 !DIPBI Residue offset um wieder auf die Residuennuemmer zu kommen
        INQUIRE(file=TRIM(filename),exist=Datei_vorhanden) ! Abfrage ob Datei vorhanden ist.
        IF( .NOT. Datei_vorhanden) THEN
            WRITE(*,*) ' Fehler: Die_Datei_ist_nicht_vorhanden. '//TRIM(filename)
            WRITE(*,*) ' ENDE'
            CALL EXIT (1)
        END IF ! Datei-file vorhanden  

        OPEN(UNIT=11,FILE=TRIM(filename),STATUS='OLD',IOSTAT=ierror)

        !Daten mit dem Residuumsindex Nr_Abstand 1 Abstand1 Nr_Abstand 2 Abstand 2 ...bis ... Nr_Abstand 6 Abstand 6
            ! Einlesen der Daten aus der Datei neighbours.ngh
            NR=1
            DO i=1,N_Residues

READ(11,*,IOSTAT=ierror) COM_Daten(NR,1),dummy,COM_Daten(NR,2),COM_Daten(NR,3),COM_Daten(NR,4),COM_Daten(NR,5),COM_Daten(NR,6),&
            & COM_Daten(NR,7),COM_Daten(NR,8),COM_Daten(NR,9),COM_Daten(NR,10),COM_Daten(NR,11),COM_Daten(NR,12),COM_Daten(NR,13)
!WRITE(*,*) COM_Daten(NR,1),dummy,COM_Daten(NR,2),COM_Daten(NR,3),COM_Daten(NR,4),COM_Daten(NR,5),COM_Daten(NR,6)            
!WRITE(*,*) COM_Daten(NR,7),COM_Daten(NR,8),COM_Daten(NR,9),COM_Daten(NR,10),COM_Daten(NR,11),COM_Daten(NR,12),COM_Daten(NR,13)
            
READ(11,*,IOSTAT=ierror) Res_Daten(NR,1),dummy,Res_Daten(NR,2),Res_Daten(NR,3),Res_Daten(NR,4),Res_Daten(NR,5),Res_Daten(NR,6),&
            & Res_Daten(NR,7),Res_Daten(NR,8),Res_Daten(NR,9),Res_Daten(NR,10),Res_Daten(NR,11),Res_Daten(NR,12),Res_Daten(NR,13)
!WRITE(*,*) 'Resids:',Res_Daten(NR,1),dummy,Res_Daten(NR,2),Res_Daten(NR,3),Res_Daten(NR,4),Res_Daten(NR,5),Res_Daten(NR,6),&
!            & Res_Daten(NR,7),Res_Daten(NR,8),Res_Daten(NR,9),Res_Daten(NR,10),Res_Daten(NR,11),Res_Daten(NR,12),Res_Daten(NR,13)	
                NR=NR+1
                IF(ierror < 0) EXIT
                IF(ierror > 0) THEN
			NR=NR-1
WRITE(*,*) 'COM:',COM_Daten(NR,1),dummy,COM_Daten(NR,2),COM_Daten(NR,3),COM_Daten(NR,4),COM_Daten(NR,5),COM_Daten(NR,6)            
WRITE(*,*) COM_Daten(NR,7),COM_Daten(NR,8),COM_Daten(NR,9),COM_Daten(NR,10),COM_Daten(NR,11),COM_Daten(NR,12),COM_Daten(NR,13)

WRITE(*,*) 'Resids:',Res_Daten(NR,1),dummy,Res_Daten(NR,2),Res_Daten(NR,3),Res_Daten(NR,4),Res_Daten(NR,5),Res_Daten(NR,6),&
            & Res_Daten(NR,7),Res_Daten(NR,8),Res_Daten(NR,9),Res_Daten(NR,10),Res_Daten(NR,11),Res_Daten(NR,12),Res_Daten(NR,13)			
                    WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(filename),' aufgetreten.'
		    WRITE(*,*) 'Anzahl der nachbarn nicht 12 ? '
                    WRITE(*,*) 'Beende Einlesen !'
                    CALL EXIT (1)
                END IF
            END DO
            CLOSE(11)
            !NR=1
            !WRITE(*,*) COM_Daten(NR,1),COM_Daten(NR,2),COM_Daten(NR,3),COM_Daten(NR,4),COM_Daten(NR,5),COM_Daten(NR,6)&
            !&,COM_Daten(NR,7),COM_Daten(NR,8),COM_Daten(NR,9),COM_Daten(NR,10),COM_Daten(NR,11),COM_Daten(NR,12),COM_Daten(NR,13)
            !WRITE(*,*) Res_Daten(NR,1),Res_Daten(NR,2),Res_Daten(NR,3),Res_Daten(NR,4),Res_Daten(NR,5),Res_Daten(NR,6)&
            !&,Res_Daten(NR,7),Res_Daten(NR,8),Res_Daten(NR,9),Res_Daten(NR,10),Res_Daten(NR,11),Res_Daten(NR,12),Res_Daten(NR,13)
            !
            !NR=N_Residues
            !WRITE(*,*) COM_Daten(NR,1),COM_Daten(NR,2),COM_Daten(NR,3),COM_Daten(NR,4),COM_Daten(NR,5),COM_Daten(NR,6)&
            !&,COM_Daten(NR,7),COM_Daten(NR,8),COM_Daten(NR,9),COM_Daten(NR,10),COM_Daten(NR,11),COM_Daten(NR,12),COM_Daten(NR,13)
            !WRITE(*,*) Res_Daten(NR,1),Res_Daten(NR,2),Res_Daten(NR,3),Res_Daten(NR,4),Res_Daten(NR,5),Res_Daten(NR,6)&
            !&,Res_Daten(NR,7),Res_Daten(NR,8),Res_Daten(NR,9),Res_Daten(NR,10),Res_Daten(NR,11),Res_Daten(NR,12),Res_Daten(NR,13)

END SUBROUTINE Read_neighbourlist


END PROGRAM main
