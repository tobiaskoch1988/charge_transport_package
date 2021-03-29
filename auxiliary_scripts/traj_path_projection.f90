!traj_path_projection 
PROGRAM traj_path_projection
IMPLICIT NONE
call main()

CONTAINS

SUBROUTINE main
! programm zur projektion der kMC-Trajetorien auf die Boxwaende
IMPLICIT NONE
Real, ALLOCATABLE, Dimension(:,:) :: koordA,COM_Resids
INTEGER*8, ALLOCATABLE, Dimension(:) :: Count_visits,Count_visits_lo,Count_visits_el,Current_Resids,last_Resids ! Zaehlt die Spruenge
INTEGER*8, ALLOCATABLE, Dimension(:) :: Count_pos_xyz ! Zaehlt das Auftreten der Resids 
Character(10), ALLOCATABLE, Dimension(:)::Charge_name_arr
INTEGER*8, ALLOCATABLE, Dimension(:,:,:) :: grid_xyz,grid_xyz_lo,grid_xyz_el !zaehlt die Spruenge
INTEGER*8, ALLOCATABLE, Dimension(:,:,:) :: grid_pos_xyz,grid_pos_xyz_lo,grid_pos_xyz_el  ! Zaehlt die Anzahl der Aufenthalt
INTEGER*8, ALLOCATABLE, Dimension(:)::N_Gesamtlines
Character(3), ALLOCATABLE, Dimension(:):: atomsorteA
Character(500), ALLOCATABLE, Dimension(:)::kmc_traj_arr
Real, Dimension(3) :: koordatom
Character(50) :: Molname
Character(10000)::ursprung,trajectory,no_box_filename,zielfile,lio
Character(500) ::foldername,extrainfo
INTEGER  :: resid
LOGICAL ::DIPBI,DIPBI_KETTE,P3HT,P3MT
LOGICAL::file_exists
Character(50) ::dummyname
Character(5) :: residue_name
Character(3) :: atom
Character(2) ::Ebene
INTEGER::i,i_max,i_min,j,t,t_max,ierror,residue_number,atom_number,MAXATM,Gesamtanzahl,ind
INTEGER::residue_number_min,residue_number_max,last_residue_number,Anzahl_Kantenpunkte
INTEGER::N_Resids,N_Charges,N_bins,Nx_bins,Ny_bins,Nz_bins
INTEGER::ind_x,ind_y,ind_z
INTEGER*8::N_kmc_frames,start_frame,frame
REAL ::x,y,z,v_x,v_y,v_z
REAL ::x_min,y_min,z_min,x_max,y_max,z_max
REAL ::x_box,y_box,z_box,sum_box
REAL ::time_kmc_all,time_kmc,time_kmc_start
LOGICAL::autoscale_on=.false.,make_gnuplot_output=.false.,make_gnu_multiplot=.false.,kmc_traj_list_read=.false.
CHARACTER(10)::f1='(3F14.8)'


x_box=139.0228  ! boxsize in Angstroenm
y_box=269.6753     ! boxsize in Angstroenm
z_box=145.8116    ! boxsize in Angstroenm

trajectory='kmc_traj_V1.xyz'
no_box_filename='no_box_equi_2000ps_500k_def2_theta_75_G0.xyz'
zielfile='Output.dat'
make_gnuplot_output=.false.
i_max=1 ! selects the calculation type
t_max=1 !maximum number of kmc_traj_*.xyz for loop!
! Gibt die Anzahl der Diskretisierungspunkte an an den Kanten (Summe x,y,z)
Anzahl_Kantenpunkte=100


IF(command_argument_count() == 5 .or. command_argument_count() == 8 .or. command_argument_count() == 9 ) THEN         !Einlesen der Flags
        CALL get_command_argument(1,trajectory)
        IF(TRIM(trajectory(len(TRIM(trajectory(:)))-3:)) == '.dat') THEN ! Selects kmc_traj_list_read for *.dat
            kmc_traj_list_read=.true.
        ELSE
            kmc_traj_list_read=.false.
        END IF

        CALL get_command_argument(2,no_box_filename)
        CALL get_command_argument(3,zielfile)
        CALL get_command_argument(4,dummyname)
        read (dummyname, '(I1)') i_max
        IF(i_max > 5 .OR. (i < 0)) THEN
                WRITE(*,*) ' Select: i=(1-5) : 1) jumps_all 2) jumps_el 3) jumps_lo 4) position_all 5) create ALL files '
                STOP
        ELSE
                WRITE(*,*) ' Mode selected: i=',i_max
                WRITE(*,*) ' Select: i=(1-5) : 1) jumps_all 2) jumps_el 3) jumps_lo 4) position_all 5) create ALL files '
        END IF

        CALL get_command_argument(5,dummyname)
        IF(TRIM(adjustl(dummyname))=='gnuplot') THEN
                make_gnuplot_output=.true.
                make_gnu_multiplot=.false.       
        ELSE IF(TRIM(adjustl(dummyname))=='multiplot') THEN
                make_gnuplot_output=.true.
                make_gnu_multiplot=.true.
        ELSE
                make_gnuplot_output=.false.
                make_gnu_multiplot=.false.
        END IF
    IF(command_argument_count() > 7 ) THEN
        WRITE(*,*) 'Reading boxsize in Angstroem'
        CALL get_command_argument(6,dummyname)
        READ(dummyname,'(F12.8)') x_box
        CALL get_command_argument(7,dummyname)
        READ(dummyname,'(F12.8)') y_box
        CALL get_command_argument(8,dummyname)
        READ(dummyname,'(F12.8)') z_box
    END IF
    
    IF(command_argument_count() == 9) THEN
        CALL get_command_argument(9,dummyname)
        READ(dummyname,'(I4)') Anzahl_Kantenpunkte
    END IF
	!CALL get_command_argument(3,calculation_method)
	!CALL get_command_argument(4,b)
ELSE
        WRITE(*,*) ' Script to analyse the kmc_traj*.xyz file using a grid and calculate the projections of the positions'//&
        &' of the charges on the xy,yz,xz plane.'
		WRITE(*,*) ' Select: 1) single kmc_traj*.xyz file OR 2) list of files kmc_traj_list.dat'
        WRITE(*,*) ' no_box*.xyz file is needed with the Hopping-Sites(COM).'
        WRITE(*,*) ' ( If *.dat file is not available a list will be created for all kmc_traj*.xyz in the current folder ) '
        Write(*,*) ' input: trajectory.xyz no_box_filename.xyz zielfilename i_modus [gnuplot/multiplot/false] '//&
                    &'optional: (x_box ,y_box ,z_box  Ang) Number_of_points [50-400] (resolution) '
		Write(*,*) '1) one file example:   kmc_traj_May_V11.xyz no_box***.xyz zielname 4 gnuplot 139.0228 269.6753 145.8118 100 '
        Write(*,*) '2) list files example:   kmc_traj_list.dat no_box***.xyz zielname 4 gnuplot 139.0228 269.6753 145.8118 100 '
        WRITE(*,*) ' Select: i_modus =(1-5) : 1) jumps_all 2) jumps_el 3) jumps_lo 4) position_all 5) create ALL files '
                
		STOP
		Read(*,*) ursprung,no_box_filename,zielfile
END IF

WRITE(*,*) 'use boxsize Ang: ',x_box,y_box,z_box
WRITE(*,*) 'Anzahl der Kantenpunkte: ',Anzahl_Kantenpunkte
WRITE(*,*) 'make_gnuplot_output: ',make_gnuplot_output
WRITE(*,*) 'make_gnu_multiplot: ',make_gnu_multiplot


lio=TRIM(zielfile)//'_lio_'//TRIM(make_lio())
!WRITE(*,*) 'lio: ',TRIM(lio)


IF(TRIM(trajectory(len(TRIM(trajectory(:)))-3:)) == '.xyz') THEN ! Selects kmc_traj_list_read for *.dat
    inquire(file=TRIM(trajectory),exist=file_exists)
    IF(file_exists) THEN
        WRITE(*,*) 'use trajectory file: '//TRIM(trajectory)
        t_max=1
    ELSE
            WRITE(*,*) 'Error: The_file_does_not_exist: '//TRIM(trajectory)
            WRITE(*,*) 'Ende'
            STOP
    END IF
ELSE IF(TRIM(trajectory(len(TRIM(trajectory(:)))-3:)) == '.dat') THEN ! Selects kmc_traj_list_read for *.dat
    CALL make_traj_file_list(trajectory,kmc_traj_arr,N_Gesamtlines,t_max)
ELSE 
    WRITE(*,*) 'Error: The_file_does_not_exist_or_is_not_.xyz_or_.dat: '//TRIM(trajectory)
    WRITE(*,*) 'Ende'
    STOP    
END IF ! trajectory file_exists ?

inquire(file=TRIM(no_box_filename),exist=file_exists)
IF(file_exists) THEN
    WRITE(*,*) 'use no_box.xyz file: '//TRIM(no_box_filename)
ELSE
    WRITE(*,*) 'Error: The_file_does_not_exist: '//TRIM(no_box_filename)
    WRITE(*,*) 'Ende'
    STOP
END IF


! Define the box: Diskretisierungspunkte
sum_box=(x_box+y_box+z_box)

Nx_bins=INT(x_box/sum_box*Anzahl_Kantenpunkte)
Ny_bins=INT(y_box/sum_box*Anzahl_Kantenpunkte)
Nz_bins=INT(z_box/sum_box*Anzahl_Kantenpunkte)

IF(Nx_bins+Ny_bins+Nz_bins > 1000) THEN
    WRITE(*,*) 'Error: Number_of_points_is_to_big_and_may_cause_problems',Nx_bins*Ny_bins*Nz_bins
    WRITE(*,*) 'Choose smaler Anzahl_Kantenpunkte',Anzahl_Kantenpunkte
    STOP
END IF




!!!!!!!!! START READIND DATA !!!!!!!!!!!!!!!!!!!!!
!!! Abschnitt zum Einlesen der no_box datei
WRITE(*,*) 'Start_Reading: '//TRIM(no_box_filename)
i=0
OPEN(UNIT=70,FILE=TRIM(no_box_filename),STATUS='OLD',IOSTAT=ierror)
READ(70,*) N_Resids
WRITE(*,*) 'N_Resids=',N_Resids
!Read(70,*) dummyname                             !Einlesen des Infokopfes im File
!WRITE(*,*) 'Info:',dummyname
! Array mit no_box.xyz Daten der Massenschwerpunkte der Hopping-Sites
ALLOCATE(COM_Resids(N_Resids,3))
ALLOCATE(Count_visits(N_Resids))
ALLOCATE(Count_visits_lo(N_Resids))
ALLOCATE(Count_visits_el(N_Resids))
ALLOCATE(Count_pos_xyz(N_Resids))
Count_visits(:)=0
Count_visits_lo(:)=0
Count_visits_el(:)=0
Count_pos_xyz(:)=0

Do j=1,N_Resids
		Read(70,*,IOSTAT=ierror) atom,COM_Resids(j,1),COM_Resids(j,2),COM_Resids(j,3)
        !WRITE(*,*)  COM_Resids(j,1),COM_Resids(j,2),COM_Resids(j,3)
		IF(ierror < 0) EXIT
		IF(ierror > 0) THEN
			WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(ursprung),' aufgetreten.'
			WRITE(*,*) 'Beende Einlesen !'
			STOP
		END IF
	!WRITE(*,*) residue_number,last_residue_number,atom
	!IF ( 1 == 1 ) THEN
			!WRITE(*,*) adjustl(TRIM(foldername))
			!WRITE(*,'(4ES16.8,A3,100A)') S_AB_HO_to_HO,J_AB_HO_to_HO,S_AB_LU_to_LU,J_AB_LU_to_LU,'   ',adjustl(TRIM(foldername))!,extrainfo
	!END IF 
        IF (i==0 .and. j==N_Resids+1) THEN
            write(*,*) ' Error: Es ist ein Fehler beim Einlesen ausgetreten'
            write(*,*) ' Vermutlich wurde nicht die entsprechende Resid-Nr gefunden',resid
            write(*,*) ' ENDE '
            STOP
        END IF  
END Do
CLOSE(70)
WRITE(*,*) 'End_Reading: '//TRIM(no_box_filename)

!!! End Read no_box.xyz
!j=1
!WRITE(*,*)  COM_Resids(j,1),COM_Resids(j,2),COM_Resids(j,3)
!j=N_Resids
!WRITE(*,*)  COM_Resids(j,1),COM_Resids(j,2),COM_Resids(j,3)




!!! Read trajectory !!!
N_kmc_frames=0
time_kmc=0
DO t=1,t_max ! Read more than one trajectory file ?
    IF ( t_max > 1 ) THEN
        trajectory=kmc_traj_arr(t)
        WRITE(*,*) 'Read file number ',t,' from ',t_max
        Gesamtanzahl=N_Gesamtlines(t)
    ELSE
        Gesamtanzahl=Read_rows_file(trajectory,lio)
    END IF
WRITE(*,*) 'Start_Reading: '//TRIM(trajectory)

ind=1 !!! Hier ist der Index, mit dem die relative Verschiebung Indizes für die Resids angepasste werden kann, falls nicht mit Resid=0 gestartet wird. Default ind=1! für start mit Resid 0 in gro-file
i=0
WRITE(*,*) 'N_lines=',Gesamtanzahl
OPEN(UNIT=77,FILE=TRIM(trajectory),STATUS='OLD',IOSTAT=ierror)
Read(77,*) N_Charges  
WRITE(*,*) 'N_Charges=',N_Charges
Read(77,*) start_frame,dummyname
! Initialisierung
IF(t==1) THEN
    ALLOCATE(Current_Resids(N_Charges))
    ALLOCATE(Last_Resids(N_Charges))
    ALLOCATE(Charge_name_arr(N_Charges))
END IF

DO j=1,N_Charges
    READ(77,*) Charge_name_arr(j),x,y,z,dummyname,Last_Resids(j)
    IF( adjustl(TRIM(Charge_name_arr(j))) == 'electron' ) THEN
                    Count_visits_el(Last_Resids(j)+ind)=Count_visits_el(Last_Resids(j)+ind)+1
    ELSE IF( adjustl(TRIM(Charge_name_arr(j))) == 'hole' ) THEN
                    Count_visits_lo(Last_Resids(j)+ind)=Count_visits_lo(Last_Resids(j)+ind)+1
    END IF ! electron or hole ?
    Count_visits(Last_Resids(j)+ind)=Count_visits(Last_Resids(j)+ind)+1 ! Startpositionen zaehlen
    Count_pos_xyz(Last_Resids(j)+ind)=Count_pos_xyz(Last_Resids(j)+ind)+1 ! Startpositionen zaehlen
    !WRITE(*,*) Charge_name_arr(j),x,y,z,TRIM(dummyname),Last_Resids(j)
END DO 
WRITE(*,*) 'Ende Reading: '//TRIM(trajectory)

!!Check if the Data in the no_box.xyz and traj.xyz and the index_shift match.


koordatom(1)=x-COM_Resids(Last_Resids(N_charges)+ind,1)
koordatom(2)=y-COM_Resids(Last_Resids(N_charges)+ind,2)
koordatom(3)=z-COM_Resids(Last_Resids(N_charges)+ind,3)
!WRITE(*,*) koordatom(:),COM_Resids(Last_Resids(N_charges)+ind,:)
IF( norm(koordatom) < 1.0E-5 ) THEN
    WRITE(*,*) 'Check ok!' 
ELSE
    WRITE(*,*) 'Error: The coorindates in the files'//TRIM(trajectory)//' and '//TRIM(no_box_filename)//&
    &' do not match, or there is an index error *.f90 file Resids start with index 0 ?'
    WRITE(*,*) 'ENDE'
    STOP
END IF



!!! START CALCULATIONS !!!

!write(*,*) 'Anzahl der Atome in .gro Datei ', Gesamtanzahl
Do i=1,INT(Gesamtanzahl/(N_Charges+2))
        Read(77,*,IOSTAT=ierror) N_Charges 
        IF(ierror < 0) EXIT 
        IF(i==1) THEN
            Read(77,*,IOSTAT=ierror) start_frame,dummyname,dummyname,time_kmc_start
        ELSE
            Read(77,*,IOSTAT=ierror) frame,dummyname,dummyname,time_kmc
        END IF
        IF(ierror < 0) EXIT
        DO j=1,N_Charges
            READ(77,*,IOSTAT=ierror) Charge_name_arr(j),x,y,z,dummyname,Current_Resids(j)
            IF(ierror < 0) EXIT
            IF(ierror > 0) THEN
                WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(trajectory),' aufgetreten.'
                WRITE(*,*) 'Beende Einlesen !'
                STOP
            END IF
            !WRITE(*,*) Charge_name_arr(j),x,y,z,TRIM(dummyname),Current_Resids(j)
            IF( Current_Resids(j) /= Last_Resids(j) ) THEN ! find jump
                IF( adjustl(TRIM(Charge_name_arr(j))) == 'electron' ) THEN
                    Count_visits_el(Current_Resids(j)+ind)=Count_visits_el(Current_Resids(j)+ind)+1
                ELSE IF( adjustl(TRIM(Charge_name_arr(j))) == 'hole' ) THEN
                    Count_visits_lo(Current_Resids(j)+ind)=Count_visits_lo(Current_Resids(j)+ind)+1
                END IF ! electron or hole ?
                Count_visits(Current_Resids(j)+ind)=Count_visits(Current_Resids(j)+ind)+1
            END IF ! jump found ?
            ! count all positions / even if charge did not jump
            Count_pos_xyz(Current_Resids(j)+ind)=Count_pos_xyz(Current_Resids(j)+ind)+1 ! Zaehle immer!
        END DO 
        ! Check changes
        
        Last_Resids(:)=Current_Resids(:)

		IF(ierror < 0) EXIT
		IF(ierror > 0) THEN
			WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(trajectory),' aufgetreten.'
			WRITE(*,*) 'Beende Einlesen !'
			STOP
		END IF        
END Do
CLOSE(77)

N_kmc_frames=N_kmc_frames+frame-start_frame
time_kmc_all=time_kmc_all+time_kmc-time_kmc_start
END DO ! Read input data from trajectory t=1, to t_max ; if to t_max /=1 then loop for more trajectory files
!!! END READ DATA !!!!


IF ( sum(Count_visits(:))-N_Charges==0 ) THEN
    WRITE(*,*) 'Error: No_jump_encountered_in_'//TRIM(trajectory)//'__check_input!'
    WRITE(*,*) 'Ende'
    STOP
ELSE
    WRITE(*,*) 'N_Jumps encountered:',sum(Count_visits(:))-N_Charges
END IF



Allocate(grid_xyz(Nx_bins,Ny_bins,Nz_bins))
Allocate(grid_xyz_lo(Nx_bins,Ny_bins,Nz_bins))
Allocate(grid_xyz_el(Nx_bins,Ny_bins,Nz_bins))
Allocate(grid_pos_xyz(Nx_bins,Ny_bins,Nz_bins))

WRITE(*,*) 'Nx_bin:',Nx_bins,'Ny_bin:',Ny_bins,'Nz_bin:',Nz_bins
WRITE(*,*) 'x_box:',x_box,'y_box',y_box,'z_box',z_box

grid_xyz=0
grid_xyz_lo=0
grid_xyz_el=0
grid_pos_xyz=0
WRITE(*,*) 'N_Resids:',N_Resids
DO i=1,N_Resids
        ind_x=int(Nx_bins*COM_Resids(i,1)/x_box)+1
        ind_y=int(Ny_bins*COM_Resids(i,2)/y_box)+1
        ind_z=int(Nz_bins*COM_Resids(i,3)/z_box)+1
        grid_xyz(ind_x,ind_y,ind_z)=Count_visits(i) ! jumps
        grid_xyz_lo(ind_x,ind_y,ind_z)=Count_visits_lo(i)
        grid_xyz_el(ind_x,ind_y,ind_z)=Count_visits_el(i)
        grid_pos_xyz(ind_x,ind_y,ind_z)=Count_pos_xyz(i) ! All
        !WRITE(*,*) COM_Resids(i,:),ind_x,ind_y,ind_z,Count_visits(i)
END DO


! Vorbereitung fuer mehr als einen gnuplot output
IF(i_max==5) THEN ! alle
    i_min=1
    i_max=4
ELSE
    i_min=i_max
END IF

!!!!! Section to print the gnuplot files and make *.pdf !!!!!!!
WRITE(*,*) 'Selected mode: i_min=',i_min,'i_max: ',i_max
dummyname=TRIM(zielfile)
! Schleife zur Auswertung und zum plotten der Daten i=1 jumps (default =gesamt) , i=2 (elektronen), i=3 (loecher) , i=4 (Alle
! Zaehlungen)
DO i=i_min,i_max
        IF( i==1) THEN
                zielfile='jumps_all_'//TRIM(dummyname)
        ELSE IF(i==2) THEN
                zielfile='jumps_lo_'//TRIM(dummyname)
                grid_xyz=grid_xyz_lo
        ELSE IF (i==3) THEN
                zielfile='jumps_el_'//TRIM(dummyname)
                grid_xyz=grid_xyz_el
        ELSE IF (i==4) THEN
                zielfile='pos_all_'//TRIM(dummyname)
                grid_xyz= grid_pos_xyz
        END IF
WRITE(*,*) 'use: ',TRIM(zielfile)                
Ebene='xy'
OPEN(UNIT=80,FILE=TRIM(Ebene)//'_proj_'//TRIM(zielfile)//'.dat',STATUS='REPLACE',IOSTAT=ierror)
WRITE(*,*) '### '//TRIM(Ebene)//'-projection for '//TRIM(trajectory)
WRITE(80,*) '### '//TRIM(Ebene)//'-projection for '//TRIM(trajectory)
! Print output xy-Ebene Projektion
z_max=0
Do ind_x=1,Nx_bins
    Do ind_y=1,Ny_bins
        z=sum(grid_xyz(ind_x,ind_y,:))
        If( z > z_max) THEN
            z_max=z
        END IF
        !WRITE(80,*) ind_x/REAL(Nx_bins)*x_box,ind_y/REAL(Ny_bins)*y_box,int(z)
    END DO ! y
END DO

Do ind_x=1,Nx_bins
    Do ind_y=1,Ny_bins
        z=sum(grid_xyz(ind_x,ind_y,:))
        WRITE(80,f1) (ind_x-1)/REAL(Nx_bins-1)*x_box,(ind_y-1)/REAL(Ny_bins-1)*y_box,int(z)/Real(z_max)
    END DO ! y
END DO
CLOSE(80)

x_min=0.0  !x_box/REAL(Nx_bins)
y_min=0.0  !y_box/REAL(Ny_bins)
z_min=0.0
z_max=1.0
call make_gnuplot_file(zielfile,Ebene,x_min,x_box,y_min,y_box,z_min,z_max,autoscale_on,make_gnuplot_output,Nx_bins,Ny_bins) 
!! Ende plot xy-Ebene


Ebene='xz'
OPEN(UNIT=80,FILE=TRIM(Ebene)//'_proj_'//TRIM(zielfile)//'.dat',STATUS='REPLACE',IOSTAT=ierror)
WRITE(*,*) '### '//TRIM(Ebene)//'-projection for '//TRIM(trajectory)
WRITE(80,*) '### '//TRIM(Ebene)//'-projection for '//TRIM(trajectory)
! Print output xz-Ebene Projektion
y_max=0
Do ind_x=1,Nx_bins
    Do ind_z=1,Nz_bins
            y=sum(grid_xyz(ind_x,:,ind_z))
            IF( y > y_max) THEN
                y_max=y
            END IF
            !WRITE(80,*) ind_x/REAL(Nx_bins)*x_box,ind_z/REAL(Nz_bins)*z_box,int(y)
     END DO ! z
END DO !x

DO ind_x=1,Nx_bins
        DO ind_z=1,Nz_bins
           y=sum(grid_xyz(ind_x,:,ind_z))
           WRITE(80,f1) (ind_x-1)/REAL(Nx_bins-1)*x_box,(ind_z-1)/REAL(Nz_bins-1)*z_box,int(y)/Real(y_max)
        END DO ! z
END DO !x
CLOSE(80)

x_min=0.0  !x_box/REAL(Nx_bins)
z_min=0.0  !z_box/REAL(Nz_bins)
y_min=0.0
y_max=1.0
call make_gnuplot_file(zielfile,Ebene,x_min,x_box,z_min,z_box,y_min,y_max,autoscale_on,make_gnuplot_output,Nx_bins,Nz_bins)
!!! Ende  plot xz-Ebene



!!! yz-Ebene
Ebene='yz'
OPEN(UNIT=80,FILE=TRIM(Ebene)//'_proj_'//TRIM(zielfile)//'.dat',STATUS='REPLACE',IOSTAT=ierror)
WRITE(*,*) '### '//TRIM(Ebene)//'-projection for '//TRIM(trajectory)
WRITE(80,*) '### '//TRIM(Ebene)//'-projection for '//TRIM(trajectory)
! Print output xz-Ebene Projektion
x_max=0
Do ind_y=1,Ny_bins
    Do ind_z=1,Nz_bins
       x=sum(grid_xyz(:,ind_y,ind_z))
       IF( x > x_max) THEN
           x_max=x
       END IF
       !WRITE(80,*) ind_y/REAL(Ny_bins)*y_box,ind_z/REAL(Nz_bins)*z_box,int(x)
     END DO ! z
END DO ! y

DO ind_y=1,Ny_bins
   DO ind_z=1,Nz_bins
      x=sum(grid_xyz(:,ind_y,ind_z))
      WRITE(80,f1) (ind_y-1)/REAL(Ny_bins-1)*y_box,(ind_z-1)/REAL(Nz_bins-1)*z_box,int(x)/Real(x_max)
   END DO ! z
END DO !y
CLOSE(80)

y_min=0.0  !y_box/REAL(Ny_bins)
z_min=0.0  !z_box/REAL(Nz_bins)
x_min=0.0
x_max=1.0

call make_gnuplot_file(zielfile,Ebene,y_min,y_box,z_min,z_box,x_min,x_max,autoscale_on,make_gnuplot_output,Ny_bins,Nz_bins)
! Ende plot yz-Ebene 

IF(make_gnu_multiplot) THEN
    CALL  make_gnu_multiplot_file(zielfile,x_box,y_box,z_box,autoscale_on,make_gnuplot_output)
END IF

END DO ! i=1,i_max (Auswahl der Daten)



WRITE(*,*) 'Nx_bin:',Nx_bins,'Ny_bin:',Ny_bins,'Nz_bin:',Nz_bins
WRITE(*,*) 'time_kmc_all [s]:',time_kmc_all
WRITE(*,*) 'N_kmc_frames:',N_kmc_frames
WRITE(*,*) 'Normal termination of traj_path_projection'

END SUBROUTINE main

SUBROUTINE make_gnuplot_file(zielfile,Ebene,x_min,x_box,y_min,y_box,z_min,z_max,autoscale_on,make_gnuplot_output,xpix,ypix) 
IMPLICIT NONE
Character(10000),INTENT(IN) ::zielfile
Character(2),INTENT(IN) ::Ebene
Real,INTENT(IN)    ::x_box,y_box,z_max
REAL, INTENT(IN)   ::x_min,y_min,z_min
INTEGER, INTENT(IN) ::xpix,ypix    !Aufloesung Anzahl der Pixel in x,y-Richtung
LOGICAL,INTENT(IN) ::autoscale_on,make_gnuplot_output
Character(10000)   :: gnuplotfilename,datenfile,bashline
INTEGER::ierror

gnuplotfilename='gnuplot_'//TRIM(Ebene)//'_proj_'//TRIM(zielfile)
datenfile=TRIM(Ebene)//'_proj_'//TRIM(zielfile)//'.dat'

!!! gnuplot setings
OPEN(UNIT=81,FILE=TRIM(gnuplotfilename),STATUS='REPLACE',IOSTAT=ierror)

WRITE(81,*) ' set view map                                                                                                      '
WRITE(81,*) ' set dgrid3d '//TRIM(str(xpix))//','//TRIM(str(ypix))//',1  '
WRITE(81,*) ' ##set pm3d interpolate 0,0                                                                              '
WRITE(81,*) ' ##splot "'//TRIM(datenfile)//'" u 1:2:3 title " title " with pm3d                                '
WRITE(81,*) ' set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0) '
WRITE(81,*) ' set xlabel "'//TRIM(Ebene(1:1))//'  / \305 " font "Helvetica,16"       '
WRITE(81,*) ' set ylabel "'//TRIM(Ebene(2:2))//'  / \305 " font "Helvetica,16"       '
WRITE(81,*) ' set zlabel " density " font "Helvetica,16"                                                                        '
WRITE(81,*) ' set xtics font "Helvetica,16"                                                                                     '
WRITE(81,*) ' set ytics font "Helvetica,16"                                                                                     '
WRITE(81,*) ' set ztics font "Helvetica,16"                                                                                     '
WRITE(81,*) ' set xrange [ ',real(x_min),':',real(x_box),']                                                                       ' ! set xrange [',int(1.0/Nx_bins*x_box),':',int(x_box),']       
WRITE(81,*) ' set yrange [ ',real(y_min),':',real(y_box),']                                                                       '
WRITE(81,*) ' set zrange [ 0.0001:',real(z_max),']                            '! set zrange [',real(z_min),':',real(z_max),'] 
WRITE(81,*) ' set logscale z                                                                                                 '
WRITE(81,*) ' set logscale cb                                                                                                 '
WRITE(81,*) ' set cbrange[0.0001:',real(z_max),']                                                                                 '
WRITE(81,*) ' set cblabel "density" font "Helvetica,16"                                                                  '
iF(autoscale_on) THEN
    WRITE(81,*) ' set autoscale '
END IF

WRITE(81,*) " set object 1 rectangle from graph 0,0  to graph 1,1 behind fillcolor rgb 'black' fillstyle solid noborder"
WRITE(81,*) ' set terminal pdf                                                                                           '
WRITE(81,*) " set output '"//TRIM(Ebene)//'_proj_'//TRIM(zielfile)//".pdf'                                               "
WRITE(81,*) ' ### 4-Ecken-Interpolation '
WRITE(81,*) ' ##set pm3d interpolate 0,0                                                                                   '
WRITE(81,*) ' ##splot "'//TRIM(datenfile)//'" u 1:2:3 title "COM" with pm3d                                                '
WRITE(81,*) ' '
WRITE(81,*) 'set style fill transparent solid 0.9 noborder'
WRITE(81,*) '## plot "'//TRIM(datenfile)//'" u 1:2:3 with circles lc palette notitle'
WRITE(81,*) '#### kleine ausgefuellte Kreise'
WRITE(81,*) 'plot "'//TRIM(datenfile)//'" u 1:2:3 with point pointtype 7 ps 0.3 lc palette notitle'
CLOSE(81)

IF(make_gnuplot_output) THEN
        write(bashline,*) ' gnuplot '//TRIM(gnuplotfilename)
        Write(*,*) trim(bashline)
        CALL execute_command_line(TRIM(bashline))
END IF

END SUBROUTINE make_gnuplot_file




SUBROUTINE make_gnu_multiplot_file(zielfile,x_box,y_box,z_box,autoscale_on,make_gnuplot_output)
Character(10000),INTENT(IN) ::zielfile
REAL, INTENT(IN)   ::x_box,y_box,z_box
LOGICAL,INTENT(IN) ::autoscale_on,make_gnuplot_output
Character(10000)   :: gnuplotfilename,bashline,pdf_filename
REAL::mpl_top,mpl_bot,mpl_left,mpl_right,mpl_height,mpl_width,mpl_dx,mpl_dy
INTEGER::mpl_ny,mpl_nx  
REAL::xsize,ysize
REAL,ALLOCATABLE,DIMENSION (:)::bot,top,left,right
INTEGER::n,ierror
mpl_top    = 0.4 !! #inch  outer top margin, title goes here
mpl_bot    = 0.7 !! #inch  outer bottom margin, x label goes here
mpl_left   = 0.9 !! #inch  outer left margin, y label goes here
mpl_right  = 0.1 !! #inch  outer right margin, y2 label goes here
mpl_height = 3.5 !! #inch  height of individual plots
mpl_width  = 3.0 !! #inch  width of individual plots
mpl_dx     = 0.2 !! #inch  inter-plot horizontal spacing
mpl_dy     = 0.3 !! #inch  inter-plot vertical spacing
mpl_ny     = 2   !! #number of rows
mpl_nx     = 2   !! #number of columns

gnuplotfilename='gnuplot_multi_'//TRIM(zielfile)
pdf_filename='multi_'//TRIM(zielfile)
!# calculate full dimensions
xsize = mpl_left+mpl_right+(mpl_width*mpl_nx)+(mpl_nx-1)*mpl_dx
ysize = mpl_top+mpl_bot+(mpl_ny*mpl_height)+(mpl_ny-1)*mpl_dy

ALLOCATE(left(mpl_nx),right(mpl_nx))
ALLOCATE(top(mpl_ny),bot(mpl_ny))

!# placement functions
!#   rows are numbered from bottom to top
DO n=1,mpl_ny
    bot(n) = (mpl_bot+(n-1)*mpl_height+(n-1)*mpl_dy)/ysize
    top(n)  = 1-((mpl_top+(mpl_ny-n)*(mpl_height+mpl_dy))/ysize)
END DO
!#   columns are numbered from left to right
DO n=1,mpl_nx 
    left(n) = (mpl_left+(n-1)*mpl_width+(n-1)*mpl_dx)/xsize
    right(n)  = 1-((mpl_right+(mpl_nx-n)*(mpl_width+mpl_dx))/xsize)
END DO

OPEN(UNIT=85,FILE=TRIM(gnuplotfilename),STATUS='REPLACE',IOSTAT=ierror)
WRITE(85,*) 'set terminal pdf enhanced color dl 2.0 size ',xsize,',',ysize
WRITE(85,*) '##set terminal pdf enhanced color dl 2.0 size ',xsize,',',ysize 
WRITE(85,*) 'set encoding iso_8859_1'
WRITE(85,*) 'set tics scale 1.5'

WRITE(85,*) "set output '"//TRIM(pdf_filename)//".pdf' "

WRITE(85,*) 'set offsets       '
WRITE(85,*) 'set autoscale fix '
WRITE(85,*) 'set size 1,1      '
WRITE(85,*) 'set nokey         '

WRITE(85,*) 'set xtics font "Helvetica,16"  '                                                                                   
WRITE(85,*) 'set ytics font "Helvetica,16"  '                                                                                  
WRITE(85,*) 'set ztics font "Helvetica,16"  '

WRITE(85,*) "set object 1 rectangle from graph 0,0  to graph 1,1 behind fillcolor rgb 'black' fillstyle solid noborder"


! # start plotting
WRITE(85,*) ' set multiplot'

WRITE(85,*) '#-----------------------------------------------'
WRITE(85,*) '# subplot  1-2                                  '
WRITE(85,*) ' #  set horizontal margins for first column     '
WRITE(85,*) ' set lmargin at screen ',left(1)                                                                                    
WRITE(85,*) ' set rmargin at screen ',right(1)                                                                                   
WRITE(85,*) ' #  set horizontal margins for second row (middle)                                                                '
WRITE(85,*) ' set tmargin at screen ',top(2)                                                                                     
WRITE(85,*) ' set bmargin at screen ',bot(2)                                                                                     
WRITE(85,*) '                                                                                                                  '
WRITE(85,*) ' set title ''                                                                                                     '
WRITE(85,*) ' set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)'
WRITE(85,*) ' #set xlabel "x  / \305 " font "Helvetica,16"                                                                     '
WRITE(85,*) ' set ylabel "z  / \305 " font "Helvetica,16"                                                                      '
WRITE(85,*) ' set zlabel " density " font "Helvetica,16"                                                                       '
WRITE(85,*) ' set xtics font "Helvetica,16"                                                                                    '
WRITE(85,*) ' set ytics font "Helvetica,16"                                                                                    '
WRITE(85,*) ' set ztics font "Helvetica,16"                                                                                    '
WRITE(85,*) ' set xrange [    0.00000000     :   ',real(x_box),'  ]                                                         '      
WRITE(85,*) ' set yrange [    0.00000000     :   ',real(z_box),'  ]                                                         '      
WRITE(85,*) ' set zrange [ 0.0001:   1.00000000     ]                                                                          '
WRITE(85,*) ' set logscale z                                                                                                   '
WRITE(85,*) ' set logscale cb                                                                                                  '
WRITE(85,*) ' set cbrange[0.0001:   1.00000000     ]                                                                           '
iF(autoscale_on) THEN
    WRITE(85,*) ' set autoscale '
END IF
WRITE(85,*) ' unset colorbox                                                                                                   '
WRITE(85,*) ' set style fill transparent solid 0.99 noborder                                                                    '
WRITE(85,*) ' ## plot "xz_proj_'//TRIM(zielfile)//'.dat" u 1:2:3 with circles lc palette notitle                                '
WRITE(85,*) ' plot "xz_proj_'//TRIM(zielfile)//'.dat" u 1:2:3 with point pointtype 7 ps 0.3 lc palette notitle           '
WRITE(85,*) ' '  
WRITE(85,*) ' ;'

WRITE(85,*) '#-----------------------------------------------              '
WRITE(85,*) '# subplot  2-2 yz-Ebene                                       '
WRITE(85,*) '#  set horizontal margins for second column                   '
WRITE(85,*) 'set lmargin at screen ',left(2),'                                                                                  '
WRITE(85,*) 'set rmargin at screen ',right(2),'                                                                                 '
WRITE(85,*) '#  set horizontal margins for second row (middle)                                                                  '
WRITE(85,*) 'set tmargin at screen ',top(2),'                                                                                   '
WRITE(85,*) 'set bmargin at screen ',bot(2),'                                                                                   '
WRITE(85,*) '                                                                                                                   '
WRITE(85,*) '  set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)'
WRITE(85,*) '  set xlabel "y  / \305 " font "Helvetica,16"                                                                      '
WRITE(85,*) '  unset ylabel ##"z  / \305 " font "Helvetica,16"                                                                  '
WRITE(85,*) '  set format y " "                                                                                                 '
WRITE(85,*) '  set zlabel " density " font "Helvetica,16"                                                                       '
WRITE(85,*) '  set xtics font "Helvetica,16"                                                                                    '
WRITE(85,*) '  set ytics font "Helvetica,16"                                                                                    '
WRITE(85,*) '  set ztics font "Helvetica,16"                                                                                    '
WRITE(85,*) '  set xrange [    0.00000000     :  ',real(y_box),'   ]  '   
WRITE(85,*) '  set yrange [    0.00000000     :  ',real(z_box),'   ]  '   
WRITE(85,*) '  set zrange [ 0.0001:   1.00000000     ]                                                                          '
WRITE(85,*) '  set logscale z                                                                                                   '
WRITE(85,*) '  set logscale cb                                                                                                  '
WRITE(85,*) '  set cbrange[0.0001:   1.00000000     ]                                                                           '
iF(autoscale_on) THEN
    WRITE(85,*) ' set autoscale '
END IF
WRITE(85,*) '  unset colorbox                                                                                                   '
WRITE(85,*) '  set format z ""                                                                                                  '
WRITE(85,*) '  #set cblabel "density" font "Helvetica,16"                                                                       '
WRITE(85,*) '                                                                                                                   '         
WRITE(85,*) '  ### 4-Ecken-Interpolation                                                                                        '
WRITE(85,*) '  ##set pm3d interpolate 0,0                                                                                       '
WRITE(85,*) '  ##splot "yz_proj_pos_all_zielA.dat" u 1:2:3 title "COM" with pm3d                                                '
WRITE(85,*) '                                                                                                                   '
WRITE(85,*) ' set style fill transparent solid 0.99 noborder                                                                     '
WRITE(85,*) ' ## plot "yz_proj_'//TRIM(zielfile)//'.dat" u 1:2:3 with circles lc palette notitle                                '
WRITE(85,*) ' #### kleine ausgefuellte Kreise                                                                                   '
WRITE(85,*) ' plot "yz_proj_'//TRIM(zielfile)//'.dat" u 1:2:3 with point pointtype 7 ps 0.3 lc palette notitle                 '
WRITE(85,*) ';                                                                                                                  '
WRITE(85,*) '                                                                                                                   '
WRITE(85,*) '#-----------------------------------------------                                                                   '
WRITE(85,*) '# subplot  1-1                                                                                                     '
WRITE(85,*) '#  set horizontal margins for first column                                                                         '
WRITE(85,*) 'set lmargin at screen ',left(1)
WRITE(85,*) 'set rmargin at screen ',right(1) 
WRITE(85,*) '#  set horizontal margins for first row (bottom)                                                                   '
WRITE(85,*) 'set tmargin at screen ',top(1)  
WRITE(85,*) 'set bmargin at screen ',bot(1) 
WRITE(85,*) '     '
WRITE(85,*) 'set title ''                                                                                                        '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) '# now set a label and tic marks for the x-axis                                                                      '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) '  set xlabel "x  / \305 " font "Helvetica,16"                                                                       '
WRITE(85,*) '  set ylabel "y  / \305 " font "Helvetica,16"                                                                       '
WRITE(85,*) '  set format y                                                                                                      '

WRITE(85,*) '  set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0) '
WRITE(85,*) '  set xrange [    0.00000000     :   ',real(x_box),'   ]                                                            '   
WRITE(85,*) '  set yrange [    0.00000000     :   ',real(y_box),'   ]                                                            '   
WRITE(85,*) '  set zrange [ 0.0001:   1.00000000     ]                                                                           '
WRITE(85,*) '  set logscale z                                                                                                    '
WRITE(85,*) '  set logscale cb                                                                                                   '
WRITE(85,*) '  set cbrange[0.0001:   1.00000000     ]                                                                            '
WRITE(85,*) '  set cblabel "density" font "Helvetica,16"                                                                         '
WRITE(85,*) '  set colorbox                                                                                                      '
WRITE(85,*) '                                                                                                                    '
iF(autoscale_on) THEN
    WRITE(85,*) ' set autoscale '
END IF
WRITE(85,*) ' set xtics font "Helvetica,16"                                                                                      '
WRITE(85,*) ' set ytics font "Helvetica,16"                                                                                      '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) 'set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0                                               '
WRITE(85,*) 'set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0                                              '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) ' set style fill transparent solid 0.99 noborder                                                                     '
WRITE(85,*) ' ## plot "xy_proj_'//TRIM(zielfile)//'.dat" u 1:2:3 with circles lc palette notitle                                 '
WRITE(85,*) ' plot "xy_proj_'//TRIM(zielfile)//'.dat" u 1:2:3 with point pointtype 7 ps 0.3 lc palette notitle                   '
WRITE(85,*) ';                                                                                                                   '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) '#-----------------------------------------------                                                                    '
WRITE(85,*) '# subplot  1-2                                                                                                      '
WRITE(85,*) '#  set horizontal margins for second column                                                                         '
WRITE(85,*) '#set lmargin at screen left(2)                                                                                      '
WRITE(85,*) '#set rmargin at screen right(2)                                                                                     '
WRITE(85,*) '#  set horizontal margins for first row (bottom)                                                                    '
WRITE(85,*) '#set tmargin at screen top(1)                                                                                       '
WRITE(85,*) '#set bmargin at screen bot(1)                                                                                       '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) '#set title ''                                                                                                       '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) '#set ylabel ""             # no label here                                                                          '
WRITE(85,*) '#set yrange [0:1.5]                                                                                                 '
WRITE(85,*) '#set format y ""           # no tic labels                                                                          '
WRITE(85,*) '#set ytics mirror 1                                                                                                 '
WRITE(85,*) '#set mytics 2                                                                                                       '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) '#set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0                                              '
WRITE(85,*) '#set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0                                             '
WRITE(85,*) '                                                                                                                    '
WRITE(85,*) '#plot          \                                                                                                    '
WRITE(85,*) '#(cos(x))**2 \                                                                                                      '
WRITE(85,*) '#axes x1y1 \                                                                                                        '
WRITE(85,*) '#title '' \                                                                                                         '
WRITE(85,*) '#with lines lt 1 lc 6 lw 2\                                                                                         '
WRITE(85,*) '#;                                                                                                                  '
WRITE(85,*) 'unset multiplot                                                                                                     '

CLOSE(85)

IF(make_gnuplot_output) THEN
        write(bashline,*) ' gnuplot '//TRIM(gnuplotfilename)
        Write(*,*) trim(bashline)
        CALL execute_command_line(TRIM(bashline))
END IF

END SUBROUTINE make_gnu_multiplot_file


SUBROUTINE make_traj_file_list(trajectory,kmc_traj_arr,N_Gesamtlines,t_max)
 ! subroutine to read the kmc_traj*.xyz filenames in an array kmc_traj_arr in the working directory
 Character(500)::trajectory
 Character(500), ALLOCATABLE, Dimension(:),INTENT(OUT)::kmc_traj_arr
 INTEGER*8, ALLOCATABLE, Dimension(:),INTENT(OUT)::N_Gesamtlines
 INTEGER, INTENT(OUT)::t_max  ! maximum number of readable kmc_traj*.xyz files
 Character(10000) ::trajectory_list_name
 Character(12000) :: bashline
 LOGICAL::file_exists
 INTEGER::i,t,N_read,ierror,N_lines
 Character(10000)::lio
 Character(10000)::filename
 trajectory_list_name=TRIM(trajectory)
 lio=TRIM(make_lio())
 inquire(file=TRIM(trajectory_list_name),exist=file_exists)
 IF(file_exists) THEN
        WRITE(*,*) 'use trajectory list file: '//TRIM(trajectory_list_name)
 ELSE ! make list 
        inquire(file=TRIM(trajectory_list_name),exist=file_exists)
        WRITE(*,*) 'Create new file kmc_traj*.dat with list of kmc_traj*.xyz:'//TRIM(trajectory_list_name)
        write(bashline,*) 'ls kmc_traj*.xyz > '//TRIM(trajectory_list_name)
        Write(*,*) trim(bashline)
        CALL execute_command_line(TRIM(bashline)) 
        inquire(file=TRIM(trajectory_list_name),exist=file_exists)
 END IF ! trajectory_list_name file_exists ?

IF(file_exists) THEN
    N_read=Read_rows_file(trajectory_list_name,lio) 
    ALLOCATE(kmc_traj_arr(N_read)) ! only the t_max usefull once will be stored afterwards.
    ALLOCATE(N_Gesamtlines(N_read)) ! saves the linenumbers 
    OPEN(UNIT=87,FILE=TRIM(trajectory_list_name),STATUS='old',IOSTAT=ierror)
ELSE
    WRITE(*,*) 'Error: The_file_does_not_exist: '//TRIM(trajectory_list_name)
    WRITE(*,*) 'Ende'
    STOP
END IF


WRITE(*,*) 'Start reading List of files kmc_traj*.xyz ...'
!! Read filenames
t_max=0
DO i=1,N_read
    READ(87,*) filename
    IF(ierror < 0) EXIT
    IF(ierror > 0) THEN
        WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(trajectory_list_name),' aufgetreten.'
        WRITE(*,*) 'Beende Einlesen !'
        STOP
    END IF    
    
    inquire(file=TRIM(filename),exist=file_exists)
    IF(file_exists) THEN
        N_lines=Read_rows_file(filename,lio)
        IF( N_lines > 100) THEN
            t_max=t_max+1
            N_Gesamtlines(t_max)=N_lines
            kmc_traj_arr(t_max)=TRIM(filename)
        ELSE
            WRITE(*,*) 'Small amount of data, file skipped: ',TRIM(filename)
        END IF
    ELSE
        WRITE(*,*) 'File skipped: ',TRIM(filename)
    END IF    
END DO

WRITE(*,*) 'Number of kmc_traj*.xyz files: ',t_max
WRITE(*,*) 'N_lines  kmc_traj*.xyz'
DO i=1,t_max
    WRITE(*,*) N_Gesamtlines(i),TRIM(kmc_traj_arr(i))
END DO
WRITE(*,*) 'Ende List of files'
END SUBROUTINE make_traj_file_list




SUBROUTINE read_coord_from_gro_resid(ursprung,resid,atomsorteA,koordA,Molname,DIPBI,DIPBI_KETTE,P3HT,P3MT)
    ! SUBROUTINE zum Einlesen von Molekuelkoordinaten aus einer gro-Datei
    ! bei gegebener der resid
    Real, ALLOCATABLE, Dimension(:,:), INTENT(INOUT) :: koordA
    Character(3), ALLOCATABLE, Dimension(:), INTENT(INOUT) :: atomsorteA
    Character(50), INTENT (OUT) :: Molname
    Character(50), INTENT (IN) :: ursprung
    INTEGER,       INTENT (IN) :: resid
    LOGICAL, INTENT (INOUT) ::DIPBI,DIPBI_KETTE,P3HT,P3MT
    Character(50)             ::dummyname
    Character(5) :: residue_name
    Character(3) :: atom
    INTEGER::i,j,ierror,residue_number,atom_number,MAXATM,Gesamtanzahl
    REAL::x,y,z,v_x,v_y,v_z
    
    OPEN(UNIT=78,FILE=TRIM(ursprung),STATUS='OLD',IOSTAT=ierror)
    Read(78,*) dummyname                             !Einlesen des Infokopfes im File
    Read(78,*) Gesamtanzahl                          !Einlesen der Atomanzahl im File
    !write(*,*) 'Anzahl der Atome in .gro Datei ', Gesamtanzahl
    
    i=0
    DO j=1,Gesamtanzahl
        Read(78,'(i5,2a5,i5,3f8.3,3f8.4)',IOSTAT=ierror) residue_number,residue_name,atom,atom_number,x,y,z,v_x,v_y,v_z
        IF(ierror < 0) EXIT
        IF(ierror > 0) THEN
            WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(ursprung),'.xyz aufgetreten.'
            WRITE(*,*) 'Beende Einlesen !'
            STOP
        END IF
        if ( residue_number == resid ) THEN
                i=i+1
            IF(i==1) THEN
                IF ('THP'==TRIM(residue_name(1:3))) THEN        
                    P3HT=.TRUE.
                    Molname="P3HT"
                    IF(P3MT) Molname="P3MT"
                END IF
                
                IF( TRIM(residue_name)=='DIPBI' ) THEN
                    DIPBI=.TRUE.
                    Molname="DIPBI"
                    IF(DIPBI_KETTE) Molname="DIPBI_KETTE"
                END IF
                    
            END IF
        END IF
        IF (i==0 .and. j==Gesamtanzahl ) THEN
        write(*,*) ' Error: Es ist ein Fehler beim Einlesen ausgetreten'
        write(*,*) ' Vermutlich wurde nicht die entsprechende Resid-Nr gefunden',resid
        write(*,*) ' ENDE '
        STOP
        END IF  
    END DO
    CLOSE(78)
    
    MAXATM=i !Anzahl der Atomzentren in gro mit der resid
    
    CALL read_coord_from_gro_resid_to_array(atomsorteA,koordA,MAXATM,ursprung,resid)

END SUBROUTINE read_coord_from_gro_resid

SUBROUTINE read_coord_from_gro_resid_to_array(atomsorteA,koordA,MAXATM,ursprung,resid)
! Einlesen aus der .gro Datei und Umrechnung in .xyz Format (Ang)
INTEGER, INTENT (IN) :: MAXATM,resid
Character(3), ALLOCATABLE, Dimension(:), INTENT(OUT) :: atomsorteA
Real, ALLOCATABLE, Dimension(:,:), INTENT(OUT) :: koordA
Character(50), INTENT (IN) :: ursprung
Character(50)              ::dummyname
Character(5) :: residue_name
Character(3) :: atom
INTEGER::i,j,ierror,residue_number,atom_number,Gesamtanzahl
REAL::x,y,z,v_x,v_y,v_z

! Einlesen
OPEN(UNIT=79,FILE=TRIM(ursprung),STATUS='OLD',IOSTAT=ierror)
Read(79,*) dummyname                             !Einlesen des Infokopfes im File
Read(79,*) Gesamtanzahl                          !Einlesen der Atomanzahl im File

ALLOCATE(atomsorteA(MAXATM))
ALLOCATE(koordA(MAXATM,3))
i=0
DO j=1,Gesamtanzahl
    Read(79,'(i5,2a5,i5,3f8.3,3f8.4)',IOSTAT=ierror) residue_number,residue_name,atom,atom_number,x,y,z,v_x,v_y,v_z
    IF(ierror < 0) EXIT
    IF(ierror > 0) THEN
        WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(ursprung),'.xyz aufgetreten.'
        WRITE(*,*) 'Beende Einlesen !'
        STOP
    END IF
    if (residue_number == resid ) THEN
        i=i+1
        atomsorteA(i)=atom
        koordA(i,1)=x
        koordA(i,2)=y
        koordA(i,3)=z
    END IF
    IF( i == MAXATM ) EXIT
    IF (i==0 .and. j==Gesamtanzahl ) THEN
       write(*,*) ' Error: Es ist ein Fehler beim Einlesen ausgetreten'
       write(*,*) ' Vermutlich wurde nicht die entsprechende Resid-Nr gefunden',resid
       write(*,*) ' ENDE '
       STOP
   END IF  
END DO
CLOSE(79)

koordA(:,:)=koordA(:,:)*10.0 ! Umrechnung von nm (.gro-Format) nach Angstroem (.xyz-Format)

END SUBROUTINE read_coord_from_gro_resid_to_array

integer function Read_rows_file(filename,lio) ! Angabe mit Endung ; lio Parameter zur Vermeidung von Verwechslungen, falls mehrere Rechnungen gleichzeitig ablaufen.
    IMPLICIT NONE
    Character(10000), INTENT(IN) :: filename,lio
    Character(12000) :: bashline
    INTEGER::linenumber
            write(bashline,*) 'wc -l '//TRIM(filename)//' > wc_'//TRIM(lio)//'.txt'
            !Write(*,*) trim(bashline)
            CALL execute_command_line(TRIM(bashline)) 
            OPEN(unit=122,file='wc_'//TRIM(lio)//'.txt') 
            READ(122,*) linenumber
            CLOSE(122)
            CALL execute_command_line('rm wc_'//TRIM(lio)//'.txt') 
    Read_rows_file=linenumber
END function Read_rows_file

Real Function norm(atom1)
    IMPLICIT NONE
    Real, Dimension(3), Intent(IN) :: atom1
    norm = SQRT(atom1(1)**2 + atom1(2)**2 + atom1(3)**2)
END Function norm


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



character(len=5) function make_lio() ! make_lio create an string with RANDOM_NUMBER, in order_to_create_unique string
       ! Creates a random string unsing clock-time 
       IMPLICIT NONE
       !   "Convert an real to string."
       character(len=10) :: k
       INTEGER::i
       REAL::HARVEST
       
    integer, allocatable :: seed(:)
    integer:: size2

    CALL init_random_seed()

    call random_seed(size=size2)
    allocate(seed(size2))
    ! set seed(:) somehow
       call random_seed(put=seed)
       call system_clock(i)
       CALL random_seed(put=seed)
       call system_clock(i)
       CALL RANDOM_NUMBER(HARVEST)
       !WRITE(*,*) HARVEST,abs(mod(HARVEST,REAL(i)))*1.E5
       write (k,'(I5)') INT(abs(mod(HARVEST,REAL(i)))*1.E5)
       !WRITE(*,*) k
       make_lio=adjustl(TRIM(k))
end function make_lio

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


END PROGRAM  traj_path_projection
