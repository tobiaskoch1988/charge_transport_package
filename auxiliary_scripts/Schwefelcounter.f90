!Schwefel_counter_gro
PROGRAM Schwefel_counter_gro

call main()

CONTAINS

SUBROUTINE main
 Real, ALLOCATABLE, Dimension(:,:) :: koordA
Character(3), ALLOCATABLE, Dimension(:):: atomsorteA
Character(50) :: Molname
Character(50) :: ursprung,zielfile
INTEGER  :: resid
LOGICAL ::DIPBI,DIPBI_KETTE,P3HT,P3MT
Character(50) ::dummyname
Character(5) :: residue_name
Character(3) :: atom
INTEGER::i,j,ierror,residue_number,atom_number,MAXATM,Gesamtanzahl
INTEGER::residue_number_min,residue_number_max,last_residue_number
REAL::x,y,z,v_x,v_y,v_z

DIPBI_KETTE=.false.
P3MT=.false.

ursprung='indexed.gro'
zielfile='Resid_Molname_Schwefel.dat'
IF(command_argument_count() == 2) THEN         !Einlesen der Flags
	CALL get_command_argument(1,ursprung)
	CALL get_command_argument(2,zielfile)
	!CALL get_command_argument(3,calculation_method)
	!CALL get_command_argument(4,b)
    ELSE
        Write(*,*) 'ursprung.gro zieldateiname.dat '
        Write(*,*) 'example: equi_2000ps_500k_DIPBI.gro Resid_Molname_Schwefel.dat '
	STOP        
	Read(*,*) ursprung,residue_number_1
END IF




OPEN(UNIT=12,FILE=TRIM(zielfile),STATUS='REPLACE',IOSTAT=ierror,action='write')
WRITE(12,*) ' ### 1) Resid 2)Mol_type 3) N_Schwefel ( 33=DIPBI) aus: '//TRIM(ursprung)//'  '//TRIM(zielfile) 
OPEN(UNIT=77,FILE=TRIM(ursprung),STATUS='OLD',IOSTAT=ierror)
    Read(77,*) dummyname                             !Einlesen des Infokopfes im File
    Read(77,*) Gesamtanzahl                          !Einlesen der Atomanzahl im File
    Read(77,'(i5,2a5,i5,3f8.3,3f8.4)',IOSTAT=ierror) residue_number_min
    DO i=2,Gesamtanzahl
	Read(77,'(i5,2a5,i5,3f8.3,3f8.4)',IOSTAT=ierror) residue_number_max
    END DO
CLOSE(77)
WRITE(*,*) 'Grenzen der Resids in '//TRIM(ursprung)//' Resid_min / max:',residue_number_min,' / ',residue_number_max

i=0
N_Schwefel=0
last_residue_number=residue_number_min
OPEN(UNIT=77,FILE=TRIM(ursprung),STATUS='OLD',IOSTAT=ierror)
Read(77,*) dummyname                             !Einlesen des Infokopfes im File
Read(77,*) Gesamtanzahl                          !Einlesen der Atomanzahl im File
!write(*,*) 'Anzahl der Atome in .gro Datei ', Gesamtanzahl
Do j=1,Gesamtanzahl
        Read(77,'(i5,2a5,i5,3f8.3,3f8.4)',IOSTAT=ierror) residue_number,residue_name,atom,atom_number,x,y,z,v_x,v_y,v_z
        IF(ierror < 0) EXIT
        IF(ierror > 0) THEN
            WRITE(*,*) 'Error: Es ist eine Fehler beim Einlesen von ',TRIM(ursprung),'.gro aufgetreten.'
            WRITE(*,*) 'Beende Einlesen !'
            STOP
        END IF
	!WRITE(*,*) residue_number,last_residue_number,atom
	IF ( residue_number == last_residue_number+1 ) THEN
    		WRITE(*,*) (residue_number-1),TRIM(Molname),N_Schwefel
    		WRITE(12,*) (residue_number-1),TRIM(Molname),N_Schwefel
		last_residue_number=residue_number
		N_Schwefel=0
		i=0
	END IF 

        !if ( residue_number == resid ) THEN
            i=i+1
            if (TRIM(adjustl(atom))=='S') THEN
                N_Schwefel= N_Schwefel+1
            END IF
            IF(i==1) THEN
                IF ('THP'==TRIM(residue_name(1:3))) THEN        
                    P3HT=.TRUE.
                    Molname="P3HT"
                    IF(P3MT) Molname="P3MT"
                END IF
                
                IF( TRIM(residue_name)=='DIPBI' ) THEN
                    DIPBI=.TRUE.
                    Molname="DIPBI"
		    N_Schwefel=33
                    IF(DIPBI_KETTE) Molname="DIPBI_KETTE"
                END IF                
            END IF
        !END IF
	

        IF (i==0 .and. j==Gesamtanzahl ) THEN
            write(*,*) ' Error: Es ist ein Fehler beim Einlesen ausgetreten'
            write(*,*) ' Vermutlich wurde nicht die entsprechende Resid-Nr gefunden',resid
            write(*,*) ' ENDE '
            STOP
        END IF  
        
END Do
CLOSE(77)

WRITE(*,*) (residue_number),TRIM(Molname),N_Schwefel
WRITE(12,*) (residue_number),TRIM(Molname),N_Schwefel
CLOSE(12)
END SUBROUTINE main


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
    
END PROGRAM    
