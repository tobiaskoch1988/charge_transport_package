PROGRAM Partailladungen
IMPLICIT NONE
call main()

CONTAINS

SUBROUTINE main()
! Verwende Konvention erst geladene Molekuele mit Resid1 und Resid2
! Resid1 charged
! Resid2 charged
! Resid1 neutral
! Resid2 neutral
! Nachbarn neutral
! lambda_out=aussere Reoraganisationsenergie, c_p=Pekar factor, E_el_i= elektrostatische Energie(je nach Einlesefile fuer Elektronen oder Loecher)
! epsilon_opt= hight(optical) frequency permittivity
! epsilon_s= low (static) frequency permittivity
! epsilon_0= vacuum permittivitiy
! V3 mit screening function
! s = screening length in screening funciton epsilon(r)
! epsilon_screening{ r_(ai,bk) }=epsilon_s-(epsilon_s-1.0)*(1+s*r+0.5*(s**2)*(r**2))*exp(-s*r)   ! r Abstand 
IMPLICIT NONE
INTEGER::i,j,k,l,ierror,ierror1,Anzahl_Geometrien,N_coord_ges
REAL,DIMENSION(3)::summe_ai,summe_aj,diff_r
REAL, PARAMETER::PI=4.D0*DATAN(1.D0),Einheitenfaktor=1.6021766208E3 !umrechenfaktor elementarladung*10^22
REAL::lambda_out,dE_el_ij,E_el_i,E_el_j,epsilon_opt,epsilon_s,epsilon_0,c_p,r_cut
REAL::r, s,epsilon_screening,dE_el_ij_screening,E_el_i_screening,E_el_j_screening ! screening function variables
CHARACTER (LEN=10000) ::dummyline,foldername,zielfile
CHARACTER (LEN=3)::atom
CHARACTER (LEN=10)::N_Geo
Real, ALLOCATABLE, Dimension(:,:)  :: pair_xyz_q
INTEGER, ALLOCATABLE, Dimension(:) :: NAtoms_list,Start_index
REAL, ALLOCATABLE, Dimension(:) :: diff_q_ai,diff_q_aj
CHARACTER (LEN=10000)::xyz_q_filename,dummy1,lio
LOGICAL::Kugelvolumen,Debug=.false.,screening_function
REAL::x_box,y_box,z_box

! Startparameter
Kugelvolumen=.true.
screening_function=.true.
r_cut=15 ! Radius der Integrationskugel in Angström

x_box=139.0228  !in Angstroem
y_box=269.6753  !in Angstroem
z_box=145.8116  !in Angstroem

epsilon_opt=3.5 ! 2.75 ! Alq3
epsilon_s=3.0  !2.82   ! Alq3
epsilon_0=8.854187817620 !E-12 C/(V*m)

s=0.13 ! in Angstroem**-1 / screening laenge in Adrienko-Paper fur Alq3 (0.3 Ang**-1) in H20

c_p=(1/real(epsilon_opt)-1/real(epsilon_s)) ! Pekar factor
foldername='ordner'
xyz_q_filename='xyz_q_mulliken.xyz'
zielfile='results_lambda_out.dat'
Anzahl_Geometrien=9



IF(command_argument_count() == 5) THEN         !Einlesen der Flags
	CALL get_command_argument(1,xyz_q_filename)
    lio=TRIM(xyz_q_filename)  ! Parameter_um erstellte temporäre Dateien oder Ordner einzigartig zu benennen.
	CALL get_command_argument(2,N_Geo)
    CALL get_command_argument(3,dummy1)        ! in Angström!
    READ(dummy1,'(1F12.4)') r_cut
    CALL get_command_argument(4,foldername)
    CALL get_command_argument(5,zielfile)
    READ(N_Geo,'(I10)')  Anzahl_Geometrien 
ELSE IF (command_argument_count() == 1 ) THEN  ! Nutze Funktion zur Bestimmung der Anzahl der Geometrien	
    CALL get_command_argument(1,xyz_q_filename)
    Anzahl_Geometrien=Read_N_Geo(xyz_q_filename,lio)
ELSE 
    WRITE(*,*) "Usage: ./Partailladungen xyz_q_filename.xyz NGeo r_cut foldername zielfile"
    WRITE(*,*) "Bsp: ./Partailladungen xyz_q_mulliken.xyz 12 foldername zielfile"
    WRITE(*,*) "NGeo  == (Anzahl der Geometrien in der xyz_q*.xyz file)"
    WRITE(*,*) "r_cut == Cutoff-Radius der Nachbarn in der Kugel in Angröm!"
    STOP
END IF

IF(abs(r_cut)==0) THEN
    Kugelvolumen=.false.
    WRITE(*,*) 'Use box volume with:',x_box,y_box,z_box
END IF
WRITE(*,*) 'Anzahl der Geometrien',Anzahl_Geometrien
N_coord_ges=Read_rows_file(xyz_q_filename,lio)-2*Anzahl_Geometrien
WRITE(*,*) 'N_coord_ges',N_coord_ges
ALLOCATE(pair_xyz_q(N_coord_ges,4))
ALLOCATE(NAtoms_list(Anzahl_Geometrien),Start_index(Anzahl_Geometrien))

NAtoms_list=0
epsilon_screening=0
OPEN(UNIT=27,FILE=TRIM(xyz_q_filename),STATUS='OLD',IOSTAT=ierror)

k=0
WRITE(*,*) 'START EINLESEN'
DO i=1,Anzahl_Geometrien
    Read(27,*,IOSTAT=ierror) NAtoms_list(i)
    Start_index(i)=k+1   
    READ(27,*) dummyline
    DO j=1,NAtoms_list(i) 
        k=k+1
        Read(27,*,IOSTAT=ierror) atom,pair_xyz_q(k,1),pair_xyz_q(k,2),pair_xyz_q(k,3),pair_xyz_q(k,4)
        IF(ierror < 0) EXIT
        IF(ierror > 0) THEN
            WRITE(*,*) 'Es ist eine Fehler beim Einlesen von ',TRIM(xyz_q_filename),'.xyz aufgetreten.'
            WRITE(28,*) 'Fehler: Es_ist_eine_Fehler_beim_Einlesen_von_',TRIM(xyz_q_filename),'.xyz_aufgetreten._ ',TRIM(zielfile)
            WRITE(*,*) 'Beende Einlesen !'
            STOP
        END IF
    END DO
END DO

IF( Debug) THEN
    j=1
    WRITE(*,*) j, pair_xyz_q(j,1),pair_xyz_q(j,2),pair_xyz_q(j,3),pair_xyz_q(j,4)

    j=2
    WRITE(*,*) j, pair_xyz_q(j,1),pair_xyz_q(j,2),pair_xyz_q(j,3),pair_xyz_q(j,4)
    j=(N_coord_ges-1)
    WRITE(*,*) j, pair_xyz_q(j,1),pair_xyz_q(j,2),pair_xyz_q(j,3),pair_xyz_q(j,4)

    j=N_coord_ges
    WRITE(*,*) j, pair_xyz_q(j,1),pair_xyz_q(j,2),pair_xyz_q(j,3),pair_xyz_q(j,4)
END IF 

WRITE(*,*) 'START BERECHNUNG: lambda_out'


ALLOCATE(diff_q_ai(NAtoms_list(1)))
ALLOCATE(diff_q_aj(NAtoms_list(2)))
!Berechne Differenzen
DO i=0,(NAtoms_list(1)-1),1
    diff_q_ai(i+1)=pair_xyz_q(Start_index(1)+i,4)-pair_xyz_q(Start_index(3)+i,4)
END DO

DO i=0,(NAtoms_list(2)-1),1
    diff_q_aj(i+1)=pair_xyz_q(Start_index(4)+i,4)-pair_xyz_q(Start_index(2)+i,4)
END DO

lambda_out=0.0
E_el_i=0.0
E_el_j=0.0
epsilon_screening=1.0
E_el_i_screening=0.0
E_el_j_screening=0.0


DO k=Start_index(5),N_coord_ges,1 ! summe b_k alles Nachbaratome der Nachbarmolekuele in der Liste
    summe_ai=0.0
    summe_aj=0.0
    ! summe ueber a_i
    DO i=0,(NAtoms_list(1)-1),1  ! summe a_i
         DO l=1,3,1           ! diff_r=r_bk - r_ai
             diff_r(l)= pair_xyz_q(k,l)-pair_xyz_q(Start_index(1)+i,l)
         END DO
         summe_ai(:)=summe_ai(:)+diff_q_ai(i+1)*diff_r(:)/norm(diff_r(:))**3
         E_el_i=E_el_i+diff_q_ai(i+1)*pair_xyz_q(k,4)/norm(diff_r(:))
         
         IF (screening_function) THEN ! screening_function
            r=norm(diff_r(:)) ! abs( r_bk - r_ai)
            epsilon_screening=epsilon_s-(epsilon_s-1.0)*(1+s*r+0.5*(s**2)*(r**2))*exp(-1.0*s*r)
            E_el_i_screening=E_el_i_screening+diff_q_ai(i+1)*pair_xyz_q(k,4)/(r*epsilon_screening)
            IF ( i==1) THEN !.or. k >= 0.998*N_coord_ges ) THEN
                write(*,*) k, epsilon_screening,E_el_i,E_el_i_screening !r,(1+s*r+0.5*(s**2)*(r**2)),exp(-1.0*s*r)
            END IF
         END IF ! screening_function
    END DO !a_i

    DO j=0,(NAtoms_list(2)-1),1 ! summe a_j
        DO l=1,3,1              ! diff_r=r_bk - r_aj
            diff_r(l)= pair_xyz_q(k,l)-pair_xyz_q(Start_index(2)+j,l)
        END DO
        summe_aj(:)=summe_aj(:)+diff_q_aj(j+1)*diff_r(:)/norm(diff_r(:))**3
        E_el_j=E_el_j+diff_q_aj(j+1)*pair_xyz_q(k,4)/norm(diff_r(:))
        
        IF (screening_function) THEN ! screening_function
            r=norm(diff_r(:)) ! abs( r_bk - r_ai)
            epsilon_screening=epsilon_s-(epsilon_s-1.0)*(1+s*r+0.5*(s**2)*(r**2))*exp(-1.0*s*r)
            E_el_j_screening=E_el_j_screening+diff_q_aj(j+1)*pair_xyz_q(k,4)/(r*epsilon_screening)
        END IF ! screening_function    
    END DO !a_j
    lambda_out=lambda_out+norm(summe_ai(:)+summe_aj(:))**2
    IF( isnan(lambda_out)) THEN ! Test fuer Bug mit NaN
        WRITE(*,*) 'Fehler: NaN'
        WRITE(*,*) 'lambda_out',lambda_out,' in ',TRIM(xyz_q_filename)
        WRITE(*,*) ' norm(diff_r(:))', norm(diff_r(:))
        WRITE(*,*) ' summe_ai(:)',summe_ai(:)
        WRITE(*,*) ' summe_aj(:)',summe_aj(:)
        WRITE(*,*) 'diff_r(:)',diff_r(:)
        WRITE(*,*) 'E_el_i',E_el_i
        WRITE(*,*) 'E_el_j',E_el_j
        WRITE(*,*) 'FEHLER: k=',k
        STOP 
    END IF    
END DO ! k-summe  D_I(r_bk) -D_F(r_bk)

IF( Kugelvolumen) THEN
    lambda_out=lambda_out*c_p/(8.0*epsilon_0*PI)*(4.0*PI*r_cut**3/(3.0*(N_coord_ges-Start_index(5))) )  ! 2.Faktor als Kugelvolumen fuer die Integration 4Pi*r^3/3* 1/(Anzahl der Nachbarn/Integrationspunkte) = INT_V_{out} dV 
    lambda_out=lambda_out*Einheitenfaktor !Umrechnung auf eV
ELSE
    lambda_out=lambda_out*c_p/(8.0*epsilon_0*PI)*(x_box*y_box*z_box/(N_coord_ges-Start_index(5)) )  ! 2.Faktor als Boxvolumen  x_box*y_box*z_box /N_coord_ges-Start_index(5)  = INT_V_{out} dV  Mit dem Mittelwertsatz der Integralrechnung
    lambda_out=lambda_out*Einheitenfaktor !Umrechnung auf eV
END IF


dE_el_ij=(E_el_i-E_el_j)/(4.0*PI*epsilon_0*epsilon_s)
dE_el_ij=dE_el_ij*Einheitenfaktor  !Umrechnung auf eV

IF ( screening_function) THEN ! use screening_function epsilon_(r)
    dE_el_ij_screening=(E_el_i_screening-E_el_j_screening)/(4.0*PI*epsilon_0)
    dE_el_ij_screening=dE_el_ij_screening*Einheitenfaktor  !Umrechnung auf eV
END IF ! screening_function

WRITE(*,*) '****************** Results ******************'
write(*,*) 'Die outer Reorganisationsenergie ist:', lambda_out,' eV'
WRITE(*,*) 'dE_el_ij =',dE_el_ij,' eV'
IF ( screening_function) THEN ! use screening_function epsilon_(r)
    WRITE(*,*) 'dE_el_ij_screening =',dE_el_ij_screening,' eV'
    !write(*,*) E_el_i_screening,E_el_j_screening
END IF ! screening_function
WRITE(*,*) '*********************************************'

OPEN(UNIT=28,FILE=zielfile,STATUS='unknown',IOSTAT=ierror1,action='write',position="append")
    IF ( .NOT. screening_function) THEN
        write(28,'(1X,2F12.8,2X,A100)') lambda_out, dE_el_ij, TRIM(foldername)
    ELSE 
        write(28,'(1X,3F12.8,2X,A100)') lambda_out, dE_el_ij, dE_el_ij_screening, TRIM(foldername) 
    END IF 
CLOSE(28)
DEALLOCATE(diff_q_aj)
DEALLOCATE(diff_q_ai)
DEALLOCATE(Start_index,NAtoms_list)
DEALLOCATE(pair_xyz_q)
CLOSE(27)
END SUBROUTINE main

integer function Read_rows_file(filename,lio) ! Angabe mit Endung ; lio Parameter zur Vermeidung von Verwechslungen, falls mehrere Rechnungen gleichzeitig ablaufen.
IMPLICIT NONE
Character(10000), INTENT(IN) :: filename,lio
Character(12000) :: bashline
INTEGER::linenumber
        write(bashline,*) 'wc -l '//TRIM(filename)//' > wc_'//TRIM(lio)//'.txt'
        Write(*,*) trim(bashline)
        CALL execute_command_line(TRIM(bashline)) 
        OPEN(unit=122,file='wc_'//TRIM(lio)//'.txt') 
        READ(122,*) linenumber
        CLOSE(122)
        CALL execute_command_line('rm wc_'//TRIM(lio)//'.txt') 
Read_rows_file=linenumber
END function Read_rows_file

integer function Read_N_Geo(filename,lio) ! Funktion um Anzahl der Geometrien einzulesen aus xyz_q files
IMPLICIT NONE
Character(10000), INTENT(IN) :: filename,lio
Character(12000) :: bashline
INTEGER::linenumber
        write(bashline,*) ' grep Coord '//TRIM(filename)//' | wc -l > N_Geo_'//TRIM(lio)//'.txt'
        Write(*,*) trim(bashline)
        CALL execute_command_line(TRIM(bashline)) 
        OPEN(unit=123,file='N_Geo_'//TRIM(lio)//'.txt') 
        READ(123,*) linenumber
        CLOSE(123)
        CALL execute_command_line('rm N_Geo_'//TRIM(lio)//'.txt') 
        Read_N_Geo=linenumber
END function Read_N_Geo

Real Function norm(atom1)
        IMPLICIT NONE
        Real, Dimension(3), Intent(IN) :: atom1
        norm = SQRT(atom1(1)**2 + atom1(2)**2 + atom1(3)**2)
END Function norm

END PROGRAM Partailladungen
