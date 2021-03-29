Program Projection
!Program for calculating transfer integrals between two monomers. Execute 
!this program in a folder where ./dimer.log, ./fort.7, ../molA/monomer.log, 
!../molB/monomer.log, ../molA/fort.7, ../molB/fort.7 are available.
!It will print out certain transfer integrals (in eV!) and writes every integral of 
!interest in fock.txt

implicit none
!overlap-matrix
Real*8, Allocatable:: overlap(:,:), fock(:,:), transmat(:,:)
!expansioncoefficients of monomer A and B and Dimer D
Real*8, Allocatable:: zeta_A(:,:), zeta_B(:,:), zeta_D(:,:)
!dimer-eigenvalues
Real*8, Allocatable:: energy(:,:), zeit
!projectionmatrices 
Real*8, Allocatable:: gamma_A(:,:), gamma_B(:,:)
!transferintegral J_AB and effective transferintegral t_AB
Real*8, Allocatable::J_AB(:,:), e_A(:,:), e_B(:,:), S_AB(:,:)
Integer:: ierror, i, j, l, m, k, p,i_max,j_max
Integer:: n_A, n_B, n_D, counter, homoA, homoB, lumoA, lumoB, n_MO=0, n1_MO, n2_MO, n_MOmax , frontier,Maximal_MO

Character (len=8):: dummy
Character (len=30):: dummy2
Real*8:: dummy3, fock_max, S_AB_max
Real*8:: HOMO_A_to_HOMO_B,LUMO_A_to_LUMO_B,S_AB_HOMO_A_to_HOMO_B,S_AB_LUMO_A_to_LUMO_B
integer::i_homo,i_lumo
Character (len=100):: buf
Character (len=12) :: FM1='(A26,ES14.6)'         !Format output 1
Character (len=21) :: FM2='(A26,ES14.6,A,ES14.6)' !Format output 2
Character (len=29) :: FM3='(A8,ES14.6,X,ES14.6,X,ES14.6)' !Format output 3
!arrays for transformation-purposes 
Real*8, Allocatable:: eigval(:), S_diag(:,:), work(:)
Integer:: Info, lwork
!Zusatz fuer viele nMOs
REAL, ALLOCATABLE, DIMENSION (:)::E_Orbitale_A,E_Orbitale_B,w_A,w_B
!INTEGER, DIMENSION (17) :: Output_nMO=(/0,1,2,4,5,10,25,50,75,100,125,150,175,200,225,240,1000/) ! letztes Element wird spaeter aus n_MOmax gesetzt
INTEGER, DIMENSION (2) :: Output_nMO=(/0,300/) 
CHARACTER(30)::filename1
INTEGER::DummyNR1,DummyNR2
REAL::k_B,T,DummyNR3
REAL::homo_dim,lumo_dim
REAL::sumA_lo,sumB_lo,sumA_el,sumB_el
REAL::J_ges_lo_w1,J_ges_lo_w2,J_ges_el_w1,J_ges_el_w2,J_ges_lo_w3,J_ges_el_w3,J_ges_lo_w4,J_ges_el_w4
LOGICAL::calc_dipro_with_max_nMO=.false.,read_g09=.false.,read_dftb_plus=.false.,DEBUG=.false.
! additional data for dftb+
INTEGER, ALLOCATABLE, DIMENSION(:,:) ::IAtom_Nneigh_Norbs
INTEGER::NAtoms,ATOM1,INEIGH,IATOM2F,icase,i_min,j_min,Ineighbours,Ineighbours_max
CHARACTER(500)::method
INTEGER, DIMENSION(3)::CELL
REAL::E_deg_cut_off


E_deg_cut_off=0.01 ! eV energy cutoff to treat the orbitals with degenerated states

IF(command_argument_count() == 1) THEN         !Einlesen der Flags
    CALL get_command_argument(1,method)
    IF( (TRIM(method) == 'g09') .OR. (TRIM(method) == 'sub_g09') )THEN
        read_g09=.true.
        read_dftb_plus=.false.
    ELSE IF((TRIM(method) == 'dftb+') .OR. (TRIM(method) == 'sub_dftb_plus') .OR. (TRIM(method) == 'sub_dftb+') )THEN
        read_dftb_plus=.true.
        read_g09=.false.
    ELSE
        WRITE(*,*) 'Error: select g09 or dftb+'
        STOP
    END IF
ELSE
        read_g09=.true. !default is g09
END IF

OPEN(Unit=23, File='Daten-V10.txt', Status='REPLACE', IOSTAT=ierror)
do p=1,size(Output_nMO)
        n_MO=Output_nMO(p)
        write(*,*) "Berechne Schritt fuer n_MO=",n_MO
        
    IF( read_g09 ) THEN ! read gausian09 output (default setting)
        Open(Unit=7, File='dimer.log', Status='OLD', IOSTAT=ierror)
        
        Open(Unit=8, File='../molA/fort.7', Status='OLD', IOSTAT=ierror)
        
        Open(Unit=9, File='../molB/fort.7', Status='OLD', IOSTAT=ierror)
        
        Open(Unit=10, File='fort.7', Status='OLD', IOSTAT=ierror)
        
        OPEN(Unit=11, File='fock.txt', Status='REPLACE', IOSTAT=ierror)
        
        
        !Reading number of basis functions for molecule A
        Open(Unit=12, File="../molA/monomer.log", Status='OLD', IOSTAT=ierror)
        
        do
            read(12,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            if (index(buf,'NBasis=') .ne. 0) then
                Read(dummy,*) n_A
                exit
            end if
        
        end do
        close(12)
        
        !looking for the indices of HOMO and LUMO
        Open(Unit=12, File="../molA/monomer.log", Status='OLD', IOSTAT=ierror)
        
        do 
                read(12,*,IOSTAT=ierror) dummy, buf
                if(ierror/=0) exit
                if(index(buf,'alpha') .ne. 0) then
                        Read(dummy,*) homoA
                        lumoA=homoA+1
                        exit
                end if
        end do
        
        !extra: Read orbital energies to construct weighting function.
        Allocate(E_Orbitale_A(n_A))
        E_Orbitale_A(:)=0.0d0
        do ! Einlesen der Orbitalenergien in umgekehrter Reihenfolge; hinterher aufsteigend in E_Orbitale_A
            read(12,*,IOSTAT=ierror) buf,dummy
            if(ierror/=0) exit
            if (index(buf,'Orbital') .ne. 0 .AND. index(dummy,'energies') .ne. 0) then
                Read(12,*) DummyNR1,DummyNR2
                DO i=1,n_A
                    read(12,*,IOSTAT=ierror) DummyNR1,dummy,E_Orbitale_A(i),DummyNR3
                    if(ierror/=0) exit
                END DO
                exit
            end if
        end do
        
        close(12)
        
        !Reading number of basis functions for molecule B
        Open(Unit=13, File="../molB/monomer.log", Status='OLD', IOSTAT=ierror)
        
        do
            read(13,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            if (index(buf,'NBasis=') .ne. 0) then
                Read(dummy,*) n_B
                exit
            end if
        end do
        
        !extra: Read orbital energies to construct weighting function.
        Allocate(E_Orbitale_B(n_B))
        E_Orbitale_B(:)=0.0d0
        do ! Einlesen der Orbitalenergien in umgekehrter Reihenfolge; hinterher aufsteigend in E_Orbitale_B
            read(13,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            if (index(buf,'Orbital') .ne. 0 .AND. index(dummy,'energies') .ne. 0) then
                Read(13,*) DummyNR1,DummyNR2
                DO i=1,n_B,1
                    read(13,*,IOSTAT=ierror) DummyNR1,dummy,E_Orbitale_B(i),DummyNR3   
                    !WRITE(*,*) E_Orbitale_A(i),E_Orbitale_B(i)
                    if(ierror/=0) exit
                END DO
                exit
            end if
        end do
        
        close(13)
        
        !looking for the indices of HOMO and LUMO
        Open(Unit=13, File="../molB/monomer.log", Status='OLD', IOSTAT=ierror)
        
        do 
                read(13,*,IOSTAT=ierror) dummy, buf
                if(ierror/=0) exit
                if(index(buf,'alpha') .ne. 0) then
                        Read(dummy,*) homoB
                        lumoB=homoB+1
                        exit
                end if
        end do
        close(13)
        !Reading number of basis functions for the dimer
        do
            read(7,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            if (index(buf,'NBasis=') .ne. 0) then
                Read(dummy,*) n_D
                exit
            end if
        end do
        close(7)
        
        
    ELSE IF (read_dftb_plus) THEN ! read data from dftb+ (formated) output
    
        IF (file_exists('eigenvec_dim.dat')) THEN
            Open(Unit=7, File='eigenvec_dim.dat', Status='OLD', IOSTAT=ierror)
        END IF 
        IF (file_exists('../molA/eigenvec_molA.dat')) THEN
            Open(Unit=8, File='../molA/eigenvec_molA.dat', Status='OLD', IOSTAT=ierror)
        END IF 
        IF (file_exists('../molB/eigenvec_molB.dat')) THEN
            Open(Unit=9, File='../molB/eigenvec_molB.dat', Status='OLD', IOSTAT=ierror)
        END IF         
        
        OPEN(Unit=11, File='fock.txt', Status='REPLACE', IOSTAT=ierror)
        
       !Reading number of basis functions for molecule A       
        do
            read(8,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            IF(index(dummy,'NBasis=') .ne. 0) then ! dummy = nbasis, for an uneven number of basis functions
                backspace(8) ! back 
                CYCLE
            END IF
            if (index(buf,'NBasis=') .ne. 0) then ! read number of basis function in molA
                Read(dummy,*) n_A
                read(8,*,IOSTAT=ierror) buf, dummy ! Read homo index
                if(ierror/=0) stop "Error: reading eigenvector molA."
                if (index(buf,'HOMO=') .ne. 0) then 
                    Read(dummy,*) homoA
                    lumoA=homoA+1
                END IF ! read HOMO
                exit
            end if
        end do ! reading: n_A,homoA,lumoA
        close(8)
        !WRITE(*,*) 'molA:',n_A,homoA,lumoA
        Open(Unit=8, File='../molA/eigenvec_molA.dat', Status='OLD', IOSTAT=ierror) ! open again and read orbital energies. (as the order in the file is mixed)
        Allocate(E_Orbitale_A(n_A))
        E_Orbitale_A=0.0d0
        do ! Einlesen der Orbitalenergien in umgekehrter Reihenfolge; hinterher aufsteigend in E_Orbitale_A
            read(8,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            IF(index(dummy,'Eigenvalues') .ne. 0) then 
                backspace(8) ! back 
                CYCLE
            END IF
            if (index(buf,'Eigenvalues') .ne. 0 .AND. index(dummy,'[a.u.]') .ne. 0) then
                !Read(8,*,IOSTAT=ierror) dummy
                DO i=1,n_A,1
                    read(8,*,IOSTAT=ierror) E_Orbitale_A(i)
                    if(ierror/=0) stop "Error: reading eigenvector molA."
                END DO
                exit
            end if
        end do ! read orbital energies in Hartee molA
        !Do i=1,n_A
        !   WRITE(*,*)  E_Orbitale_A(i)
        !END Do 
        WRITE(*,*) 'molA:',n_A,homoA,lumoA

       !Reading number of basis functions for molecule B       
        do
            read(9,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            IF(index(dummy,'NBasis=') .ne. 0) then ! dummy = nbasis, for an uneven number of basis functions
                backspace(9) ! back 
                CYCLE
            END IF
            if (index(buf,'NBasis=') .ne. 0) then ! read number of basis function in molB
                Read(dummy,*) n_B
                read(9,*,IOSTAT=ierror) buf, dummy ! Read homo index
                if(ierror/=0) stop "Error: reading eigenvector molB."
                if (index(buf,'HOMO=') .ne. 0) then 
                    Read(dummy,*) homoB
                    lumoB=homoB+1
                END IF ! read HOMO
                exit
            end if
        end do ! reading: n_B,homoB,lumoB
        close(9)
        !WRITE(*,*) 'molB:',n_B,homoB,lumoB
        Open(Unit=9, File='../molB/eigenvec_molB.dat', Status='OLD', IOSTAT=ierror) ! open again and read orbital energies. (as the order in the file is mixed)
        Allocate(E_Orbitale_B(n_B))
        E_Orbitale_B=0.0d0
        do ! Einlesen der Orbitalenergien in umgekehrter Reihenfolge; hinterher aufsteigend in E_Orbitale_B
            read(9,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            IF(index(dummy,'Eigenvalues') .ne. 0) then 
                backspace(9) ! back 
                CYCLE
            END IF
            if (index(buf,'Eigenvalues') .ne. 0 .AND. index(dummy,'[a.u.]') .ne. 0) then
                !Read(9,*,IOSTAT=ierror) dummy
                DO i=1,n_B,1
                    read(9,*,IOSTAT=ierror) E_Orbitale_B(i)
                    if(ierror/=0) stop "Error: reading eigenvector molB."
                END DO
                exit
            end if
        end do ! read orbital energies in Hartee molB
        !Do i=1,n_B
        !   WRITE(*,*)  E_Orbitale_B(i)
        !END Do 
        WRITE(*,*) 'molB:',n_B,homoB,lumoB
        IF ( n_A == 0 .OR. n_B ==0 ) THEN
                WRITE(*,*) 'Error: Number of basis functions is zero, check input:',n_A,n_B
                STOP
        END IF

       !Reading number of basis functions for dimer: n_D      
        do
            read(7,*,IOSTAT=ierror) buf, dummy
            if(ierror/=0) exit
            IF(index(dummy,'NBasis=') .ne. 0) then ! dummy = nbasis, for an uneven number of basis functions
                backspace(7) ! back 
                CYCLE
            END IF
            if (index(buf,'NBasis=') .ne. 0) then ! read number of basis function in dimer
                Read(dummy,*) n_D
                read(7,*,IOSTAT=ierror) buf, dummy ! Read homo index
                if(ierror/=0) stop "Error: reading eigenvector dimer"
                if (index(buf,'HOMO=') .ne. 0) then 
                    Read(dummy,*) homo_dim
                    lumo_dim=homo_dim+1
                END IF ! read HOMO
                exit
            end if
        end do ! reading: n_D
        CLOSE(7)
    ELSE
        WRITE(*,*) 'Error: You need to select an input format for DIPRO calcualtion (g09 or DFTB+)'
        STOP
    END IF ! read g09 or dftb+

!allocating arrays with known number of orbitals
Allocate(energy(n_D,n_D))	
Allocate(overlap(n_D,n_D))
Allocate(zeta_A(n_A,n_D))
Allocate(zeta_B(n_B,n_D))
Allocate(zeta_D(n_D,n_D))
Allocate(gamma_A(n_A,n_D))
Allocate(gamma_B(n_B,n_D))
Allocate(J_AB(n_A,n_B))
Allocate(S_AB(n_A,n_B))
Allocate(e_A(n_A,n_A))
Allocate(e_B(n_B,n_B))
Allocate(fock(n_D,n_D))
Allocate(transmat(n_D,n_D))


energy=0.0d0
zeta_A=0.0d0
zeta_B=0.0d0

IF( read_g09 ) THEN ! read gausian09 output (default setting)
    !Read overlap
    Open(Unit=7, File='dimer.log', Status='OLD', IOSTAT=ierror)
    do
        Read (7,'(A)',IOSTAT=ierror) buf
        if (ierror /= 0) stop "Error reading overlap matrix."
        !reading the file until the Overlap is shown
        if ( index(buf, '* Overlap') .ne. 0) exit
    end do
    
    l=ceiling(real(n_D)/real(5))
    Do i=1,l
    Read(7,'(I17)') j
    
        Do m=j,n_D
            if (m<=j+4) then
                Read(7,*) dummy, overlap(m, j:m)
            else
                Read(7,*) dummy, overlap(m,j:(j+4))
            End if
        End Do
    End Do
    
    !overlap is symmetric, so mirror the lower triangle-matrix
    Do i=1,n_D
        overlap(i,i:n_D)=overlap(i:n_D,i)
    end Do
    
    
    
    j=ceiling(real(n_A)/real(5))
    !Reading coefficents of monomer A out of fort.7
    !first line is FORMAT of coefficients
    Read(8,'(A8)') dummy
    counter=0
    Do i=1,n_A
        Read(8,'(A18, D15.8)') dummy2, energy(i,i)
        
        Do l=1,j
            if (l/=j) Then
                Read(8,dummy) zeta_A(i,(5*(l-1)+1):(5*l))
            Else
                Read(8,dummy) zeta_A(i,(j*5-5+1):(n_A))
            End if
        End Do
    End Do
    
    
    j=ceiling(real(n_B)/real(5))
    !Reading coefficents of monomer B out of fort.7
    Read(9,'(A8)') dummy
    counter=0
    Do i=1,n_B
        Read(9,'(A18, D15.8)') dummy2, energy(i,i)
        Do l=1,j
            if (l/=j) Then
    
                Read(9,dummy) zeta_B(i,n_A+(5*(l-1)+1):n_A+(5*l))
            Else
                Read(9,dummy) zeta_B(i,n_A+(j*5-5+1):n_A+(n_B))
    
            End if
        End Do
    End Do
    
    
    !Reading coefficients and eigenvalues of dimer out of fort.7
    Read(10,'(A8)') dummy
    
    Do i=1,n_D
    Read(10,'(A18, D15.8)') dummy2, energy(i,i)
    j=ceiling(real(n_D)/real(5))
    
        Do l=1,j
            if (l/=j) Then
    
            Read(10,dummy) zeta_D(i,(5*(l-1)+1):(5*l))
        Else
            Read(10,dummy) zeta_D(i,(j*5-5+1):(n_D))
    
            End if
        End Do
    End Do
    CLOSE(10)
    CLOSE(9)
    CLOSE(8)
    CLOSE(7)
ELSE IF (read_dftb_plus) THEN ! read data from dftb+ (formated) output
    !!!!!! Reading DFTB+ data !!!!!!!!!!!!!!!!!!!!
    ! build 1-Matrix
    overlap=0.0d0 
    DO i=1,n_D
        overlap(i,i)=1.0d0
    END DO
    !!! reading overlap S_AB from DFTB+ format for sparse storage S_AB-Matrix (overreal.dat WriteRealHS =Yes) !!!
    IF (file_exists('./overreal.dat')) THEN
        Open(Unit=10, File='./overreal.dat', Status='OLD', IOSTAT=ierror)
    END IF ! overlap matrix in real format exists? 
        
    icase=1 ! select cases
    Ineighbours_max=1E5 !Maximale Anzahl an besuchten Nachbarn (Wird spaeter als sum_i(Nneigh_i) gesetzt.
    Ineighbours=0 ! counts the number of neighbours (gives the loop exit)
    do    ! read data from overreal.dat (DFTB+ real overlap matrix S) 
        Read (10,'(A)',IOSTAT=ierror) buf ! read #-lines
        !WRITE(*,*) TRIM(buf)
        if (index(buf,'#') .ne. 0) then 
            CONTINUE
        END IF
        if (ierror /= 0) THEN
            WRITE(*,*) "Error reading overlap matrix."
            stop 
        END IF
        IF (icase==1) THEN ! read NAtoms
            Read (10,*,IOSTAT=ierror) NAtoms
            if (ierror /= 0) stop "Error reading overlap matrix."
            WRITE(*,*) 'NAtoms=',NAtoms
            icase=icase+1
            CYCLE
        END IF ! read NAtoms
        
        IF(icase==2) THEN
            Allocate(IAtom_Nneigh_Norbs(NAtoms,4)) ! Number of Atoms, Number of Neighbour including the atom I, Number of orbitals Norbs, 4= sum of N_Orbs to the index
            IAtom_Nneigh_Norbs=0 ! Initializing
            Do j=1,NAtoms
                Read (10,*,IOSTAT=ierror) IAtom_Nneigh_Norbs(j,1),IAtom_Nneigh_Norbs(j,2),IAtom_Nneigh_Norbs(j,3)
                !WRITE(*,*) IAtom_Nneigh_Norbs(j,1),IAtom_Nneigh_Norbs(j,2),IAtom_Nneigh_Norbs(j,3)
                IF (j==1) THEN
                    IAtom_Nneigh_Norbs(j,4)=1
                ELSE
                    IAtom_Nneigh_Norbs(j,4)=IAtom_Nneigh_Norbs(j-1,4)+IAtom_Nneigh_Norbs(j-1,3)  !
                END IF
                if (ierror /= 0) THEN
                    WRITE(*,*) "Error reading overlap matrix."
                    stop 
                END IF ! ierror
            END DO ! read 
            Ineighbours_max=sum(IAtom_Nneigh_Norbs(:,2))
            WRITE(*,*) 'Ineighbours_max=',Ineighbours_max
            icase=icase+1
            CYCLE
        END IF ! icase=2 read IAtom_Nneigh_Norbs
        
        IF(icase==3) THEN ! read: IATOM1    INEIGH   IATOM2F  ICELL(1)  ICELL(2)  ICELL(3)
            Read (10,*,IOSTAT=ierror) ATOM1,INEIGH,IATOM2F,CELL(1),CELL(2),CELL(3)
            !WRITE(*,*) ATOM1,INEIGH,IATOM2F,CELL(1),CELL(2),CELL(3)
            if (ierror /= 0) stop "Error reading overlap matrix."
            icase=4
            CYCLE
        END IF ! case 3
        
        IF(icase==4) THEN ! read overlapp matrix part
            !WRITE(*,*) ATOM1,INEIGH,IATOM2F,CELL(1),CELL(2),CELL(3)
            i_min= IAtom_Nneigh_Norbs(ATOM1,4)             !SUM(IAtom_Nneigh_Norbs(1:(ATOM1),3))+1 ! sum_1^ATOM1 NOrbs + 1
            IF (ATOM1==NAtoms) THEN ! check end of array, array will get out of range otherwise 
                i_max=IAtom_Nneigh_Norbs(ATOM1,4)+IAtom_Nneigh_Norbs(ATOM1,3)-1
            ELSE
                i_max=IAtom_Nneigh_Norbs(ATOM1+1,4)-1
            END IF
            !WRITE(*,*) 'SUM_ORBS i_min,i_max:',i_min,i_max
            j_min=IAtom_Nneigh_Norbs(IATOM2F,4)        !SUM(IAtom_Nneigh_Norbs(1:(IATOM2F),3))+1 ! sum_1^AATOMF NOrbs + 1
            IF(IATOM2F==NAtoms) THEN
                j_max=IAtom_Nneigh_Norbs(IATOM2F,4)+IAtom_Nneigh_Norbs(IATOM2F,3)-1
            ELSE
                j_max=IAtom_Nneigh_Norbs(IATOM2F+1,4)-1
            END IF             
            !WRITE(*,*) 'SUM_ORBS j_min,j_max:',j_min,j_max
            DO i=i_min,i_max,1
                Read (10,*,IOSTAT=ierror) overlap(i,j_min:j_max)
                if (ierror /= 0) stop "Error reading overlap matrix."
                !WRITE(*,*) overlap(i,j_min:j_max)
            END DO
            icase=3
            Ineighbours=Ineighbours+1
            if ( Ineighbours == Ineighbours_max) exit  !reading the file until the maximum number of neighbours is reached: Ineighbours_max
            CYCLE
        END IF ! case 4        
    end do ! read data from overreal.dat 
    CLOSE(10)
    !!! end reading overlap S_AB from DFTB+ format !!!
     WRITE(*,*) 'Ineighbours=', Ineighbours

    !overlap is symmetric, so mirror the lower triangle-matrix
    Do i=1,n_D
        overlap(i,i:n_D)=overlap(i:n_D,i)
    end Do
    
    WRITE(*,*) '________End reading dftb+ overlap_______'
    
    !!!Start reading coefficents of monomer A from dftb+:     
    i=0 !count eigenvectors
    do ! Einlesen eigenvectoren monomer A
        read(8,*,IOSTAT=ierror) buf
        !WRITE(*,*) buf
        if(ierror/=0) EXIT
        if (index(buf,'Eigenvector:') .ne. 0 ) then
            i=i+1
            !Read(8,*,IOSTAT=ierror) dummy
            DO j=1,n_A,1 !read eigenvector
                read(8,*,IOSTAT=ierror) zeta_A(i,j)
                !WRITE(*,*) zeta_A(i,j)
                if(ierror/=0) STOP "Error: reading eigenvectors: eigenvec_molA.dat"
            END DO
            IF (i== n_A) EXIT
        end if
    end do ! Einlesen eigenvectoren
    CLOSE(8)

    WRITE(*,*) 'n_A=',n_A

    
    !!!Start reading coefficents of monomer B from dftb+:  
    i=0 !count eigenvectors
    do   ! Einlesen eigenvectoren monomer B
        read(9,*,IOSTAT=ierror) buf
        !WRITE(*,*) buf
        if(ierror/=0) EXIT
        if (index(buf,'Eigenvector:') .ne. 0 ) then
            i=i+1
            !Read(9,*,IOSTAT=ierror) dummy
            DO j=1,n_B,1 !read eigenvector
                read(9,*,IOSTAT=ierror) zeta_B(i,n_A+j)
                !WRITE(*,*) zeta_B(i,n_A+j)
                if(ierror/=0) STOP "Error: reading eigenvectors: eigenvec_molB.dat"
            END DO
            IF (i== n_B) EXIT
        end if
    end do ! Einlesen eigenvectoren monomer B 
    CLOSE(9)
    WRITE(*,*) 'n_B=',n_B
  
    
    !!!Start reading dimer data from dftb+
    IF (file_exists('eigenvec_dim.dat')) THEN !open again
            Open(Unit=7, File='eigenvec_dim.dat', Status='OLD', IOSTAT=ierror)
    END IF 
    do ! Einlesen der Orbitalenergien in energy(i,i) (in Hartree)
            read(7,*,IOSTAT=ierror) buf, dummy
            !WRITE(*,*) buf, dummy
            if(ierror/=0) exit

            if (index(buf,'Eigenvalues') .ne. 0 .AND. index(dummy,'[a.u.]') .ne. 0) then
                DO i=1,n_D,1
                    read(7,*,IOSTAT=ierror) energy(i,i)
                    if(ierror/=0) STOP "Error: reading eigenvectors: eigenvec_dim.dat"
                END DO
                exit
            end if
    end do ! read orbital energies in Hartee molA
    !Do i=1,n_D
    !    WRITE(*,*) energy(i,i)
    !END Do 
    !WRITE(*,*) 'n_D=',n_D
       
    !!!Start reading coefficents of the dimer from dftb+: 
    i=0 !count eigenvectors
    do   ! Einlesen eigenvectoren dimer
        read(7,*,IOSTAT=ierror) buf
        !WRITE(*,*) buf
        if(ierror/=0) EXIT
        if (index(buf,'Eigenvector:') .ne. 0 ) then
            i=i+1
            !Read(7,*,IOSTAT=ierror) dummy
            DO j=1,n_D,1 !read eigenvector
                read(7,*,IOSTAT=ierror) zeta_D(i,j)
                !WRITE(*,*) zeta_D(i,j)
                if(ierror/=0) STOP "Error: reading eigenvectors: eigenvec_dim.dat"
            END DO
            IF (i== n_D) EXIT
        end if
    end do ! Einlesen eigenvectoren dimer 
    CLOSE(7)
    !WRITE(*,*) 'n_D=',n_D
    !WRITE(*,*) 'Ineighbours=', Ineighbours
        
    
ELSE ! end read dftb+ data
        WRITE(*,*) 'Error: You need to select an input format for DIPRO calcualtion (g09 or DFTB+)'
        STOP
END IF ! read g09 or dftb+


WRITE(*,*) 'START DIPRO METHODE'

  
!!!!!!! Start with the DIPRO calculations
!calculating gamma
gamma_A=matmul(zeta_A,(matmul(overlap,transpose(zeta_D))))
gamma_B=matmul(zeta_B,(matmul(overlap,transpose(zeta_D))))

!calculation of J_AB
J_AB=matmul(gamma_A,matmul(energy,transpose(gamma_B)))

!calculation of side-energies
e_A=matmul(gamma_A,matmul(energy,transpose(gamma_A)))
e_B=matmul(gamma_B,matmul(energy,transpose(gamma_B)))

!calculation of overlap of monomers					!80
S_AB=matmul(gamma_A,transpose(gamma_B))

!forming the complete fock-matrix
fock(1:n_A,1:n_A)=e_A
fock(n_A+1:n_D,n_A+1:n_D)=e_B
fock(1:n_A,n_A+1:n_D)=J_AB
fock(n_A+1:n_D,1:n_A)=transpose(J_AB)

!calculating the overlap-matrix for transformation to effective fock-matrix
transmat(1:n_A,1:n_A)=matmul(gamma_A,transpose(gamma_A))
transmat(n_A+1:n_D,n_A+1:n_D)=matmul(gamma_B,transpose(gamma_B))
transmat(1:n_A,n_A+1:n_D)=S_AB
transmat(n_A+1:n_D,1:n_A)=transpose(S_AB)

!converting Hartree to eV!
fock=fock*27.21138386d0
  
write (filename1,'(I0)') n_MO
Open(Unit=14, File='Ergebnis-'//TRIM(filename1)//'MO-V10.txt' , Status='REPLACE', IOSTAT=ierror)
Open(Unit=15, File='Matrix_'//TRIM(filename1)//'MO-V10.txt' , Status='REPLACE', IOSTAT=ierror)
Write(14,*) 'HomoA:', homoA
Write(14,*) 'LumoA:', lumoA
Write(14,*) 'HomoB:', homoB
Write(14,*) 'LumoB:', lumoB
!printing certain matrix-elements to check validity
		Write(14,*) 'untransformed:'
		Write(14,FM1) 'HOMO_A to HOMO_A', fock(homoA,homoA)
		Write(14,FM1) 'HOMO_A to LUMO_A', fock(homoA,lumoA)
		Write(14,FM1) 'HOMO_A to HOMO_B', fock(homoA,homoB+n_A)
		Write(14,FM1) 'HOMO_A to LUMO_B', fock(homoA,lumoB+n_A)
		Write(14,FM1) 'LUMO_A to LUMO_B', fock(lumoA,lumoB+n_A)

!Writing out lower left part of the fock-matrix (monomerA to monomerB-transfer)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Deallocate(energy)
Deallocate(zeta_A)
Deallocate(zeta_B)
Deallocate(zeta_D)
IF (read_dftb_plus) THEN
    Deallocate(IAtom_Nneigh_Norbs)
END IF 

!Now we take just the elements we need for HOMO and LUMO +/-n_MO
!while taking energy as a dummy-matrix
!Determine the maximum Value for n_MO as the smallest gap for the frontier space!
if( homoA < n_A/2 ) then ! determine n1_MO for molA
	n1_MO=homoA-1
else
	n1_MO=n_A-homoA-1
end if

if( homoB < n_B/2 ) then ! determine n2_MO for molB
        n2_MO=homoB-1
else
        n2_MO=n_B-homoB-1
end if

IF(n1_MO < n2_MO) THEN ! select smaller range as n_MO
        n_MO=n1_MO
ELSE
        n_MO=n2_MO
END IF

IF(n_MO < 0) THEN
     WRITE(*,*) ' n_MO < 0; Error in evaluation of n_MO. Check selected space of frontier orbitals.'
     WRITE(*,*) ' STOP DIPRO program'
     STOP
END IF
Maximal_MO=n_MO !set last element in Output_nMO to the maximum nMO value
Output_nMO(size(Output_nMO))=Maximal_MO
n_MO=Output_nMO(p)                ! calculate iteration p ! selects n_MO here!!!!
!WRITE(*,*) "Gebe die Anzahl der Orbitale zusaetzlich zum HOMO an von 0 bis ",n_MO
!READ(*,*) n_MO
WRITE(*,*) "Anzahl der HOMO-", n_MO, " und Anzahl der LUMO+",n_MO
IF( (n_MO .GT. Maximal_MO) )THEN
    WRITE(*,*) 'Error: The_selected n_MO=',n_MO,' is too big'
    WRITE(*,*) 'The limit is:',Maximal_MO
    IF(.not. calc_dipro_with_max_nMO) THEN ! calculate with 
        n_MO=Maximal_MO
        calc_dipro_with_max_nMO=.true.
        WRITE(*,*) 'Selected nMO limit:',n_MO
    ELSE
        WRITE(*,*) 'STOP DIPRO program'
        STOP
    END IF
END IF
n_MOmax=4*(n_MO+1) ! maximum dimension of matrix

Allocate(energy(n_MOmax,n_MOmax))
energy(1:(n_MOmax/2),1:(n_MOmax/2))=transmat(homoA-n_MO:lumoA+n_MO,homoA-n_MO:lumoA+n_MO)
energy((n_MOmax/2+1):n_MOmax,(n_MOmax/2+1):n_MOmax)=transmat((n_A+homoB)-n_MO:(n_A+lumoB)+n_MO,(n_A+homoB)-n_MO:(n_A+lumoB)+n_MO)
energy((n_MOmax/2+1):n_MOmax,1:(n_MOmax/2))=transmat((n_A+homoB)-n_MO:(n_A+lumoB)+n_MO,homoA-n_MO:lumoA+n_MO)
energy(1:(n_MOmax/2),(n_MOmax/2+1):n_MOmax)=transpose(energy((n_MOmax/2+1):n_MOmax,1:(n_MOmax/2+1)))
Deallocate(transmat)
Allocate(transmat(n_MOmax,n_MOmax))
transmat=energy

! Now the diagonalization of S and transformation to effective fock-matrix will be done
lwork=3*n_MOmax
Allocate(eigval(n_MOmax))
Allocate(work(lwork))
Allocate(S_diag(n_MOmax,n_MOmax))

!Initializing the diagonal-matrix and matrix to be diagonalized

CALL DSYEV('V','U', n_MOmax , transmat, n_MOmax , eigval, work, lwork, info)

IF ( info == 0) then
        Write(*,*) ' Diagonalisation successful '
        else if ( info < 0) then 
                Write(*,*) ' The argument',info,' had an illigal value'        
        else
                Write(*,*) info, ' > 0 => the algorithm failed to converge'
                Write(*,*) 'off-diagonal elements of an intermediate tridiagonal form did not converge to zero.'
END IF     


!! check for small eigenvalues, to prevent division by zero in S_diag calculation, which leads to NaN for alle Elements in transmat. 
j=0
Do i=1,n_MOmax
	IF ( abs(eigval(i)) < 1.0E-11 ) THEN
		IF(j==0) THEN ! Warning is only printed once.		
			write(*,*) ' Caution, small eigenvalue in eigval(i) encountered, check the results: ',eigval(i),' i=',i
			write(*,*) ' use: ',sign(1.0D-11,eigval(i))
			j=1
		END IF 
		eigval(i)=sign(1.0D-11,eigval(i))
	END IF 
end do ! check small eigenvalues

S_diag=0.0d0
Do i=1,n_MOmax
	S_diag(i,i)=1.0d0/sqrt(eigval(i))
End Do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Having calculated all the eigenvalues and eigenvectors, one can finally build the transformation-matrix


transmat=matmul(transmat,matmul(S_diag,transpose(transmat)))


!finally calculating the effective fock-matrix:
!takiing the relevant elements of fock matrix analogous to transmat
energy(1:(n_MOmax/2),1:(n_MOmax/2))=fock(homoA-n_MO:lumoA+n_MO,homoA-n_MO:lumoA+n_MO)
energy((n_MOmax/2+1):n_MOmax,(n_MOmax/2+1):n_MOmax)=fock((n_A+homoB)-n_MO:(n_A+lumoB)+n_MO,(n_A+homoB)-n_MO:(n_A+lumoB)+n_MO)
energy((n_MOmax/2+1):n_MOmax,1:(n_MOmax/2))=fock((n_A+homoB)-n_MO:(n_A+lumoB)+n_MO,homoA-n_MO:lumoA+n_MO)
energy(1:(n_MOmax/2),(n_MOmax/2+1):n_MOmax)=transpose(energy((n_MOmax/2+1):n_MOmax,1:(n_MOmax/2)))
Deallocate(fock)
Allocate(fock(n_MOmax,n_MOmax))
fock=energy

!  \mathbf{H}^{\text{eff}}= \mathbf{S}^{- \frac{1}{2}} \mathbf{H}\mathbf{S}^{- \frac{1}{2}}
fock=matmul(transmat,matmul(fock, transmat))

!!! start check energetic degeneration of frontier orbitals in A and B
! Check if the frontier orbital energies are degenerated
IF(allocated(E_Orbitale_A) .AND. allocated(E_Orbitale_B)) THEN
        WRITE(*,*) 'A HOMO , HOMO-1',E_Orbitale_A(homoA),E_Orbitale_A(homoA-1)
        WRITE(*,*) 'A LUMO , LUMO+1',E_Orbitale_A(lumoA),E_Orbitale_A(lumoA+1)
        WRITE(*,*) 'B HOMO , HOMO-1',E_Orbitale_B(homoB),E_Orbitale_B(homoB-1)
        WRITE(*,*) 'B LUMO , LUMO+1',E_Orbitale_B(lumoB),E_Orbitale_B(lumoB+1)
        
        HOMO_A_to_HOMO_B=fock((n_MO+1),(n_MOmax/2+1+n_MO))**2
        S_AB_HOMO_A_to_HOMO_B=S_AB(homoA,homoB)**2
        i_homo=1
        IF(abs(E_Orbitale_A(homoA)-E_Orbitale_A(homoA-1)) < E_deg_cut_off) THEN
                ! A HOMO_A to HOMO_B + HOMO_A -1 to HOMO_B
                HOMO_A_to_HOMO_B=HOMO_A_to_HOMO_B+fock((n_MO),(n_MOmax/2+1+n_MO))**2
                S_AB_HOMO_A_to_HOMO_B=S_AB_HOMO_A_to_HOMO_B+S_AB(homoA-1,homoB)**2
                i_homo=i_homo+1
        END IF
        IF(abs(E_Orbitale_B(homoB)-E_Orbitale_B(homoB-1)) < E_deg_cut_off) THEN
                ! A HOMO_A to HOMO_B + HOMO_A  to HOMO_B -1
                HOMO_A_to_HOMO_B=HOMO_A_to_HOMO_B+fock((n_MO+1),(n_MOmax/2+n_MO))**2
                S_AB_HOMO_A_to_HOMO_B=S_AB_HOMO_A_to_HOMO_B+S_AB(homoA,homoB-1)**2
                i_homo=i_homo+1
        END IF
        IF (i_homo ==3 ) THEN
                ! HOMO_A to HOMO_B + HOMO_A -1 to HOMO_B -1
                HOMO_A_to_HOMO_B=HOMO_A_to_HOMO_B+fock((n_MO),(n_MOmax/2+n_MO))**2
                S_AB_HOMO_A_to_HOMO_B=S_AB_HOMO_A_to_HOMO_B+S_AB(homoA-1,homoB-1)**2
                i_homo=i_homo+1
        END IF
        ! summe J_AB = sqrt( 1/N* sum_{i,j}  |J_ij|*2)  mit i,j in {homo, homo-1}
        HOMO_A_to_HOMO_B=sqrt(1.0d0/(i_homo)*HOMO_A_to_HOMO_B)
        S_AB_HOMO_A_to_HOMO_B=sqrt(S_AB_HOMO_A_to_HOMO_B/Real(i_homo))
        IF (i_homo /= 1) THEN
                WRITE(*,*) 'Caution, degenerated HOMO, HOMO-1 orbital energies found'
                WRITE(*,*) 'Replace J_AB:',fock((n_MO+1),(n_MOmax/2+1+n_MO)),' by ',HOMO_A_to_HOMO_B
                WRITE(*,*) 'and     S_AB:',S_AB(homoA,homoB),' by ',S_AB_HOMO_A_to_HOMO_B
                fock((n_MO+1),(n_MOmax/2+1+n_MO))=HOMO_A_to_HOMO_B
                S_AB(homoA,homoB)=S_AB_HOMO_A_to_HOMO_B
        END IF



        ! check degeneration of lumo, lumo+1 Orbitals
        LUMO_A_to_LUMO_B=fock((n_MO+2),(n_MOmax/2+2+n_MO))**2
        S_AB_LUMO_A_to_LUMO_B=S_AB(lumoA,lumoB)**2
        i_lumo=1
        IF(abs(E_Orbitale_A(lumoA)-E_Orbitale_A(lumoA+1)) < E_deg_cut_off) THEN
                ! A LUMO_A to LUMO_B + LUMO_A +1 to LUMO_B
                LUMO_A_to_LUMO_B=LUMO_A_to_LUMO_B+fock((n_MO+3),(n_MOmax/2+2+n_MO))**2
                S_AB_LUMO_A_to_LUMO_B=S_AB_LUMO_A_to_LUMO_B+S_AB(lumoA+1,lumoB)**2
                i_lumo=i_lumo+1
        END IF
        IF(abs(E_Orbitale_B(lumoB)-E_Orbitale_B(lumoB+1)) < E_deg_cut_off) THEN
                ! A LUMO_A to LUMO_B + LUMO_A  to LUMO_B +1
                LUMO_A_to_LUMO_B=LUMO_A_to_LUMO_B+fock((n_MO+2),(n_MOmax/2+3+n_MO))**2
                S_AB_LUMO_A_to_LUMO_B=S_AB_LUMO_A_to_LUMO_B+S_AB(lumoA,lumoB+1)**2
                i_lumo=i_lumo+1
        END IF
        IF (i_lumo == 3) THEN
                ! LUMO_A to LUMO_B + LUMO_A +1 to LUMO_B +1
                LUMO_A_to_LUMO_B=LUMO_A_to_LUMO_B+fock((n_MO+3),(n_MOmax/2+3+n_MO))**2
                S_AB_LUMO_A_to_LUMO_B=S_AB_LUMO_A_to_LUMO_B+S_AB(lumoA+1,lumoB+1)**2
                i_lumo=i_lumo+1
        END IF
        ! summe J_AB = sqrt( 1/N* sum_{i,j}  |J_ij|*2)  mit i,j in {lumo, lumo+1}
        LUMO_A_to_LUMO_B=sqrt(1.0d0/(i_lumo)*LUMO_A_to_LUMO_B)
        S_AB_LUMO_A_to_LUMO_B=sqrt(S_AB_LUMO_A_to_LUMO_B/Real(i_lumo)) 
        IF (i_lumo /= 1) THEN
                WRITE(*,*) 'Caution, degenerated LUMO, LUMO+1 orbital energies found'
                WRITE(*,*) 'Replace J_AB:',fock((n_MO+2),(n_MOmax/2+2+n_MO)),' by ',LUMO_A_to_LUMO_B
                WRITE(*,*) 'and     S_AB:',S_AB(lumoA,lumoB),' by ',S_AB_LUMO_A_to_LUMO_B
                fock((n_MO+2),(n_MOmax/2+2+n_MO))=LUMO_A_to_LUMO_B
                S_AB(lumoA,lumoB)=S_AB_LUMO_A_to_LUMO_B
        END IF
END IF 
!!! end check energetic degeneration of frontier orbitals in A and B


!printing certain matrix-elements to check validity
		Write(14,*) 'transformed'
		Write(14,FM1) 'HOMO_A to HOMO_A  ', fock((n_MO+1),(n_MO+1))
		Write(14,FM1) 'HOMO_A to LUMO_A  ', fock((n_MO+1),(n_MO+2))
                Write(14,FM1) '1 HOMO_A to HOMO_B', fock((n_MO+1),(n_MOmax/2+1+n_MO))
		Write(14,FM1) 'HOMO_A to LUMO_B  ', fock((n_MO+1),(n_MOmax/2+2+n_MO))
		Write(14,FM1) '2 LUMO_A to LUMO_B', fock((n_MO+2),(n_MOmax/2+2+n_MO))
                Write(14,FM1) 'LUMO_A to LUMO_A  ', fock((n_MO+2),(n_MO+2))
                Write(14,FM1) 'HOMO_B to HOMO_B  ', fock((n_MOmax/2+1+n_MO),(n_MOmax/2+1+n_MO))
                Write(14,FM1) 'LUMO_B to LUMO_B  ', fock((n_MOmax/2+2+n_MO),(n_MOmax/2+2+n_MO))
                
            if( n_MO /= 0) then                
                Write(14,*) '*** mixed transitions ***'
                Write(14,FM1) 'HOMO_A   to HOMO_B-1  ',fock((n_MO+1),(n_MOmax/2+n_MO))
                Write(14,FM1) 'HOMO_A-1 to HOMO_B    ',fock((n_MO),(n_MOmax/2+1+n_MO))
                Write(14,FM1) 'LUMO_A   to LUMO_B+1  ',fock((n_MO+2),(n_MOmax/2+3+n_MO))
                Write(14,FM1) 'LUMO_A+1 to LUMO_B    ',fock((n_MO+3),(n_MOmax/2+2+n_MO))
            end if

            if( n_MO >= 2) then
                Write(14,*) '*** n_MO >= 2 ***'
                Write(14,FM1) 'HOMO_A-2 to HOMO_B    ',fock((n_MO-1),(n_MOmax/2+1+n_MO))
                Write(14,FM1) 'LUMO_A+2 to LUMO_B+1  ',fock((n_MO+4),(n_MOmax/2+3+n_MO))
             end if
!Writing out the effective matrix in the file fock.txt
!Locate the biggest matrix element
!write(*,*) 'Start to print fock matrix '
!fock_max=0.0d0
!i_max=0
!j_max=0
!Do i=1,n_MOmax
	!Do j=1,n_MOmax
		!Write(11,*)i,j, fock(i,j)
		!if(abs(fock(i,j)) > fock_max ) then
			!fock_max=abs(fock(i,j))		
			!i_max=i
			!j_max=j
		!end if
!	End Do
!End Do
 ! Write(14,*) '*** biggest fock-matrix element ***'
 ! Write(14,*) i_max,' to ',j_max,abs(fock(i_max,j_max))

IF(n_MO>=4) THEN
! Print Overlapp of the Monomer-Orbitals
  Write(14,*) '*** overlapp of orbitals S_AB from monomer A and monomer B ***'
  Write(14,FM1) 'S_AB HOMO_A   HOMO_B  ',   S_AB(homoA,homoB)
  Write(14,FM1) 'S_AB HOMO_A-1 HOMO_B-1',   S_AB((homoA-1),(homoB-1))
  Write(14,FM1) 'S_AB HOMO_A-1 HOMO_B  ',   S_AB((homoA-1),(homoB))
  Write(14,FM1) 'S_AB HOMO_A   HOMO_B-1',   S_AB(homoA,(homoB-1))
  Write(14,FM1) 'S_AB LUMO_A   LUMO_B  ',   S_AB(lumoA,lumoB)  
  Write(14,FM1) 'S_AB LUMO_A+1 LUMO_B+1',   S_AB((lumoA+1),(lumoB+1))   
  Write(14,FM1) 'S_AB LUMO_A+1 LUMO_B  ',   S_AB((lumoA+1),(lumoB))  
  Write(14,FM1) 'S_AB LUMO_A   LUMO_B+1',   S_AB((lumoA),(lumoB+1)) 

Write(14,*) '*** overlapp of orbitals |S_AB| holes from monomer A and monomer B and |J ab|***'

Write(14,FM2) '|S_AB| HOMO_A   HOMO_B     ',abs(S_AB(homoA,homoB))        ,' ',abs(fock((n_MO+1),(n_MOmax/2+1+n_MO)) )
Write(14,FM2) '|S_AB| HOMO_A   HOMO_B-1   ',abs(S_AB(homoA,(homoB-1))    ),' ',abs(fock((n_MO+1),(n_MOmax/2+n_MO))   )
Write(14,FM2) '|S_AB| HOMO_A-1 HOMO_B     ',abs(S_AB((homoA-1),(homoB))  ),' ',abs(fock((n_MO),(n_MOmax/2+1+n_MO))   )
Write(14,FM2) '|S_AB| HOMO_A-1 HOMO_B-1   ',abs(S_AB((homoA-1),(homoB-1))),' ',abs(fock((n_MO),(n_MOmax/2+n_MO))     )
         
Write(14,FM2) '|S_AB| HOMO_A-2 HOMO_B     ',abs(S_AB(homoA-2,homoB))      ,' ',abs(fock((n_MO-1),(n_MOmax/2+1+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-2 HOMO_B-1   ',abs(S_AB(homoA-2,homoB-1))    ,' ',abs(fock((n_MO-1),(n_MOmax/2+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-2 HOMO_B-2   ',abs(S_AB(homoA-2,homoB-2))    ,' ',abs(fock((n_MO-1),(n_MOmax/2-1+n_MO)))
         
Write(14,FM2) '|S_AB| HOMO_A-3 HOMO_B     ',abs(S_AB(homoA-3,homoB))      ,' ',abs(fock((n_MO-2),(n_MOmax/2+1+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-3 HOMO_B-1   ',abs(S_AB(homoA-3,homoB-1))    ,' ',abs(fock((n_MO-2),(n_MOmax/2+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-3 HOMO_B-2   ',abs(S_AB(homoA-3,homoB-2))    ,' ',abs(fock((n_MO-2),(n_MOmax/2-1+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-3 HOMO_B-3   ',abs(S_AB(homoA-3,homoB-3))    ,' ',abs(fock((n_MO-2),(n_MOmax/2-2+n_MO)))
         
Write(14,FM2) '|S_AB| HOMO_A-4 HOMO_B     ',abs(S_AB(homoA-4,homoB))      ,' ',abs(fock((n_MO-3),(n_MOmax/2+1+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-4 HOMO_B-1   ',abs(S_AB(homoA-4,homoB-1))    ,' ',abs(fock((n_MO-3),(n_MOmax/2+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-4 HOMO_B-2   ',abs(S_AB(homoA-4,homoB-2))    ,' ',abs(fock((n_MO-3),(n_MOmax/2-1+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-4 HOMO_B-3   ',abs(S_AB(homoA-4,homoB-3))    ,' ',abs(fock((n_MO-3),(n_MOmax/2-2+n_MO)))
Write(14,FM2) '|S_AB| HOMO_A-4 HOMO_B-4   ',abs(S_AB(homoA-4,homoB-4))    ,' ',abs(fock((n_MO-3),(n_MOmax/2-3+n_MO)))
         
Write(14,*) '*** overlapp of orbitals |S_AB| elektrons from monomer A and monomer B and |J ab|***'                                                                     
         
Write(14,FM2) '|S_AB| LUMO_A   LUMO_B     ',abs(S_AB(lumoA,lumoB)        ),' ',abs(fock((n_MO+2),(n_MOmax/2+2+n_MO)) )
Write(14,FM2) '|S_AB| LUMO_A   LUMO_B+1   ',abs(S_AB((lumoA),(lumoB+1))  ),' ',abs(fock((n_MO+2),(n_MOmax/2+3+n_MO)) )

Write(14,FM2) '|S_AB| LUMO_A+1 LUMO_B     ',abs(S_AB((lumoA+1),(lumoB))  ),' ',abs(fock((n_MO+3),(n_MOmax/2+2+n_MO)) )
Write(14,FM2) '|S_AB| LUMO_A+1 LUMO_B+1   ',abs(S_AB((lumoA+1),(lumoB+1))),' ',abs(fock((n_MO+3),(n_MOmax/2+3+n_MO)) )
                                                                      
Write(14,FM2) '|S_AB| LUMO_A+2 LUMO_B     ',abs(S_AB(lumoA+2,lumoB))      ,' ',abs(fock((n_MO+4),(n_MOmax/2+2+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+2 LUMO_B+1   ',abs(S_AB(lumoA+2,lumoB+1))    ,' ',abs(fock((n_MO+4),(n_MOmax/2+3+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+2 LUMO_B+2   ',abs(S_AB(lumoA+2,lumoB+2))    ,' ',abs(fock((n_MO+4),(n_MOmax/2+4+n_MO)))
                                                                     
Write(14,FM2) '|S_AB| LUMO_A+3 LUMO_B     ',abs(S_AB(lumoA+3,lumoB))      ,' ',abs(fock((n_MO+5),(n_MOmax/2+2+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+3 LUMO_B+1   ',abs(S_AB(lumoA+3,lumoB+1))    ,' ',abs(fock((n_MO+5),(n_MOmax/2+3+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+3 LUMO_B+2   ',abs(S_AB(lumoA+3,lumoB+2))    ,' ',abs(fock((n_MO+5),(n_MOmax/2+4+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+3 LUMO_B+3   ',abs(S_AB(lumoA+3,lumoB+3))    ,' ',abs(fock((n_MO+5),(n_MOmax/2+5+n_MO)))
         
Write(14,FM2) '|S_AB| LUMO_A+4 LUMO_B     ',abs(S_AB(lumoA+4,lumoB))      ,' ',abs(fock((n_MO+6),(n_MOmax/2+2+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+4 LUMO_B+1   ',abs(S_AB(lumoA+4,lumoB+1))    ,' ',abs(fock((n_MO+6),(n_MOmax/2+3+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+4 LUMO_B+2   ',abs(S_AB(lumoA+4,lumoB+2))    ,' ',abs(fock((n_MO+6),(n_MOmax/2+4+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+4 LUMO_B+3   ',abs(S_AB(lumoA+4,lumoB+3))    ,' ',abs(fock((n_MO+6),(n_MOmax/2+5+n_MO)))
Write(14,FM2) '|S_AB| LUMO_A+4 LUMO_B+4   ',abs(S_AB(lumoA+4,lumoB+4))    ,' ',abs(fock((n_MO+6),(n_MOmax/2+6+n_MO)))





! Find maximum overlapp in frontier orbitals homo-3 to lumo+3 from monomers
! for hole transport
frontier=10 !homo+/-frontier
i_max=0
j_max=0
S_AB_max=0.0
Do i=(homoA-frontier),(homoA),1
	Do j=(homoB-frontier),(homoB),1		
		if(abs(S_AB(i,j)) > S_AB_max ) then
			S_AB_max=abs(S_AB(i,j))		
			i_max=i
			j_max=j
		end if
	End Do
End Do

  Write(14,*) '*** frontier monomer orbital maximum |S_AB| *** number frontier=',frontier
  Write(14,'(2X,A12,1X,I3,1X,A11,1X,I3,1X,A8,1X,ES14.6)') 'S_AB HOMO_A ',(i_max-homoA),' HOMO_B ',(j_max-homoB),' |S_AB|=',S_AB_max
  Write(14,'(2X,A12,1X,I3,1X,A11,1X,I3,1X,A8,1X,ES14.6)') 'J_AB HOMO_A ',(i_max-homoA),' to HOMO_B ',(j_max-homoB),&
  ' |J_AB|=',fock((n_MO+1+i_max-homoA),(n_MOmax/2+1+n_MO+j_max-homoB))

! for elektron transport
i_max=0
j_max=0
S_AB_max=0.0
Do i=(lumoA),(lumoA+frontier),1
	Do j=(lumoB),(lumoB+frontier),1		
		if(abs(S_AB(i,j)) > S_AB_max ) then
			S_AB_max=abs(S_AB(i,j))		
			i_max=i
			j_max=j
		end if
	End Do
End Do

  Write(14,*) '*** frontier monomer orbital maximum |S_AB| *** '
  Write(14,'(2X,A12,1X,I3,1X,A11,1X,I3,1X,A8,1X,ES14.6)') 'S_AB LUMO_A+',(i_max-lumoA),' LUMO_B+',(j_max-lumoB),' |S_AB|=',S_AB_max
  Write(14,'(2X,A12,1X,I3,1X,A11,1X,I3,1X,A8,1X,ES14.6)') 'J_AB LUMO_A+',(i_max-lumoA),' to LUMO_B+',(j_max-lumoB),&
  ' |J_AB|=',fock((n_MO+2+i_max-lumoA),(n_MOmax/2+2+n_MO+j_max-lumoB))

END IF !!! END  PRINTING for (IF nMO>=4)


! Gewichtets |J_AB|

!write(*,*) 'E_Orbitale_A (homoA)', E_Orbitale_A(homoA), 'homoA=',homoA
!write(*,*) 'E_Orbitale_B (homoB)', E_Orbitale_B(homoB), 'homoB=',homoB

ALLOCATE(w_A(n_MOmax/2))
ALLOCATE(w_B(n_MOmax/2))

T=300            ! K
k_B=8.617332E-5  !eV/K
sumA_lo=0.0
sumB_lo=0.0
sumA_el=0.0
sumB_el=0.0
J_ges_lo_w1=0.0
J_ges_lo_w2=0.0
J_ges_lo_w3=0.0

J_ges_el_w1=0.0
J_ges_el_w2=0.0
J_ges_el_w3=0.0
J_ges_lo_w4=0.0
J_ges_el_w4=0.0

!converting Hartree to eV!
E_Orbitale_A=E_Orbitale_A*27.21138386
E_Orbitale_B=E_Orbitale_B*27.21138386


DO i=1,n_MOmax/4 
    !L\F6cher
    w_A(i)=exp(-abs(E_Orbitale_A(homoA-(i-1))-E_Orbitale_A(homoA))/(k_B*T))
    w_B(i)=exp(-abs(E_Orbitale_B(homoB-(i-1))-E_Orbitale_B(homoB))/(k_B*T))
    sumA_lo=sumA_lo+w_A(i)
    sumB_lo=sumB_lo+w_B(i)
    !Elektronen
    w_A(i+n_MOmax/4)=exp(-abs(E_Orbitale_A(lumoA+(i-1))-E_Orbitale_A(lumoA))/(k_B*T))
    w_B(i+n_MOmax/4)=exp(-abs(E_Orbitale_B(lumoB+(i-1))-E_Orbitale_B(lumoB))/(k_B*T))
    sumA_el=sumA_el+w_A(i+n_MOmax/4)
    sumB_el=sumB_el+w_B(i+n_MOmax/4)
    !write(*,*) 'wA_lo=',w_A(i),'wB_lo=',w_B(i),'wA_el=',w_A(i+n_MOmax/4),'wB_el=',w_B(i+n_MOmax/4)
END DO

DO i=1,n_MOmax/4 !Normierung
    w_A(i)=w_A(i)/sumA_lo
    w_B(i)=w_B(i)/sumB_lo
    w_A(i+n_MOmax/4)=w_A(i+n_MOmax/4)/sumA_el
    w_B(i+n_MOmax/4)=w_B(i+n_MOmax/4)/sumB_el
    !write(*,*) 'NORMIERT: wA_lo=',w_A(i),'wB_lo=',w_B(i),'wA_el=',w_A(i+n_MOmax/4),'wB_el=',w_B(i+n_MOmax/4)
END do

!WRITE(*,*) sumA_el,sumB_el,sumA_lo,sumB_lo
! |J_{1}^{ges}| = \frac{1}{n_{\text{MO}}+1}  \sum_{i=0}^{n_\text{MO}}  \sum_{j=0}^{n_\text{MO}} w_{A} |J_{ij}| 
! Hier mit index-Verschiebung i=>i+1 und n_MOmax/4=n_MO+1
DO i=1,n_MOmax/4
    DO j=1,n_MOmax/4                 !Loecher fock(homoA-i,homoB-j)       
        J_ges_lo_w1=J_ges_lo_w1+w_A(i)*abs(fock((n_MO+1)-(i-1),(n_MOmax/2+1+n_MO)-(j-1)))
        J_ges_lo_w2=J_ges_lo_w2+(w_A(i)*w_B(j))*abs(fock((n_MO+1)-(i-1),(n_MOmax/2+1+n_MO)-(j-1)))
        J_ges_lo_w3=J_ges_lo_w3+w_A(i)*abs(fock((n_MO+1)-(i-1),(n_MOmax/2+1+n_MO)-(j-1)))**2
        J_ges_lo_w4=J_ges_lo_w4+(w_A(i)*w_B(j))*abs(fock((n_MO+1)-(i-1),(n_MOmax/2+1+n_MO)-(j-1)))**2
                                    !Elektronen fock(lumoA+i,lumoB+j)  
        J_ges_el_w1=J_ges_el_w1+w_A(i+n_MOmax/4)*abs(fock((n_MO+2)+(i-1),(n_MOmax/2+2+n_MO)+(j-1)))
        J_ges_el_w2=J_ges_el_w2+(w_A(i+n_MOmax/4)*w_B(j+n_MOmax/4))*abs(fock((n_MO+2)+(i-1),(n_MOmax/2+2+n_MO)+(j-1)) )
        J_ges_el_w3=J_ges_el_w3+w_A(i+n_MOmax/4)*abs(fock((n_MO+2)+(i-1),(n_MOmax/2+2+n_MO)+(j-1)))**2
        J_ges_el_w4=J_ges_el_w4+(w_A(i+n_MOmax/4)*w_B(j+n_MOmax/4))*abs(fock((n_MO+2)+(i-1),(n_MOmax/2+2+n_MO)+(j-1)))**2
     END DO
END DO
J_ges_lo_w1=J_ges_lo_w1/(n_MOmax/4)
J_ges_el_w1=J_ges_el_w1/(n_MOmax/4)

J_ges_lo_w3=J_ges_lo_w3/(n_MOmax/4)
J_ges_el_w3=J_ges_el_w3/(n_MOmax/4)

J_ges_lo_w4=J_ges_lo_w4/(n_MOmax/4)
J_ges_el_w4=J_ges_el_w4/(n_MOmax/4)

write(23,*)'nMO=',Output_nMO(p), ' HOMO_A_to_HOMO_B ', abs(fock((n_MO+1),(n_MOmax/2+1+n_MO)) ),&
&'J_ges_lo_w1=',J_ges_lo_w1,'J_ges_lo_w2=', J_ges_lo_w2,' LUMO_A_to_LUMO_B ',&
&abs(fock((n_MO+2),(n_MOmax/2+2+n_MO)) ),'J_ges_el_w1=',J_ges_el_w1, ' J_ges_el_w2= ',J_ges_el_w2,&
&' J_ges_lo_w3= ',J_ges_lo_w3, ' J_ges_el_w3= ',J_ges_el_w3,&
&' J_ges_lo_w4= ',J_ges_lo_w4, ' J_ges_el_w4= ',J_ges_el_w4




IF( n_MO >= 9 ) THEN
    Write(15,*) '*** overlapp of orbitals |S_AB| elektrons from monomer A and monomer B and |J ab|**2 und w_A(i)*|J ab|**2 ***'                                                                     
!    
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+0 ',abs(S_AB(lumoA,lumoB))    ,' ',abs(fock((n_MO+2),(n_MOmax/2+2+n_MO)))**2 !,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+1 ',abs(S_AB(lumoA,lumoB+1))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+3+n_MO)))**2 !,w_A(2+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+2 ',abs(S_AB(lumoA,lumoB+2))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+4+n_MO)))**2 !,w_A(3+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+3 ',abs(S_AB(lumoA,lumoB+3))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+5+n_MO)))**2 !,w_A(4+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+4 ',abs(S_AB(lumoA,lumoB+4))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+6+n_MO)))**2 !,w_A(5+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+5 ',abs(S_AB(lumoA,lumoB+5))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+7+n_MO)))**2 !,w_A(6+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+6 ',abs(S_AB(lumoA,lumoB+6))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+8+n_MO)))**2 !,w_A(7+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+7 ',abs(S_AB(lumoA,lumoB+7))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+9+n_MO)))**2 !,w_A(8+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+8 ',abs(S_AB(lumoA,lumoB+8))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+10+n_MO)))**2!,w_A(9+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+0 LUMO_B+9 ',abs(S_AB(lumoA,lumoB+9))  ,' ',abs(fock((n_MO+2),(n_MOmax/2+11+n_MO)))**2!,w_A(10+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+11+n_MO)))**2
!                                                                                                 
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+0 ',abs(S_AB(lumoA+1,lumoB))  ,' ',abs(fock((n_MO+3),(n_MOmax/2+2+n_MO)))**2 !,w_A(1+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+1 ',abs(S_AB(lumoA+1,lumoB+1)),' ',abs(fock((n_MO+3),(n_MOmax/2+3+n_MO)))**2 !,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+2 ',abs(S_AB(lumoA+1,lumoB+2)),' ',abs(fock((n_MO+3),(n_MOmax/2+4+n_MO)))**2 !,w_A(3+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+3 ',abs(S_AB(lumoA+1,lumoB+3)),' ',abs(fock((n_MO+3),(n_MOmax/2+5+n_MO)))**2 !,w_A(4+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+4 ',abs(S_AB(lumoA+1,lumoB+4)),' ',abs(fock((n_MO+3),(n_MOmax/2+6+n_MO)))**2 !,w_A(5+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+5 ',abs(S_AB(lumoA+1,lumoB+5)),' ',abs(fock((n_MO+3),(n_MOmax/2+7+n_MO)))**2 !,w_A(6+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+6 ',abs(S_AB(lumoA+1,lumoB+6)),' ',abs(fock((n_MO+3),(n_MOmax/2+8+n_MO)))**2 !,w_A(7+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+7 ',abs(S_AB(lumoA+1,lumoB+7)),' ',abs(fock((n_MO+3),(n_MOmax/2+9+n_MO)))**2 !,w_A(8+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+8 ',abs(S_AB(lumoA+1,lumoB+8)),' ',abs(fock((n_MO+3),(n_MOmax/2+10+n_MO)))**2!,w_A(9+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+1 LUMO_B+9 ',abs(S_AB(lumoA+1,lumoB+9)),' ',abs(fock((n_MO+3),(n_MOmax/2+11+n_MO)))**2!,w_A(10+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+11+n_MO)))**2
!                                                                                                     
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+0 ',abs(S_AB(lumoA+2,lumoB))  ,' ',abs(fock((n_MO+4),(n_MOmax/2+2+n_MO)))**2 !,w_A(1+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+1 ',abs(S_AB(lumoA+2,lumoB+1)),' ',abs(fock((n_MO+4),(n_MOmax/2+3+n_MO)))**2 !,w_A(2+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+2 ',abs(S_AB(lumoA+2,lumoB+2)),' ',abs(fock((n_MO+4),(n_MOmax/2+4+n_MO)))**2 !,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+3 ',abs(S_AB(lumoA+2,lumoB+3)),' ',abs(fock((n_MO+4),(n_MOmax/2+5+n_MO)))**2 !,w_A(4+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+4 ',abs(S_AB(lumoA+2,lumoB+4)),' ',abs(fock((n_MO+4),(n_MOmax/2+6+n_MO)))**2 !,w_A(5+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+5 ',abs(S_AB(lumoA+2,lumoB+5)),' ',abs(fock((n_MO+4),(n_MOmax/2+7+n_MO)))**2 !,w_A(6+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+6 ',abs(S_AB(lumoA+2,lumoB+6)),' ',abs(fock((n_MO+4),(n_MOmax/2+8+n_MO)))**2 !,w_A(7+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+7 ',abs(S_AB(lumoA+2,lumoB+7)),' ',abs(fock((n_MO+4),(n_MOmax/2+9+n_MO)))**2 !,w_A(8+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+8 ',abs(S_AB(lumoA+2,lumoB+8)),' ',abs(fock((n_MO+4),(n_MOmax/2+10+n_MO)))**2!,w_A(9+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+2 LUMO_B+9 ',abs(S_AB(lumoA+2,lumoB+9)),' ',abs(fock((n_MO+4),(n_MOmax/2+11+n_MO)))**2!,w_A(10+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+11+n_MO)))**2
!                                                                                                              
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+0 ',abs(S_AB(lumoA+3,lumoB))  ,' ',abs(fock((n_MO+5),(n_MOmax/2+2+n_MO)))**2 !,w_A(1+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+1 ',abs(S_AB(lumoA+3,lumoB+1)),' ',abs(fock((n_MO+5),(n_MOmax/2+3+n_MO)))**2 !,w_A(2+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+2 ',abs(S_AB(lumoA+3,lumoB+2)),' ',abs(fock((n_MO+5),(n_MOmax/2+4+n_MO)))**2 !,w_A(3+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+3 ',abs(S_AB(lumoA+3,lumoB+3)),' ',abs(fock((n_MO+5),(n_MOmax/2+5+n_MO)))**2 !,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+4 ',abs(S_AB(lumoA+3,lumoB+4)),' ',abs(fock((n_MO+5),(n_MOmax/2+6+n_MO)))**2 !,w_A(5+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+5 ',abs(S_AB(lumoA+3,lumoB+5)),' ',abs(fock((n_MO+5),(n_MOmax/2+7+n_MO)))**2 !,w_A(6+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+6 ',abs(S_AB(lumoA+3,lumoB+6)),' ',abs(fock((n_MO+5),(n_MOmax/2+8+n_MO)))**2 !,w_A(7+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+7 ',abs(S_AB(lumoA+3,lumoB+7)),' ',abs(fock((n_MO+5),(n_MOmax/2+9+n_MO)))**2 !,w_A(8+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+8 ',abs(S_AB(lumoA+3,lumoB+8)),' ',abs(fock((n_MO+5),(n_MOmax/2+10+n_MO)))**2!,w_A(9+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+3 LUMO_B+9 ',abs(S_AB(lumoA+3,lumoB+9)),' ',abs(fock((n_MO+5),(n_MOmax/2+11+n_MO)))**2!,w_A(10+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+11+n_MO)))**2
!                                                                                                  
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+0 ',abs(S_AB(lumoA+4,lumoB))  ,' ',abs(fock((n_MO+6),(n_MOmax/2+2+n_MO)))**2 !,w_A(1+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+1 ',abs(S_AB(lumoA+4,lumoB+1)),' ',abs(fock((n_MO+6),(n_MOmax/2+3+n_MO)))**2 !,w_A(2+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+2 ',abs(S_AB(lumoA+4,lumoB+2)),' ',abs(fock((n_MO+6),(n_MOmax/2+4+n_MO)))**2 !,w_A(3+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+3 ',abs(S_AB(lumoA+4,lumoB+3)),' ',abs(fock((n_MO+6),(n_MOmax/2+5+n_MO)))**2 !,w_A(4+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+4 ',abs(S_AB(lumoA+4,lumoB+4)),' ',abs(fock((n_MO+6),(n_MOmax/2+6+n_MO)))**2 !,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+5 ',abs(S_AB(lumoA+4,lumoB+5)),' ',abs(fock((n_MO+6),(n_MOmax/2+7+n_MO)))**2 !,w_A(6+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+6 ',abs(S_AB(lumoA+4,lumoB+6)),' ',abs(fock((n_MO+6),(n_MOmax/2+8+n_MO)))**2 !,w_A(7+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+7 ',abs(S_AB(lumoA+4,lumoB+7)),' ',abs(fock((n_MO+6),(n_MOmax/2+9+n_MO)))**2 !,w_A(8+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+8 ',abs(S_AB(lumoA+4,lumoB+8)),' ',abs(fock((n_MO+6),(n_MOmax/2+10+n_MO)))**2!,w_A(9+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+4 LUMO_B+9 ',abs(S_AB(lumoA+4,lumoB+9)),' ',abs(fock((n_MO+6),(n_MOmax/2+11+n_MO)))**2!,w_A(10+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+11+n_MO)))**2  
!                                                                                                       
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+0 ',abs(S_AB(lumoA+5,lumoB))  ,' ',abs(fock((n_MO+7),(n_MOmax/2+2+n_MO)))**2 !,w_A(1+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+1 ',abs(S_AB(lumoA+5,lumoB+1)),' ',abs(fock((n_MO+7),(n_MOmax/2+3+n_MO)))**2 !,w_A(2+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+2 ',abs(S_AB(lumoA+5,lumoB+2)),' ',abs(fock((n_MO+7),(n_MOmax/2+4+n_MO)))**2 !,w_A(3+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+3 ',abs(S_AB(lumoA+5,lumoB+3)),' ',abs(fock((n_MO+7),(n_MOmax/2+5+n_MO)))**2 !,w_A(4+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+4 ',abs(S_AB(lumoA+5,lumoB+4)),' ',abs(fock((n_MO+7),(n_MOmax/2+6+n_MO)))**2 !,w_A(5+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+5 ',abs(S_AB(lumoA+5,lumoB+5)),' ',abs(fock((n_MO+7),(n_MOmax/2+7+n_MO)))**2 !,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+6 ',abs(S_AB(lumoA+5,lumoB+6)),' ',abs(fock((n_MO+7),(n_MOmax/2+8+n_MO)))**2 !,w_A(7+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+7 ',abs(S_AB(lumoA+5,lumoB+7)),' ',abs(fock((n_MO+7),(n_MOmax/2+9+n_MO)))**2 !,w_A(8+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+8 ',abs(S_AB(lumoA+5,lumoB+8)),' ',abs(fock((n_MO+7),(n_MOmax/2+10+n_MO)))**2!,w_A(9+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+5 LUMO_B+9 ',abs(S_AB(lumoA+5,lumoB+9)),' ',abs(fock((n_MO+7),(n_MOmax/2+11+n_MO)))**2!,w_A(10+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+11+n_MO)))**2
!                                                                                                      
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+0 ',abs(S_AB(lumoA+6,lumoB))  ,' ',abs(fock((n_MO+8),(n_MOmax/2+2+n_MO)))**2 !,w_A(1+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+1 ',abs(S_AB(lumoA+6,lumoB+1)),' ',abs(fock((n_MO+8),(n_MOmax/2+3+n_MO)))**2 !,w_A(2+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+2 ',abs(S_AB(lumoA+6,lumoB+2)),' ',abs(fock((n_MO+8),(n_MOmax/2+4+n_MO)))**2 !,w_A(3+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+3 ',abs(S_AB(lumoA+6,lumoB+3)),' ',abs(fock((n_MO+8),(n_MOmax/2+5+n_MO)))**2 !,w_A(4+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+4 ',abs(S_AB(lumoA+6,lumoB+4)),' ',abs(fock((n_MO+8),(n_MOmax/2+6+n_MO)))**2 !,w_A(5+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+5 ',abs(S_AB(lumoA+6,lumoB+5)),' ',abs(fock((n_MO+8),(n_MOmax/2+7+n_MO)))**2 !,w_A(6+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+6 ',abs(S_AB(lumoA+6,lumoB+6)),' ',abs(fock((n_MO+8),(n_MOmax/2+8+n_MO)))**2 !,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+7 ',abs(S_AB(lumoA+6,lumoB+7)),' ',abs(fock((n_MO+8),(n_MOmax/2+9+n_MO)))**2 !,w_A(8+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+8 ',abs(S_AB(lumoA+6,lumoB+8)),' ',abs(fock((n_MO+8),(n_MOmax/2+10+n_MO)))**2!,w_A(9+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+6 LUMO_B+9 ',abs(S_AB(lumoA+6,lumoB+9)),' ',abs(fock((n_MO+8),(n_MOmax/2+11+n_MO)))**2!,w_A(10+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+11+n_MO)))**2
!                                                                                                          
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+0 ',abs(S_AB(lumoA+7,lumoB))  ,' ',abs(fock((n_MO+9),(n_MOmax/2+2+n_MO)))**2 !,w_A(1+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+1 ',abs(S_AB(lumoA+7,lumoB+1)),' ',abs(fock((n_MO+9),(n_MOmax/2+3+n_MO)))**2 !,w_A(2+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+2 ',abs(S_AB(lumoA+7,lumoB+2)),' ',abs(fock((n_MO+9),(n_MOmax/2+4+n_MO)))**2 !,w_A(3+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+3 ',abs(S_AB(lumoA+7,lumoB+3)),' ',abs(fock((n_MO+9),(n_MOmax/2+5+n_MO)))**2 !,w_A(4+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+4 ',abs(S_AB(lumoA+7,lumoB+4)),' ',abs(fock((n_MO+9),(n_MOmax/2+6+n_MO)))**2 !,w_A(5+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+5 ',abs(S_AB(lumoA+7,lumoB+5)),' ',abs(fock((n_MO+9),(n_MOmax/2+7+n_MO)))**2 !,w_A(6+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+6 ',abs(S_AB(lumoA+7,lumoB+6)),' ',abs(fock((n_MO+9),(n_MOmax/2+8+n_MO)))**2 !,w_A(7+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+7 ',abs(S_AB(lumoA+7,lumoB+7)),' ',abs(fock((n_MO+9),(n_MOmax/2+9+n_MO)))**2 !,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+8 ',abs(S_AB(lumoA+7,lumoB+8)),' ',abs(fock((n_MO+9),(n_MOmax/2+10+n_MO)))**2!,w_A(9+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+7 LUMO_B+9 ',abs(S_AB(lumoA+7,lumoB+9)),' ',abs(fock((n_MO+9),(n_MOmax/2+11+n_MO)))**2!,w_A(10+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+11+n_MO)))**2
!                                                                                              
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+0 ',abs(S_AB(lumoA+8,lumoB))  ,' ',abs(fock((n_MO+10),(n_MOmax/2+2+n_MO)))**2!,,w_A(1+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+1 ',abs(S_AB(lumoA+8,lumoB+1)),' ',abs(fock((n_MO+10),(n_MOmax/2+3+n_MO)))**2!,,w_A(2+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+2 ',abs(S_AB(lumoA+8,lumoB+2)),' ',abs(fock((n_MO+10),(n_MOmax/2+4+n_MO)))**2!,,w_A(3+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+3 ',abs(S_AB(lumoA+8,lumoB+3)),' ',abs(fock((n_MO+10),(n_MOmax/2+5+n_MO)))**2!,,w_A(4+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+4 ',abs(S_AB(lumoA+8,lumoB+4)),' ',abs(fock((n_MO+10),(n_MOmax/2+6+n_MO)))**2!,,w_A(5+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+5 ',abs(S_AB(lumoA+8,lumoB+5)),' ',abs(fock((n_MO+10),(n_MOmax/2+7+n_MO)))**2!,,w_A(6+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+6 ',abs(S_AB(lumoA+8,lumoB+6)),' ',abs(fock((n_MO+10),(n_MOmax/2+8+n_MO)))**2!,,w_A(7+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+7 ',abs(S_AB(lumoA+8,lumoB+7)),' ',abs(fock((n_MO+10),(n_MOmax/2+9+n_MO)))**2!,,w_A(8+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+8 ',abs(S_AB(lumoA+8,lumoB+8)),' ',abs(fock((n_MO+10),(n_MOmax/2+10+n_MO)))**2!,,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+8 LUMO_B+9 ',abs(S_AB(lumoA+8,lumoB+9)),' ',abs(fock((n_MO+10),(n_MOmax/2+11+n_MO)))**2!,,w_A(10+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+11+n_MO)))**2
!                                                                                                        
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+0 ',abs(S_AB(lumoA+9,lumoB))  ,' ',abs(fock((n_MO+11),(n_MOmax/2+2+n_MO)))**2!,,w_A(1+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+2+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+1 ',abs(S_AB(lumoA+9,lumoB+1)),' ',abs(fock((n_MO+11),(n_MOmax/2+3+n_MO)))**2!,,w_A(2+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+3+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+2 ',abs(S_AB(lumoA+9,lumoB+2)),' ',abs(fock((n_MO+11),(n_MOmax/2+4+n_MO)))**2!,,w_A(3+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+4+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+3 ',abs(S_AB(lumoA+9,lumoB+3)),' ',abs(fock((n_MO+11),(n_MOmax/2+5+n_MO)))**2!,,w_A(4+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+5+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+4 ',abs(S_AB(lumoA+9,lumoB+4)),' ',abs(fock((n_MO+11),(n_MOmax/2+6+n_MO)))**2!,,w_A(5+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+6+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+5 ',abs(S_AB(lumoA+9,lumoB+5)),' ',abs(fock((n_MO+11),(n_MOmax/2+7+n_MO)))**2!,,w_A(6+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+7+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+6 ',abs(S_AB(lumoA+9,lumoB+6)),' ',abs(fock((n_MO+11),(n_MOmax/2+8+n_MO)))**2!,,w_A(7+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+8+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+7 ',abs(S_AB(lumoA+9,lumoB+7)),' ',abs(fock((n_MO+11),(n_MOmax/2+9+n_MO)))**2!,,w_A(8+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+9+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+8 ',abs(S_AB(lumoA+9,lumoB+8)),' ',abs(fock((n_MO+11),(n_MOmax/2+10+n_MO)))**2!,,w_A(9+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+10+n_MO)))**2
!Write(15,FM2)'|S_AB| LUMO_A+9 LUMO_B+9 ',abs(S_AB(lumoA+9,lumoB+9)),' ',abs(fock((n_MO+11),(n_MOmax/2+11+n_MO)))**2!,,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+11+n_MO)))**2
!
!Gewichtet 
IF (DEBUG) THEN
        DO i=1,n_MOmax/4
		IF ((w_A(i+n_MOmax/4) > 1.0E-20) .and. ( w_B(i+n_MOmax/4) > 1.0E-20)) THEN 
                	WRITE(*,*) 'w_A =',w_A(i+n_MOmax/4),'w_B=',w_B(i+n_MOmax/4)
		END IF
        END DO 
END IF
!########
Write(15,*) 'Elektrontransport monomer A to monomer B with overlap |S_AB| transition matrix element |J ab|**2 und w_A(i)*|J ab|**2'
Write(15,FM3)'A+0 B+0 ',abs(S_AB(lumoA,lumoB))    ,abs(fock((n_MO+2),(n_MOmax/2+2+n_MO)))**2 &
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+0 B+1 ',abs(S_AB(lumoA,lumoB+1))  ,abs(fock((n_MO+2),(n_MOmax/2+3+n_MO)))**2 &
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+0 B+2 ',abs(S_AB(lumoA,lumoB+2))  ,abs(fock((n_MO+2),(n_MOmax/2+4+n_MO)))**2 &
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+0 B+3 ',abs(S_AB(lumoA,lumoB+3))  ,abs(fock((n_MO+2),(n_MOmax/2+5+n_MO)))**2 &
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+0 B+4 ',abs(S_AB(lumoA,lumoB+4))  ,abs(fock((n_MO+2),(n_MOmax/2+6+n_MO)))**2 &
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+0 B+5 ',abs(S_AB(lumoA,lumoB+5))  ,abs(fock((n_MO+2),(n_MOmax/2+7+n_MO)))**2 &
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+0 B+6 ',abs(S_AB(lumoA,lumoB+6))  ,abs(fock((n_MO+2),(n_MOmax/2+8+n_MO)))**2 &
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+0 B+7 ',abs(S_AB(lumoA,lumoB+7))  ,abs(fock((n_MO+2),(n_MOmax/2+9+n_MO)))**2 &
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+0 B+8 ',abs(S_AB(lumoA,lumoB+8))  ,abs(fock((n_MO+2),(n_MOmax/2+10+n_MO)))**2&
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+0 B+9 ',abs(S_AB(lumoA,lumoB+9))  ,abs(fock((n_MO+2),(n_MOmax/2+11+n_MO)))**2&
,w_A(1+n_MOmax/4)*abs(fock((n_MO+2),(n_MOmax/2+11+n_MO)))**2
                                                                                            
Write(15,FM3)'A+1 B+0 ',abs(S_AB(lumoA+1,lumoB))  ,abs(fock((n_MO+3),(n_MOmax/2+2+n_MO)))**2 &
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+1 B+1 ',abs(S_AB(lumoA+1,lumoB+1)),abs(fock((n_MO+3),(n_MOmax/2+3+n_MO)))**2 &
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+1 B+2 ',abs(S_AB(lumoA+1,lumoB+2)),abs(fock((n_MO+3),(n_MOmax/2+4+n_MO)))**2 &
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+1 B+3 ',abs(S_AB(lumoA+1,lumoB+3)),abs(fock((n_MO+3),(n_MOmax/2+5+n_MO)))**2 &
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+1 B+4 ',abs(S_AB(lumoA+1,lumoB+4)),abs(fock((n_MO+3),(n_MOmax/2+6+n_MO)))**2 &
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+1 B+5 ',abs(S_AB(lumoA+1,lumoB+5)),abs(fock((n_MO+3),(n_MOmax/2+7+n_MO)))**2 &
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+1 B+6 ',abs(S_AB(lumoA+1,lumoB+6)),abs(fock((n_MO+3),(n_MOmax/2+8+n_MO)))**2 &
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+1 B+7 ',abs(S_AB(lumoA+1,lumoB+7)),abs(fock((n_MO+3),(n_MOmax/2+9+n_MO)))**2 &
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+1 B+8 ',abs(S_AB(lumoA+1,lumoB+8)),abs(fock((n_MO+3),(n_MOmax/2+10+n_MO)))**2&
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+1 B+9 ',abs(S_AB(lumoA+1,lumoB+9)),abs(fock((n_MO+3),(n_MOmax/2+11+n_MO)))**2&
,w_A(2+n_MOmax/4)*abs(fock((n_MO+3),(n_MOmax/2+11+n_MO)))**2
                                                                                    
Write(15,FM3)'A+2 B+0 ',abs(S_AB(lumoA+2,lumoB))  ,abs(fock((n_MO+4),(n_MOmax/2+2+n_MO)))**2 &
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+2 B+1 ',abs(S_AB(lumoA+2,lumoB+1)),abs(fock((n_MO+4),(n_MOmax/2+3+n_MO)))**2 &
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+2 B+2 ',abs(S_AB(lumoA+2,lumoB+2)),abs(fock((n_MO+4),(n_MOmax/2+4+n_MO)))**2 &
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+2 B+3 ',abs(S_AB(lumoA+2,lumoB+3)),abs(fock((n_MO+4),(n_MOmax/2+5+n_MO)))**2 &
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+2 B+4 ',abs(S_AB(lumoA+2,lumoB+4)),abs(fock((n_MO+4),(n_MOmax/2+6+n_MO)))**2 &
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+2 B+5 ',abs(S_AB(lumoA+2,lumoB+5)),abs(fock((n_MO+4),(n_MOmax/2+7+n_MO)))**2 &
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+2 B+6 ',abs(S_AB(lumoA+2,lumoB+6)),abs(fock((n_MO+4),(n_MOmax/2+8+n_MO)))**2 &
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+2 B+7 ',abs(S_AB(lumoA+2,lumoB+7)),abs(fock((n_MO+4),(n_MOmax/2+9+n_MO)))**2 &
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+2 B+8 ',abs(S_AB(lumoA+2,lumoB+8)),abs(fock((n_MO+4),(n_MOmax/2+10+n_MO)))**2&
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+2 B+9 ',abs(S_AB(lumoA+2,lumoB+9)),abs(fock((n_MO+4),(n_MOmax/2+11+n_MO)))**2&
,w_A(3+n_MOmax/4)*abs(fock((n_MO+4),(n_MOmax/2+11+n_MO)))**2
                                                                            
Write(15,FM3)'A+3 B+0 ',abs(S_AB(lumoA+3,lumoB))  ,abs(fock((n_MO+5),(n_MOmax/2+2+n_MO)))**2 &
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+3 B+1 ',abs(S_AB(lumoA+3,lumoB+1)),abs(fock((n_MO+5),(n_MOmax/2+3+n_MO)))**2 &
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+3 B+2 ',abs(S_AB(lumoA+3,lumoB+2)),abs(fock((n_MO+5),(n_MOmax/2+4+n_MO)))**2 &
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+3 B+3 ',abs(S_AB(lumoA+3,lumoB+3)),abs(fock((n_MO+5),(n_MOmax/2+5+n_MO)))**2 &
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+3 B+4 ',abs(S_AB(lumoA+3,lumoB+4)),abs(fock((n_MO+5),(n_MOmax/2+6+n_MO)))**2 &
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+3 B+5 ',abs(S_AB(lumoA+3,lumoB+5)),abs(fock((n_MO+5),(n_MOmax/2+7+n_MO)))**2 &
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+3 B+6 ',abs(S_AB(lumoA+3,lumoB+6)),abs(fock((n_MO+5),(n_MOmax/2+8+n_MO)))**2 &
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+3 B+7 ',abs(S_AB(lumoA+3,lumoB+7)),abs(fock((n_MO+5),(n_MOmax/2+9+n_MO)))**2 &
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+3 B+8 ',abs(S_AB(lumoA+3,lumoB+8)),abs(fock((n_MO+5),(n_MOmax/2+10+n_MO)))**2&
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+3 B+9 ',abs(S_AB(lumoA+3,lumoB+9)),abs(fock((n_MO+5),(n_MOmax/2+11+n_MO)))**2&
,w_A(4+n_MOmax/4)*abs(fock((n_MO+5),(n_MOmax/2+11+n_MO)))**2
                                                                                           
Write(15,FM3)'A+4 B+0 ',abs(S_AB(lumoA+4,lumoB))  ,abs(fock((n_MO+6),(n_MOmax/2+2+n_MO)))**2 &
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+4 B+1 ',abs(S_AB(lumoA+4,lumoB+1)),abs(fock((n_MO+6),(n_MOmax/2+3+n_MO)))**2 &
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+4 B+2 ',abs(S_AB(lumoA+4,lumoB+2)),abs(fock((n_MO+6),(n_MOmax/2+4+n_MO)))**2 &
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+4 B+3 ',abs(S_AB(lumoA+4,lumoB+3)),abs(fock((n_MO+6),(n_MOmax/2+5+n_MO)))**2 &
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+4 B+4 ',abs(S_AB(lumoA+4,lumoB+4)),abs(fock((n_MO+6),(n_MOmax/2+6+n_MO)))**2 &
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+4 B+5 ',abs(S_AB(lumoA+4,lumoB+5)),abs(fock((n_MO+6),(n_MOmax/2+7+n_MO)))**2 &
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+4 B+6 ',abs(S_AB(lumoA+4,lumoB+6)),abs(fock((n_MO+6),(n_MOmax/2+8+n_MO)))**2 &
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+4 B+7 ',abs(S_AB(lumoA+4,lumoB+7)),abs(fock((n_MO+6),(n_MOmax/2+9+n_MO)))**2 &
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+4 B+8 ',abs(S_AB(lumoA+4,lumoB+8)),abs(fock((n_MO+6),(n_MOmax/2+10+n_MO)))**2&
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+4 B+9 ',abs(S_AB(lumoA+4,lumoB+9)),abs(fock((n_MO+6),(n_MOmax/2+11+n_MO)))**2&
,w_A(5+n_MOmax/4)*abs(fock((n_MO+6),(n_MOmax/2+11+n_MO)))**2  
                                                                   
Write(15,FM3)'A+5 B+0 ',abs(S_AB(lumoA+5,lumoB))  ,abs(fock((n_MO+7),(n_MOmax/2+2+n_MO)))**2 &
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+5 B+1 ',abs(S_AB(lumoA+5,lumoB+1)),abs(fock((n_MO+7),(n_MOmax/2+3+n_MO)))**2 &
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+5 B+2 ',abs(S_AB(lumoA+5,lumoB+2)),abs(fock((n_MO+7),(n_MOmax/2+4+n_MO)))**2 &
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+5 B+3 ',abs(S_AB(lumoA+5,lumoB+3)),abs(fock((n_MO+7),(n_MOmax/2+5+n_MO)))**2 &
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+5 B+4 ',abs(S_AB(lumoA+5,lumoB+4)),abs(fock((n_MO+7),(n_MOmax/2+6+n_MO)))**2 &
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+5 B+5 ',abs(S_AB(lumoA+5,lumoB+5)),abs(fock((n_MO+7),(n_MOmax/2+7+n_MO)))**2 &
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+5 B+6 ',abs(S_AB(lumoA+5,lumoB+6)),abs(fock((n_MO+7),(n_MOmax/2+8+n_MO)))**2 &
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+5 B+7 ',abs(S_AB(lumoA+5,lumoB+7)),abs(fock((n_MO+7),(n_MOmax/2+9+n_MO)))**2 &
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+5 B+8 ',abs(S_AB(lumoA+5,lumoB+8)),abs(fock((n_MO+7),(n_MOmax/2+10+n_MO)))**2&
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+5 B+9 ',abs(S_AB(lumoA+5,lumoB+9)),abs(fock((n_MO+7),(n_MOmax/2+11+n_MO)))**2&
,w_A(6+n_MOmax/4)*abs(fock((n_MO+7),(n_MOmax/2+11+n_MO)))**2
                                                                                     
Write(15,FM3)'A+6 B+0 ',abs(S_AB(lumoA+6,lumoB))  ,abs(fock((n_MO+8),(n_MOmax/2+2+n_MO)))**2 &
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+6 B+1 ',abs(S_AB(lumoA+6,lumoB+1)),abs(fock((n_MO+8),(n_MOmax/2+3+n_MO)))**2 &
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+6 B+2 ',abs(S_AB(lumoA+6,lumoB+2)),abs(fock((n_MO+8),(n_MOmax/2+4+n_MO)))**2 &
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+6 B+3 ',abs(S_AB(lumoA+6,lumoB+3)),abs(fock((n_MO+8),(n_MOmax/2+5+n_MO)))**2 &
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+6 B+4 ',abs(S_AB(lumoA+6,lumoB+4)),abs(fock((n_MO+8),(n_MOmax/2+6+n_MO)))**2 &
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+6 B+5 ',abs(S_AB(lumoA+6,lumoB+5)),abs(fock((n_MO+8),(n_MOmax/2+7+n_MO)))**2 &
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+6 B+6 ',abs(S_AB(lumoA+6,lumoB+6)),abs(fock((n_MO+8),(n_MOmax/2+8+n_MO)))**2 &
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+6 B+7 ',abs(S_AB(lumoA+6,lumoB+7)),abs(fock((n_MO+8),(n_MOmax/2+9+n_MO)))**2 &
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+6 B+8 ',abs(S_AB(lumoA+6,lumoB+8)),abs(fock((n_MO+8),(n_MOmax/2+10+n_MO)))**2&
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+6 B+9 ',abs(S_AB(lumoA+6,lumoB+9)),abs(fock((n_MO+8),(n_MOmax/2+11+n_MO)))**2&
,w_A(7+n_MOmax/4)*abs(fock((n_MO+8),(n_MOmax/2+11+n_MO)))**2
                                                                              
Write(15,FM3)'A+7 B+0 ',abs(S_AB(lumoA+7,lumoB))  ,abs(fock((n_MO+9),(n_MOmax/2+2+n_MO)))**2 &
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+7 B+1 ',abs(S_AB(lumoA+7,lumoB+1)),abs(fock((n_MO+9),(n_MOmax/2+3+n_MO)))**2 &
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+7 B+2 ',abs(S_AB(lumoA+7,lumoB+2)),abs(fock((n_MO+9),(n_MOmax/2+4+n_MO)))**2 &
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+7 B+3 ',abs(S_AB(lumoA+7,lumoB+3)),abs(fock((n_MO+9),(n_MOmax/2+5+n_MO)))**2 &
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+7 B+4 ',abs(S_AB(lumoA+7,lumoB+4)),abs(fock((n_MO+9),(n_MOmax/2+6+n_MO)))**2 &
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+7 B+5 ',abs(S_AB(lumoA+7,lumoB+5)),abs(fock((n_MO+9),(n_MOmax/2+7+n_MO)))**2 &
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+7 B+6 ',abs(S_AB(lumoA+7,lumoB+6)),abs(fock((n_MO+9),(n_MOmax/2+8+n_MO)))**2 &
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+7 B+7 ',abs(S_AB(lumoA+7,lumoB+7)),abs(fock((n_MO+9),(n_MOmax/2+9+n_MO)))**2 &
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+7 B+8 ',abs(S_AB(lumoA+7,lumoB+8)),abs(fock((n_MO+9),(n_MOmax/2+10+n_MO)))**2&
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+7 B+9 ',abs(S_AB(lumoA+7,lumoB+9)),abs(fock((n_MO+9),(n_MOmax/2+11+n_MO)))**2&
,w_A(8+n_MOmax/4)*abs(fock((n_MO+9),(n_MOmax/2+11+n_MO)))**2
                                                                                
Write(15,FM3)'A+8 B+0 ',abs(S_AB(lumoA+8,lumoB))  ,abs(fock((n_MO+10),(n_MOmax/2+2+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+8 B+1 ',abs(S_AB(lumoA+8,lumoB+1)),abs(fock((n_MO+10),(n_MOmax/2+3+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+8 B+2 ',abs(S_AB(lumoA+8,lumoB+2)),abs(fock((n_MO+10),(n_MOmax/2+4+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+8 B+3 ',abs(S_AB(lumoA+8,lumoB+3)),abs(fock((n_MO+10),(n_MOmax/2+5+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+8 B+4 ',abs(S_AB(lumoA+8,lumoB+4)),abs(fock((n_MO+10),(n_MOmax/2+6+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+8 B+5 ',abs(S_AB(lumoA+8,lumoB+5)),abs(fock((n_MO+10),(n_MOmax/2+7+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+8 B+6 ',abs(S_AB(lumoA+8,lumoB+6)),abs(fock((n_MO+10),(n_MOmax/2+8+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+8 B+7 ',abs(S_AB(lumoA+8,lumoB+7)),abs(fock((n_MO+10),(n_MOmax/2+9+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+8 B+8 ',abs(S_AB(lumoA+8,lumoB+8)),abs(fock((n_MO+10),(n_MOmax/2+10+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+8 B+9 ',abs(S_AB(lumoA+8,lumoB+9)),abs(fock((n_MO+10),(n_MOmax/2+11+n_MO)))**2&
,w_A(9+n_MOmax/4)*abs(fock((n_MO+10),(n_MOmax/2+11+n_MO)))**2

Write(15,FM3)'A+9 B+0 ',abs(S_AB(lumoA+9,lumoB))  ,abs(fock((n_MO+11),(n_MOmax/2+2+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+2+n_MO)))**2
Write(15,FM3)'A+9 B+1 ',abs(S_AB(lumoA+9,lumoB+1)),abs(fock((n_MO+11),(n_MOmax/2+3+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+3+n_MO)))**2
Write(15,FM3)'A+9 B+2 ',abs(S_AB(lumoA+9,lumoB+2)),abs(fock((n_MO+11),(n_MOmax/2+4+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+4+n_MO)))**2
Write(15,FM3)'A+9 B+3 ',abs(S_AB(lumoA+9,lumoB+3)),abs(fock((n_MO+11),(n_MOmax/2+5+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+5+n_MO)))**2
Write(15,FM3)'A+9 B+4 ',abs(S_AB(lumoA+9,lumoB+4)),abs(fock((n_MO+11),(n_MOmax/2+6+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+6+n_MO)))**2
Write(15,FM3)'A+9 B+5 ',abs(S_AB(lumoA+9,lumoB+5)),abs(fock((n_MO+11),(n_MOmax/2+7+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+7+n_MO)))**2
Write(15,FM3)'A+9 B+6 ',abs(S_AB(lumoA+9,lumoB+6)),abs(fock((n_MO+11),(n_MOmax/2+8+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+8+n_MO)))**2
Write(15,FM3)'A+9 B+7 ',abs(S_AB(lumoA+9,lumoB+7)),abs(fock((n_MO+11),(n_MOmax/2+9+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+9+n_MO)))**2
Write(15,FM3)'A+9 B+8 ',abs(S_AB(lumoA+9,lumoB+8)),abs(fock((n_MO+11),(n_MOmax/2+10+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+10+n_MO)))**2
Write(15,FM3)'A+9 B+9 ',abs(S_AB(lumoA+9,lumoB+9)),abs(fock((n_MO+11),(n_MOmax/2+11+n_MO)))**2&
,w_A(10+n_MOmax/4)*abs(fock((n_MO+11),(n_MOmax/2+11+n_MO)))**2
                  


END IF ! n_MO > 9 Drucken der Daten in Datei 15





DEALLOCATE(w_A)
DEALLOCATE(w_B)
!i_max=0
!i_max=0
!j_max=0
!S_AB_max=0.0
!Do i=(homoA-1),(lumoA+1),1
!	Do j=(homoB-1),(lumoB+1),1		
!		if(abs(S_AB(i,j)) > S_AB_max ) then
!			S_AB_max=abs(S_AB(i,j))		
!			i_max=i
!			j_max=j
!		end if
!	End Do
!End Do





CLOSE(15)
CLOSE(14)
CLOSE(11)
Deallocate(energy)	
Deallocate(overlap)
Deallocate(gamma_A)
Deallocate(gamma_B)
Deallocate(J_AB)
Deallocate(S_AB)
Deallocate(e_A)
Deallocate(e_B)
Deallocate(fock)
Deallocate(transmat)
Deallocate(eigval)
Deallocate(work)
Deallocate(S_diag)

DEAllocate(E_Orbitale_A)
DEAllocate(E_Orbitale_B)

write(*,*) "Durchlauf Nr",p,"ist beendet"
end do

CLOSE(23)

CONTAINS

LOGICAL FUNCTION file_exists(filename)
    implicit NONE
    Character (LEN=*), INTENT(IN) ::filename
    LOGICAL:: Datei_vorhanden
        inquire(file=TRIM(filename),exist=Datei_vorhanden) ! Abfrage ob gro-Datei vorhanden ist.
        IF( .NOT. Datei_vorhanden) THEN 
            WRITE(*,*) ' Error: The_file_does_not_exist. '//TRIM(filename)
            WRITE(*,*) ' STOP'
            STOP
        END IF ! file exists?  
        file_exists=Datei_vorhanden
END FUNCTION ! file_exists

End Program
