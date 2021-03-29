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
Integer:: n_A, n_B, n_D, counter, homoA, homoB, lumoA, lumoB, n_MO=0, n1_MO, n2_MO, n_MOmax , frontier

Character (len=8):: dummy
Character (len=30):: dummy2
Real*8:: dummy3, fock_max, S_AB_max
Character (len=100):: buf
Character (len=12) :: FM1='(A26,ES14.6)'         !Format output 1
Character (len=21) :: FM2='(A26,ES14.6,A,ES14.6)' !Format output 2

!arrays for transformation-purposes 
Real*8, Allocatable:: eigval(:), S_diag(:,:), work(:)
Integer:: Info, lwork

!Zusatz fuer viele nMOs
INTEGER, DIMENSION (2) :: Output_nMO=(/0,1/) ! letztes Element wird spaeter aus n_MOmax gesetzt
CHARACTER(30)::filename1


do p=1,size(Output_nMO)
n_MO=Output_nMO(p)
write(*,*) "Berechne Schritt fuer n_MO=",n_MO


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

!Read overlap
Open(Unit=7, File='dimer.log', Status='OLD', IOSTAT=ierror)
do
Read (7,'(A)') buf
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
	
	!looking for the indices of HOMO and LUMO
!	if(energy(i,i)<0) Then
!    	counter=counter+1
!        	else
!        homoA=counter
!        lumoA=counter+1
!    End if
    
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
    
	!looking for the indices of HOMO and LUMO
!	if (energy(i,i)<0.0d0) Then
!    	counter=counter+1
!        	else
!        homoB=counter
!        lumoB=counter+1
!    End if
    
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
CLOSE(7)
CLOSE(8)
CLOSE(9)
CLOSE(10)

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
Open(Unit=14, File='Ergebnis-'//TRIM(filename1)//'MO-V7.txt' , Status='REPLACE', IOSTAT=ierror)

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
     WRITE(*,*) ' STOP program'
     STOP
END IF

Output_nMO(size(Output_nMO))=n_MO !set last element in Output_nMO to the maximum nMO value
n_MO=Output_nMO(p)                ! calculate iteration p
!WRITE(*,*) "Gebe die Anzahl der Orbitale zusaetzlich zum HOMO an von 0 bis ",n_MO
!READ(*,*) n_MO
WRITE(*,*) "Anzahl der HOMO-", n_MO, " und Anzahl der LUMO-",n_MO
n_MOmax=4*(n_MO+1)
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


S_diag=0.0d0
Do i=1,n_MOmax
	S_diag(i,i)=1.0d0/sqrt(eigval(i))
End Do
Write(*,*) eigval(1)
Write(*,*) eigval(2)
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

fock=matmul(transmat,matmul(fock, transmat))



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
Do i=1,n_MOmax
	Do j=1,n_MOmax
		!Write(11,*)i,j, fock(i,j)
		!if(abs(fock(i,j)) > fock_max ) then
			!fock_max=abs(fock(i,j))		
			!i_max=i
			!j_max=j
		!end if
	End Do
End Do
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
frontier=3 !homo+/-frontier
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

  Write(14,*) '*** frontier monomer orbital maximum |S_AB| *** '
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

write(*,*) "Durchlauf Nr",p,"ist beendet"
end do



End Program
