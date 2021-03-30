Program Projection
!Program for calculating transfer integrals between two monomers. Execute 
!this in a folder where ./dimer.log, ./fort.7, ../molA/monomer.log, 
!../molB/monomer.log, ../molA/fort.7, ../molB/fort.7 are available.
!It will print out certain transfer integrals (in eV!) and writes every integral of 
!interest in Ergebnis.txt
!This version does only consider HOMO/LUMO, no degenerations included

implicit none
!overlap-matrix
Real*8, Allocatable:: overlap(:,:), fock(:,:), transmat(:,:)
!expansioncoefficients of monomer A and B and Dimer D
Real*8, Allocatable:: zeta_A(:,:), zeta_B(:,:), zeta_D(:,:)
!dimer-eigenvalues
Real*8, Allocatable:: energy(:,:)
!projectionmatrices 
Real*8, Allocatable:: gamma_A(:,:), gamma_B(:,:)
!transferintegral J_AB
Real*8, Allocatable::J_AB(:,:), e_A(:,:), e_B(:,:), S_AB(:,:)
Integer:: ierror, i, j, l, m, k
Integer:: n_A, n_B, n_D, counter, homoA, homoB, lumoA, lumoB

Character (len=8):: dummy
Character (len=30):: dummy2
Real*8:: dummy3
Character (len=100):: buf


!arrays for transformation-purposes (DSYEV) 
Real*8, Allocatable:: eigval(:), S_diag(:,:), work(:)
Integer:: Info, lwork



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

!searching for the indices of HOMO and LUMO
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

!Output-Format of .log is 5 numbers per row
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

Open(Unit=14, File='Ergebnis.txt', Status='REPLACE', IOSTAT=ierror)

Write(14,*) 'HomoA:', homoA
Write(14,*) 'LumoA:', lumoA
Write(14,*) 'HomoB:', homoB
Write(14,*) 'LumoB:', lumoB
!printing certain matrix-elements to check validity
		Write(14,*) 'untransformed:'
		Write(14,*) 'HOMO_A to HOMO_A', fock(homoA,homoA)
		Write(14,*) 'HOMO_A to LUMO_A', fock(homoA,lumoA)
		Write(14,*) 'HOMO_A to HOMO_B', fock(homoA,homoB+n_A)
		Write(14,*) 'HOMO_A to LUMO_B', fock(homoA,lumoB+n_A)
		Write(14,*) 'LUMO_A to LUMO_B', fock(lumoA, lumoB+n_A)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Deallocate(energy)
Deallocate(zeta_A)
Deallocate(zeta_B)
Deallocate(zeta_D)

!Now we take just the elements we need for HOMO and LUMO
!while taking energy as a dummy-matrix

Allocate(energy(4,4))
energy(1:2,1:2)=transmat(homoA:lumoA,homoA:lumoA)
energy(3:4,3:4)=transmat((n_A+homoB):(n_A+lumoB),(n_A+homoB):(n_A+lumoB))
energy(3:4,1:2)=transmat((n_A+homoB):(n_A+lumoB),homoA:lumoA)
energy(1:2,3:4)=transpose(energy(3:4,1:2))
Deallocate(transmat)
Allocate(transmat(4,4))
transmat=energy

! Now the diagonalization of S and transformation to effective fock-matrix will be done
lwork=3*4
Allocate(eigval(4))
Allocate(work(lwork))
Allocate(S_diag(4,4))

!Initializing the diagonal-matrix and matrix to be diagonalized

CALL DSYEV('V','U', 4 , transmat, 4 , eigval, work, lwork, info)

S_diag=0.0d0
Do i=1,4
	S_diag(i,i)=1.0d0/sqrt(eigval(i))
End Do
Write(*,*) eigval(1)
Write(*,*) eigval(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Having calculated all the eigenvalues and eigenvectors, one can finally build the transformation-matrix


transmat=matmul(transmat,matmul(S_diag,transpose(transmat)))


!finally calculating the effective fock-matrix:
!takiing the relevant elements of fock matrix analogous to transmat
energy(1:2,1:2)=fock(homoA:lumoA,homoA:lumoA)
energy(3:4,3:4)=fock((n_A+homoB):(n_A+lumoB),(n_A+homoB):(n_A+lumoB))
energy(3:4,1:2)=fock((n_A+homoB):(n_A+lumoB),homoA:lumoA)
energy(1:2,3:4)=transpose(energy(3:4,1:2))
Deallocate(fock)
Allocate(fock(4,4))
fock=energy

fock=matmul(transmat,matmul(fock, transmat))

!printing certain matrix-elements to check validity
		Write(14,*) 'transformed'
		Write(14,*) 'HOMO_A to HOMO_A', fock(1,1)
		Write(14,*) 'HOMO_A to LUMO_A', fock(1,2)
		Write(14,*) 'HOMO_A to HOMO_B', fock(1,3)
		Write(14,*) 'HOMO_A to LUMO_B', fock(1,4)
		Write(14,*) 'LUMO_A to LUMO_B', fock(2,4)
                Write(14,*) 'LUMO_A to LUMO_A', fock(2,2)
                Write(14,*) 'HOMO_B to HOMO_B', fock(3,3)
                Write(14,*) 'LUMO_B to LUMO_B', fock(4,4)
End Program
