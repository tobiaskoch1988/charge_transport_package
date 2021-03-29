#!/usr/bin/python2.7
# @author t_koch08
### the idea is to map the neighbours in the definition in A to the neighbours in definition in B and use the smallest distance in A 
### in order to define the distance between the Resid and neighbours in B
### one needs to supply the sorted_neighbours_0deg.ngh file which is created for SCCS-dihedral = 0 deg 
### and resid_cut0_to_resid_r_cut_${cut_angle}.0.dat file which relates the resids in the 0 deg representation to the representation in ${cut_angle}
##
## Variabels:
## ind_A = index for Resid in molekule A
## Resid_B = Resid in molekule B
## Resid_NB_A = Resid for the neighbour (NB) in Molekuel A
## Resid_NB_B = Resid for the neighbour (NB) in Molekuel B
## dist_NB_A = distance to a neighbour (NB) to Resid_A in Molekuel A 
import os,sys,numpy as np, re

if 1==1:																									#check the input
	SCRIPT=sys.argv[0].split('/')
	NARGS=len(sys.argv)-1
	THELP=['-h','-help']

	if NARGS<2 or (NARGS>0 and sys.argv[1] in THELP):
		print""
		print"parameterexpected:",SCRIPT[len(SCRIPT)-1],
		print"[sorted_neighbours_0deg.ngh] [resid_cut0_to_resid_r_cut_30.0.dat] "
		print""
		sys.exit(0)

	if (sys.argv[1][-4:] != ".ngh"):
		print""
		print" wrong input format"
		print"[sorted_neighbours_0deg.ngh] [resid_cut0_to_resid_r_cut_30.0.dat] "
		print""
		sys.exit(0)
	
	if (sys.argv[2][-4:] != ".dat"):
		print""
		print" wrong input format"
		print"[sorted_neighbours_0deg.ngh] [resid_cut0_to_resid_r_cut_30.0.dat] "
		print""
		sys.exit(0)

#	if (sys.argv[3][-4:] != ".gro"):
#		print""
#		print" wrong input format"
#		print"[sorted_neighbours_0deg.ngh] [resid_cut0_to_resid_r_cut.dat] "
#		print""
#		sys.exit(0)


neighbours_in=open(sys.argv[1])
resid_cut0_to_resid_r_cut=open(sys.argv[2])
#gro_file=open(sys.argv[3])
print 'create new_neighbour list with referece data from ',sys.argv[1],' and ',sys.argv[2]
numberofneighbours=12

##### starts counting at 0!
##### contains in the i comp.: [index, distance, distance to 1st jump partner, distance to 2nd ... , usw.]
##### contains in the i comp.: [index, index, index of 1. jump partner, index of second jump partner, ...]
distList = []

###########################
## use selection sort to sort elements in alist, elements are changed in b_list simultaniously 
## make sure alist and b_list have the same length.
def selectionSort(alist,b_list):
	if (len(alist) == len(b_list) ):
   		for fillslot in range(len(alist)-1,0,-1):
			positionOfMax=0
       			for location in range(1,fillslot+1):
           			if alist[location]>alist[positionOfMax]:
               				positionOfMax = location

       			temp = alist[fillslot]
       			alist[fillslot] = alist[positionOfMax]
       			alist[positionOfMax] = temp

       			temp_b = b_list[fillslot]
       			b_list[fillslot] = b_list[positionOfMax]
       			b_list[positionOfMax] = temp_b
	else:
		print 'Error in : selectionSort'
		print 'Lists do not have the same size, EXIT'
		print alist
		print b_list
		sys.exit(1)
##########################

i=0
### read neighbours
while 1:
	line1 = neighbours_in.readline()																#to avoid smaug killing
	if not line1:
		break	
	distList.append(line1.split()[0:])
	##### cast list to float
	if ( i % 2 == 0): #even 
		distList[i] = [float(a) for a in distList[i]]
	else: #odd
		distList[i] = [int(a) for a in distList[i]]
	i = i+1
neighbours_in.close()
### end read neighbours

N_Resids_cut0=int(i/2)
#### contains list with the relation, how residA (rings theta=0) is mapped to residB ( conection of rings)
residA_in_residB=[]
#### new_neighbour_list contains the new resids and distances
new_neighbour_list=[]

### start reading resid_cut0_to_resid_r_cut.dat ; which should carry residA(theta=0) and residB(theta=r_cut) 
i=0
while 1:
	line1 = resid_cut0_to_resid_r_cut.readline()															#to avoid smaug killing
	if not line1:
		break
	line = line1.split()
	if line1[0] != '#':
		residA_in_residB.append([])
		for l in range(len(line)):
			residA_in_residB[i].append(int(line[l]))
		i=i+1

#print residA_in_residB[0][:]
#print residA_in_residB[i-1][:],residA_in_residB[i-1][0],residA_in_residB[i-1][1]
#print new_neighbour_list[i-1]
print 'distList first distances ',distList[0]
print 'distList first resids    ',distList[1]
resid_offset=distList[1][1]
print 'resid_offset:',resid_offset
print "N_Resids_cut0:",N_Resids_cut0
N_Resids_B=int(max(residA_in_residB[:][:])[1])-min(residA_in_residB[:][:])[1]+1
print "N_Resids_B:",N_Resids_B
residA_in_residB[:][1]
resid_offset_B=min(residA_in_residB[:][:])[1]
for i in range(N_Resids_B):
	Resid_B=i+resid_offset_B
	new_neighbour_list.append([])
	new_neighbour_list.append([])
	#new_neighbour_list.append([Resid_B,0])
	#new_neighbour_list.append([Resid_B,Resid_B])

### the idea is to map the neighbours in the definition in A to the neighbours in definition in B and use the smallest distance in A 
### in order to define the distance between the Resid and neighbours in B
## ind_A = index for Resid in molekule A
## Resid_B = Resid in molekule B
## Resid_NB_A = Resid for the neighbour (NB) in Molekuel A
## Resid_NB_B = Resid for the neighbour (NB) in Molekuel B
## dist_NB_A = distance to a neighbour (NB) to Resid_A in Molekuel A 
##N_Resids_cut0=20 #### Loeschen nur zum Test auf 2
for ind_A in range(N_Resids_cut0): 
	j=0
	### get Resid_B in target molecule
	Resid_B=residA_in_residB[(distList[(2*ind_A+1)][0]-resid_offset)][1]
	#print 'Resid_B',Resid_B,'Resid_A',residA_in_residB[(distList[(2*ind_A+1)][0]-resid_offset)][0]

	### run over all neighbours
	for Resid_NB_A in (distList[(2*ind_A+1)][2:]): 
		dist_NB_A=distList[(2*ind_A)][2+j]
		### get the resid of the neighbour(NB) to Resid_B in target molecule
		Resid_NB_B=residA_in_residB[(Resid_NB_A-resid_offset)][1] 
		j=j+1
		
		#print 'indA=',ind_A,' Resid_B=',Resid_B,' dist_NB_A=',dist_NB_A,' Resid_NB_A=',Resid_NB_A,' Resid_NB_B=',Resid_NB_B
		if (Resid_B != Resid_NB_B):
			# test if pairs are already available in new_neighbour_list
			pair_already_in_new_list= False
			### k count index new_neighbours
			k=0
			for New_neighbour in new_neighbour_list[2*(Resid_B-resid_offset)+1][:]:
				#print 'New_neighbour',New_neighbour			
				if New_neighbour == Resid_NB_B:
					pair_already_in_new_list = True
					dist_NB_B=new_neighbour_list[2*(Resid_B-resid_offset)][k]
					### replace if distance is smaller
					if (dist_NB_B > dist_NB_A):
						new_neighbour_list[2*(Resid_B-resid_offset)][k]=dist_NB_A
						#print 'Replace: ',dist_NB_B ,' by ', dist_NB_A
					#print 'k=',k,' Resid_NB_B=',Resid_NB_B,'New_neigh ',new_neighbour_list[2*(Resid_B-resid_offset)+1][k],new_neighbour_list[2*(Resid_B-resid_offset)][k]
				k=k+1

			if ( not pair_already_in_new_list):
				### add to new_neighbour_list 
				new_neighbour_list[2*(Resid_B-resid_offset)].append(dist_NB_A)
				new_neighbour_list[2*(Resid_B-resid_offset)+1].append(Resid_NB_B)
				
			######### check data and add in reverse order, if needed
				rev_pair_already_in_new_list = False
				### k_rev count index rev_new_neighbours
				k_rev=0
				for rev_neighbour in new_neighbour_list[2*(Resid_NB_B-resid_offset)+1][:]:
					#print 'New_neighbour',New_neighbour			
					if rev_neighbour == Resid_B:
						rev_pair_already_in_new_list = True
						rev_dist_NB_B=new_neighbour_list[2*(Resid_NB_B-resid_offset)][k_rev]
						### replace if distance is smaller
						if (rev_dist_NB_B > dist_NB_A):
							new_neighbour_list[2*(Resid_NB_B-resid_offset)][k_rev]=dist_NB_A
					k_rev=k_rev+1
				if ( not rev_pair_already_in_new_list):
					new_neighbour_list[2*(Resid_NB_B-resid_offset)].append(dist_NB_A)
					new_neighbour_list[2*(Resid_NB_B-resid_offset)+1].append(Resid_B)	

			### end check data and add in reverse order, if needed


### sort the list
sorted_list=open("new_sorted_neighbours.ngh",'w')
##with new_neighbour_list as r:
##    for line in sorted(r, key=sortkey_natural):
##        sorted_list.write(line)
##
######## close files
for i in range(N_Resids_B):
	Resid_B=i+resid_offset_B
	list_temp_dist=new_neighbour_list[2*i]
	list_temp_resid=new_neighbour_list[2*i+1]
	selectionSort(list_temp_dist,list_temp_resid)
	sorted_list.write(str(Resid_B)+'	'+str(0)+'	'+str('	'.join(map(str, list_temp_dist[:(numberofneighbours)])))+ '\n')
	sorted_list.write(str(Resid_B)+'	'+str(Resid_B)+'	'+str('	'.join(map(str, list_temp_resid[:(numberofneighbours)])))+ '\n')
	
sorted_list.close()
print 'File created: new_sorted_neighbours.ngh'
print 'Normal termination'
