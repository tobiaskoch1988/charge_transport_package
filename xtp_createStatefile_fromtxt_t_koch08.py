#!/usr/bin/env python
import sqlite3
import sys
import numpy as np
import os
import shutil
import lxml.etree as lxml
import argparse as ap
import numpy.linalg as lg
import os.path




parser=ap.ArgumentParser(description="Creating state file from a pairs and segment file")
parser.add_argument("-f","--statefile",required=True,type=str, help="statefile to create/insert")
parser.add_argument('-b',"--boxfile",default="",type=str, help="xml file to read boxdata from")
parser.add_argument('-p',"--pairfile",default="",type=str, help="pairfile to parse data from")
parser.add_argument('-s',"--segfile",default="",type=str, help="segment file to parse data from")
parser.add_argument('-c',"--create", action='store_const', const=1,default=0,help="Creates statefile")
parser.add_argument("-t","--type", choices=["h","e","s","t"],type=str,help="Specify which kind of state to import, e,h,s,t")

args=parser.parse_args()


def writeSegtypetoSql(cursor,segtype,segtypeid):
	c.execute('''INSERT INTO segmentTypes(frame,top,id,name,basis,orbfile,torbnrs,coordfile,canRigid) VALUES('0','0','{}','{}','NOT USED','NOT USED','NOT USED','NOT USED',1 ) '''.format(segtypeid,segtype))

def writeMoltoSql(cursor,name,molid):
	c.execute('''INSERT INTO molecules(frame,top,id,name,type) VALUES('0','0','{}','{}','{}') '''.format(molid,name,name))



def writePairtoSql(cursor,pair,statetype):
	has=[0,0,0,0]
	Jeff=[0,0,0,0]
	# k_AB_el '   k_AB_lo	    k_BA_el           k_BA_lo 
	#print 'tmp:',pair[:]
	rate=[0,0,0,0] 
	if args.type=="e":
		has[0]=1
		Jeff[0]=pair[6]
		rate[0]=pair[8]
		rate[1]=pair[10]
	elif args.type=="h":
		has[1]=1
		Jeff[1]=pair[7]
		rate[2]=pair[9]
		rate[3]=pair[11]
	elif args.type=="s":
		has[2]=1
		Jeff[2]=pair[6]
	elif args.type=="t":
		has[3]=1
		Jeff[3]=pair[6]

	c.execute('''INSERT INTO pairs(frame,top,id,seg1,seg2,drX,drY,drZ,lOe,lOh,lOs,lOt,has_e,has_h,has_s,has_t,rate12e,rate21e,rate12h,rate21h,rate12s,rate21s,rate12t,rate21t,Jeff2e,Jeff2h,Jeff2s,Jeff2t,type) VALUES('0','0','{}','{}','{}','{}','{}','{}',0,0,0,0,'{}','{}','{}','{}','{}','{}','{}','{}',0,0,0,0,'{}','{}','{}','{}',0 ) '''.format(int(pair[0]),int(pair[1]),int(pair[2]),pair[3],pair[4],pair[5],has[0],has[1],has[2],has[3],rate[0],rate[1],rate[2],rate[3],Jeff[0],Jeff[1],Jeff[2],Jeff[3]))
	#print int(pair[0]),int(pair[1]),int(pair[2]),pair[3],pair[4],pair[5],has[0],has[1],has[2],has[3],rate[0],rate[1],rate[2],rate[3],Jeff[0],Jeff[1],Jeff[2],Jeff[3]
	

def writeSegtoSql(cursor,seg,statetype,moltype):
	has=[0,0,0,0]
	Eseg=[0,0,0,0]
	if args.type=="e":
		has[0]=1
		Eseg[0]=seg[5]
	elif args.type=="h":
		has[1]=1
		Eseg[1]=seg[5]
	elif args.type=="s":
		has[2]=1
		Eseg[2]=seg[5]
	elif args.type=="t":
		has[3]=1
		Eseg[3]=seg[5]
	c.execute('''INSERT INTO segments (frame,top,id,name,type,mol,posX,posY,posZ,UnCnNe,UnCnNh,UcNcCe,UcNcCh,UcCnNe,UcCnNh,UxXnNs,UxXnNt,UxNxXs,UxNxXt,UxXnNs,UxXnNt,eAnion,eNeutral,eCation,eSinglet,eTriplet,has_e,has_h,has_s,has_t,occPe,occPh,occPs,occPt) VALUES ('0','0','{}','{}','{}','{}','{}','{}','{}','0','0','0','0','0','0','0','0','0','0','0','0','{}','0','{}','{}','{}','{}','{}','{}','{}','0','0','0','0')'''.format(seg[0],seg[1],moltype,seg[0],seg[2],seg[3],seg[4],Eseg[0],Eseg[1],Eseg[2],Eseg[3],has[0],has[1],has[2],has[3],))


if args.pairfile!="" or args.segfile!="":
	if not args.type:
		print "Specify type to read in, -t"
		sys.exit()


types=["e","h","s","t"]
print "Running xtp_createStatefile" 



if os.path.isfile(args.statefile):
	if args.create:
		print"File {} already exists".format(args.statefile)
		sys.exit()
elif not args.create:
	print"File {} does not exist exist, run with -c to create".format(args.statefile)
	sys.exit()
print 'Create new:',args.statefile
con = sqlite3.connect(args.statefile)
c=con.cursor()

if args.create:
	if not os.path.isfile(args.boxfile):
		print("Boxfile is missing")
		sys.exit()
	parser=lxml.XMLParser(remove_comments=True)
	tree = lxml.parse(args.boxfile,parser)
	root = tree.getroot() 
	box=[]

	box.append(float(root.find("boxX").text))
	box.append(float(root.find("boxY").text))
	box.append(float(root.find("boxZ").text))
	if(root.find("boxXY")!=None):
		box.append(float(root.find("boxXY").text))
		box.append(float(root.find("boxXZ").text))
		box.append(float(root.find("boxYZ").text))
	else:
		box.append(0)
		box.append(0)
		box.append(0)


	c.execute('''CREATE TABLE atoms (_id INTEGER PRIMARY KEY AUTOINCREMENT,frame INT NOT NULL,top INT NOT NULL,id INT NOT NULL,name TEXT NOT NULL,type INT NOT NULL,mol INT NOT NULL,seg INT NOT NULL,frag INT NOT NULL,resnr INT NOT NULL,resname TEXT NOT NULL,posX REAL NOT NULL,posY REAL NOT NULL,posZ REAL NOT NULL,weight REAL NOT NULL,element TEXT NOT NULL,qmid INT NOT NULL,qmPosX REAL NOT NULL,qmPosY REAL NOT NULL,qmPosZ REAL NOT NULL)''')
	c.execute('''CREATE TABLE fragments (_id INTEGER PRIMARY KEY AUTOINCREMENT,frame INT NOT NULL,top INT NOT NULL,id INT NOT NULL,name TEXT NOT NULL,type TEXT NOT NULL,mol INT NOT NULL,seg INT NOT NULL,posX REAL NOT NULL,posY REAL NOT NULL,posZ REAL NOT NULL,symmetry INT NOT NULL,leg1 INT NOT NULL,leg2 INT NOT NULL,leg3 INT NOT NULL)''')
	c.execute('''CREATE TABLE frames (_id INTEGER PRIMARY KEY AUTOINCREMENT,id INT NOT NULL,time REAL NOT NULL,step INT NOT NULL,box11 REAL NOT NULL,box12 REAL NOT NULL,box13 REAL NOT NULL,box21 REAL NOT NULL,box22 REAL NOT NULL,box23 REAL NOT NULL,box31 REAL NOT NULL,box32 REAL NOT NULL,box33 REAL NOT NULL,canRigid INT NOT NULL)''')
	c.execute('''CREATE TABLE molecules (_id INTEGER PRIMARY KEY AUTOINCREMENT,frame INT NOT NULL,top INT NOT NULL,id INT NOT NULL,name TEXT NOT NULL,type TEXT NOT NULL)''')
	c.execute('''CREATE TABLE pairs (_id INTEGER PRIMARY KEY AUTOINCREMENT,frame INT NOT NULL,top INT NOT NULL,id INT NOT NULL,seg1 INT NOT NULL,seg2 INT NOT NULL,drX REAL NOT NULL,drY REAL NOT NULL,drZ REAL NOT NULL,lOe REAL DEFAULT 0,lOh REAL DEFAULT 0,lOs REAL DEFAULT 0,lOt REAL DEFAULT 0,has_e INT DEFAULT 0,has_h INT DEFAULT 0,has_s INT DEFAULT 0,has_t INT DEFAULT 0,rate12e REAL DEFAULT 0,rate21e REAL DEFAULT 0,rate12h REAL DEFAULT 0,rate21h REAL DEFAULT 0,rate12s REAL DEFAULT 0,rate21s REAL DEFAULT 0,rate12t REAL DEFAULT 0,rate21t REAL DEFAULT 0,Jeff2e REAL DEFAULT 0,Jeff2h REAL DEFAULT 0,Jeff2s REAL DEFAULT 0,Jeff2t REAL DEFAULT 0,type INT DEFAULT 0)''')
	c.execute('''CREATE TABLE segmentTypes (_id INTEGER PRIMARY KEY AUTOINCREMENT,frame INT NOT NULL,top INT NOT NULL,id INT NOT NULL,name TEXT NOT NULL,basis TEXT NOT NULL,orbfile TEXT NOT NULL,torbnrs TEXT NOT NULL,coordfile TEXT NOT NULL,canRigid INT NOT NULL)''')
	c.execute('''CREATE TABLE segments (_id INTEGER PRIMARY KEY AUTOINCREMENT,frame INT NOT NULL,top INT NOT NULL,id INT NOT NULL,name TEXT NOT NULL,type TEXT NOT NULL,mol INT NOT NULL,posX REAL NOT NULL,posY REAL NOT NULL,posZ REAL NOT NULL,UnCnNe REAL DEFAULT 0,UnCnNh REAL DEFAULT 0,UcNcCe REAL DEFAULT 0,UcNcCh REAL DEFAULT 0,UcCnNe REAL DEFAULT 0,UcCnNh REAL DEFAULT 0,UnXnNs REAL DEFAULT 0,UnXnNt REAL DEFAULT 0,UxNxXs REAL DEFAULT 0,UxNxXt REAL DEFAULT 0,UxXnNs REAL DEFAULT 0,UxXnNt REAL DEFAULT 0,eAnion REAL DEFAULT 0,eNeutral REAL DEFAULT 0,eCation REAL DEFAULT 0,eSinglet REAL DEFAULT 0,eTriplet REAL DEFAULT 0,has_e INT DEFAULT 0,has_h INT DEFAULT 0,has_s INT DEFAULT 0,has_t INT DEFAULT 0,occPe REAL DEFAULT -1,occPh REAL DEFAULT -1,occPs REAL DEFAULT -1,occPt REAL DEFAULT -1)''')
	c.execute('''CREATE TABLE superExchange (_id INTEGER PRIMARY KEY AUTOINCREMENT,frame INT NOT NULL,top INT NOT NULL,type TEXT NOT NULL)''')
	c.execute('''INSERT INTO frames (id,time,step,box11,box12,box13,box21,box22,box23,box31,box32,box33,canRigid) VALUES('0','0','0','{}','{}','{}','{}','{}','{}','{}','{}','{}',1 ) '''.format(box[0],box[3],box[4],box[3],box[1],box[5],box[4],box[5],box[2]))


if args.pairfile!="":
	with open(args.pairfile,"r") as f:
		lines=f.readlines()
		for line in lines:
			if "pair" in line or "#" in line:
				continue
			toc=line.split()
			pair=np.array(toc,dtype=float)
			writePairtoSql(c,pair,args.type)


if args.segfile!="":	
	types=[]
	with open(args.segfile,"r") as f:
		lines=f.readlines()
		for line in lines:
			if "ID" in line or "#" in line:
				continue
			toc=line.split()
			segtype=toc[1]
			if segtype not in types:
				types.append(segtype)
				writeSegtypetoSql(c,segtype,len(types))
			seg=[int(toc[0]),segtype,float(toc[2]),float(toc[3]),float(toc[4]),float(toc[5])]
			writeSegtoSql(c,seg,args.type,len(types))
			writeMoltoSql(c,segtype,int(toc[0]))


con.commit()
con.close()

