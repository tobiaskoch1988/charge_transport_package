#!/usr/bin/python
#Programm zum Erzeugen von Gauss-Funktionen
##Module einbinden math und cmath
import math, sys
import numpy as np 
import os
import re
#Startparameter
Inputdateiname="Input"
Outputdateiname=Inputdateiname
Inputdateiname=str(Inputdateiname)+".dat"

#OutputAnfrage =1
Output=1
Zusatz=""   #String fuer Outputdateinamen     
LorentzPeaks=1 # 1 == on: Loretzkurve fuer die Spektrenverbreiterung mit s=HalfWidth  


xMin=0
xMax=800
xStep=8000
FWHM=20  ##0.05 fuer Gauss
s=FWHM/2.0 #  s ist Half-Width fuer LorentzPeaks s=FWHM/2 gibt den dem molden Output !!

#Normierung des Spektrums auf 1 der Maximalintensitaet
NormSpec=1
cm2nm=0                    #Umrechnung cm**-1 nach nm
eV2nm=0
nm2eV=0


if 1==1:
    SCRIPT=sys.argv[0].split('/')
    NARGS=len(sys.argv)-1
    THELP=['-h','-help']

if (NARGS>0):
    if NARGS >1:
        Inputdateiname=sys.argv[1]
        #print sys.argv[1]
        #print sys.argv[1][-4:]
    if ( not ( (sys.argv[1][-4:] == ".dat") or (sys.argv[1][-4:] == ".DAT")) ):
        print " "
        print " wrong input format: *.dat"
        print " "
        sys.exit(0)
    if ( NARGS == 2):
        Outputdateiname=sys.argv[2]
        if (Outputdateiname[-4:] == '.dat'):
            Outputdateiname=Outputdateiname[:-4]
    elif ( NARGS > 3 ): 
        print ' '
        print 'Error: to many arguments '
        print 'Use: '+str(sys.argv[0])+' input.dat output.dat '
        sys.exit(0)

if(LorentzPeaks==1):
    Outputdateiname=str(Outputdateiname)+'_lorentz_s'+str(s)+'nm'

#Einlesen der Daten aus Datei
if (os.path.isfile(Inputdateiname)):
        print 'use:'+str(Inputdateiname)
else:
    print 'Error: File does not exist: '+str(Inputdateiname)
    sys.exit(0)


#Einlesen der Daten aus Datei
with open(Inputdateiname, "r") as f:
    ZeilenINputges = len(f.readlines())-1
    f.close()
    print("Zeileninputanzahl", ZeilenINputges)

SpecDatei = open (Inputdateiname, "r")
print( "Eingelesenes File: ", SpecDatei.name)


Daten = np.zeros((ZeilenINputges,3), float)    # Erzeugung einer N x3 

for i in range(int(ZeilenINputges)):
    str2=SpecDatei.readline()
    str2= re.sub("D","E", str2)
    
    Daten[i,0]=np.float(str2.split()[0])
    Daten[i,1]=np.float(str2.split()[1])
    #Daten[i,2]=np.float(str2.split()[2])
    #print( "Split1", str2.split()[0]) 
    #print( "Split2", str2.split()[1])
    #print( "Split3", str2.split()[2])
      
    
    #Daten[i,:] =list(float(x) for x in str2.split()[0:])  
SpecDatei.close()

# Umrechnung von cm*(-1) in nm

if(cm2nm==1):
    for i in range(int(ZeilenINputges)):
        Daten[i,0]=1.0E+7/Daten[i,0]
elif(eV2nm==1):
    for i in range(int(ZeilenINputges)):
        Daten[i,0]=1239.84192920042/Daten[i,0]
elif(eV2nm==1):
    for i in range(int(ZeilenINputges)):
        Daten[i,0]=1239.84192920042/Daten[i,0]



#Erzeugung eines Array in dem die Gauss-Funktionen gerechnet werden
xArray=np.zeros(int(xStep+1),float)
SpecGes=np.zeros(int(xStep+1),float)
SpecArray=np.zeros((int(ZeilenINputges+1),int(xStep+1)), float)  
sigma=FWHM/(2.0*math.sqrt(2*math.log(2)) )
dx=(xMax-xMin)/(1.0*(xStep))
print("dx=",dx)
print("sigma=",sigma)

for j in range(xStep+1):
    x=xMin+j*dx
    xArray[j]=x
    
for i in range(ZeilenINputges):
    x0=Daten[i,0]
    for j in range(xStep+1):
        x=xArray[j]
        if (LorentzPeaks==1):
            y=Daten[i,1]*s/(s**2+(x-x0)**2)/np.pi
            SpecArray[i,j]=y
            SpecGes[j]=SpecGes[j]+y
            y=0
        else:
            ### Gaussian fit
            y=math.exp(-(x-x0)**2/(2.0*sigma**2) )/( 1.0*sigma*math.sqrt(2*math.pi) )*Daten[i,1]
            SpecArray[i,j]=y
            SpecGes[j]=SpecGes[j]+y
            #SpecArray[int(ZeilenINputges),j]=SpecArray[int(ZeilenINputges+1),j]+y
            y=0


Maxi=max(SpecGes[:])
print( "Maxi = ", Maxi)
for j in range(int(xStep+1)):
    if(SpecGes[j]==Maxi): 
        Jmax=j
        print("Jmax= ",Jmax)
print( "An Maximum an Stelle j=", Jmax, " Wellenlaenge in nm =", xArray[Jmax], "Maximum Intensitaet=", SpecGes[Jmax])      
  
if(NormSpec==1):          #Normierung auf 1
    SpecGes[:]=SpecGes[:]/Maxi
    Zusatz=Zusatz+"-norm"
    print( "Das Spektrum wurde NORMIERT !!!" )



if(Output==1):
        # Ausgabe der Matrix in die Datei
        Dateiname=str(Outputdateiname+Zusatz+".dat")  
        Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
        
      
        for j in range(int(xStep+1)):
            Datout.write(str("  "+"{:12.3f}".format(round(xArray[j],3))+"  "+"{:12.6f}".format(round(SpecGes[j],6))+"\n"))  ##+" "+"{:12.6f}".format(round(Daten[i,2],6))+"\n"))
        Datout.write(str("    "+"\n"))
        Datout.close()
        print("Folgende Datei wurde erzeugt: ",Dateiname)

print("ENDE")
