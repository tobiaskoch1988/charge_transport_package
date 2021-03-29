#!/usr/bin/python3
##Module einbinden math und cmath
import math
import numpy as np 
import os
#Startparameter festlegen!


xyzInputfilename="K9K.xyz"
# Einlesen der Filename.xyz 
FileXYZ = open (xyzInputfilename, "r")

print( "Name of the file: ", FileXYZ.name)
str1=FileXYZ.readline()
N_Elemente=int(str1.split()[0])
print ( "Atome_Anzahl ist:", N_Elemente)
str2=FileXYZ.readline()
print ( "Info:", str2)

Element = [""] * N_Elemente                      # Erzeugt einer leeren Liste mit N_Elementen 
Koordinaten = np.zeros((N_Elemente,3), float)    # Erzeugung einer N x3 Matrix mit Nullen

for i in range(N_Elemente):
    str3=FileXYZ.readline()
    Element[i]=str3.split()[0]
    Koordinaten[i,:]=list(float(x) for x in str3.split()[1:])

#close the FileXYZ
FileXYZ.close() 


#Startparameter festlegen AUFPASSEN!
Nproc=72                    #Anzahl der Prozessoren Palma 12 Smaug 8
Memory=str("64GB")           #Auswahl Memory
Molekuelname=str("K9K")               #Name des Molekueles
Basissatz=str("SDD")                     #Basissatz    6-31G* SDD             
Funktional=str("PBE0_SDD") 
##str("LC-BLYP")  ##str("CAM-B3LYP")                      #Funktional Becke3LYP in G09: PBE1PBE
Funktionalzeile=str("PBE1PBE")  ##"LC-BLYP"  ##str("PBE1PBE")              #str("Becke3LYP")                    
Modelname=str("D3BJ_THF_UFF")         ###str("D3BJ-PCM-CHCL3-Bondi-V1")                        #"PCM-Bondi" #"PCM-Pauling"
Modelzeile=str(" SCF=(XQC,MaxCycle=512,MaxConventionalCycles=300)  EmpiricalDispersion=GD3BJ  SCRF(PCM,Solvent=Tetrahydrofuran,Read)  ")  ###str(" SCRF(PCM,Solvent=Chloroform,Read)  ")
### Modelzeile=str(" Guess=(Huckel)  SCF=(XQC,MaxCycle=512,DIIS,MaxConventionalCycles=300)  EmpiricalDispersion=GD3BJ  SCRF(PCM,Solvent=Tetrahydrofuran,Read)  density=current ")  ###str(" SCRF(PCM,Solvent=Chloroform,Read)  ")

##%"SCRF(PCM,Solvent=Toluene,Read)"

##str("SCRF(PCM,Solvent=Acetonitrile,Read)")

##str("SCRF(PCM,Solvent=Tetrahydrofuran,Read)")  

##"SCRF(PCM,Solvent=Generic,Read)"      #str("SCRF(PCM,Solvent=Tetrahydrofuran,Read)")  #" "        #   #" " Modellzeile
Temperatur=300      
Anhangszeile=str("Radii=UFF")      ###str("Radii=UFF")+"\n"+str("eps=3.2")+"\n"+str("epsinf=2.2216")+"\n"+str("SolventName=PMMA2")+"\n"+str("QConv=Tight")+"\n"+str("GeomView")+"\n"
##str("Radii=UFF")
##str("Radii=Bondi ") 

#" "   "Radii=UFF"    str("Radii=Bondi")  #Zeile unter den xyz Daten Bsp: Radii=Bondi ; Radii=UFF, Radii=Pauling sonst " "
Charge=0   #2        #Ladung fuer normales Molekuel=0
ZustandA=str("S0")
ZustandB=str("T1")                            # S1 oder T1

Info=str("  Geometrie aus:")+str(xyzInputfilename)  # String fuer weitere Informationen ueber die Rechnung

# Zur Auswahl der Berechnungsdateien =1 Setzen
OptA=1                                         # =1 Optimierung Anfordern 
OptAtight=1
OptAverytight=0
OptB=1
OptBtight=1
OptBverytight=0

frA=1
frAtight=1
frAverytight=0
frB=1
frBtight=1
frBverytight=0

EmSpec=1
EmSpectight=1
AbsSpec=0



# Start Output
if ZustandA=="S0":
    MA=1
elif ZustandA=="T1":
    MA=3   
else:
    print( "Es ist Fehler bei der Zustandsbennenung aufgetreten")

S1Modus=str(" ")
if ZustandB=="S1":
   MB=1
   S1Modus=str(" TD(NSTATES=5,Root=1,Singlets) ")     # Zusatz bei Optimierung im S1
elif ZustandB=="T1":
    MB=3
else: 
    print( "Es ist Fehler bei der Zustandsbennenung aufgetreten")    
if( ZustandA==ZustandB):
    print( "Fehler: Anfangs und Endzustand sind gleich :-(") 

if Modelname=="Vac":
    Anhangszeile=" "
    Modelzeile=" " 


if(OptA==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+".chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Opt SCF=(XQC) "+str(Modelzeile)+" NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"   "+str(Info)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str("    "+str(Charge)+"   "+str(MA)+"\n") )
    
    for j in range( len(Koordinaten)):
        Datout.write(str(str(Element[j])+"  "+"{:12.6f}".format(round(Koordinaten[j,0],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,1],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,2],6))+"\n"))
        
            
    Datout.write( str("    "+"\n"))
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n"))
    Datout.close()
    print("Folgende Datei wurde erzeugt: ",Dateiname)


if(OptAtight==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+".chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Opt=Tight Int=UltraFine SCF=(XQC,tight) "+str(Modelzeile)+" NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+"   "+str(Info)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str("    "+str(Charge)+"   "+str(MA)+"\n") )
    
    for j in range( len(Koordinaten)):
        Datout.write(str(str(Element[j])+"  "+"{:12.6f}".format(round(Koordinaten[j,0],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,1],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,2],6))+"\n"))
            
    Datout.write( str("    "+"\n"))
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n"))
    Datout.close()
    print("Folgende Datei wurde erzeugt: ",Dateiname)

if(OptAverytight==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+".chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Opt=VeryTight Int=UltraFine SCF=(XQC,tight,TightLinEq) "+str(Modelzeile)+" NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+"   "+str(Info)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str("    "+str(Charge)+"   "+str(MA)+"\n") )
    
    for j in range( len(Koordinaten)):
        Datout.write(str(str(Element[j])+"  "+"{:12.6f}".format(round(Koordinaten[j,0],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,1],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,2],6))+"\n"))
            
    Datout.write( str("    "+"\n"))
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n"))
    Datout.close()
    print("Folgende Datei wurde erzeugt: ",Dateiname)


if(OptB==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+".chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Opt "+str(S1Modus)+" SCF=(XQC) "+str(Modelzeile)+"   NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"   "+str(Info)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str("    "+str(Charge)+"   "+str(MB)+"\n") )
    
    for j in range( len(Koordinaten)):
        Datout.write(str(str(Element[j])+"  "+"{:12.6f}".format(round(Koordinaten[j,0],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,1],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,2],6))+"\n"))
            
    Datout.write( str("    "+"\n"))
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n"))
    Datout.close()
    print("Folgende Datei wurde erzeugt: ",Dateiname)


if(OptBtight==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+".chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Opt=Tight "+str(S1Modus)+" Int=UltraFine SCF=(XQC,tight) "+str(Modelzeile)+" NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+"   "+str(Info)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str("    "+str(Charge)+"   "+str(MB)+"\n") )
    
    for j in range( len(Koordinaten)):
        Datout.write(str(str(Element[j])+"  "+"{:12.6f}".format(round(Koordinaten[j,0],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,1],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,2],6))+"\n"))
            
    Datout.write( str("    "+"\n"))
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n"))
    Datout.close()
    print("Folgende Datei wurde erzeugt: ",Dateiname)

   
if(OptBverytight==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+".chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Opt=VeryTight "+str(S1Modus)+" Int=UltraFine SCF=(XQC,tight,TightLinEq) "+str(Modelzeile)+" Current= Density NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+"   "+str(Info)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str("    "+str(Charge)+"   "+str(MB)+"\n") )
    
    for j in range( len(Koordinaten)):
        Datout.write(str(str(Element[j])+"  "+"{:12.6f}".format(round(Koordinaten[j,0],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,1],6))+"  "+"{:12.6f}".format(round(Koordinaten[j,2],6))+"\n"))
            
    Datout.write( str("    "+"\n"))
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n"))
    Datout.close()
    print("Folgende Datei wurde erzeugt: ",Dateiname)



   
if(frA==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str("fr"+"-"+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K.chk"+"\n"))
    Datout.write(str("%OldChk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+".chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=SaveNM Temperature="+str(Temperatur)+" "+str(Modelzeile)+" Geom=AllChk Guess=Read NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    if(Anhangszeile==" "):
      Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"   "+str(Info)+"\n"))
      Datout.write( str("    "+"\n") )
    else:
     Datout.write( str(str(Anhangszeile)+"\n"))
     Datout.write( str("    "+"\n") )
     Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"   "+str(Info)+"\n"))
     Datout.write( str("    "+"\n") )
    
    Datout.close() 
    print("Folgende Datei wurde erzeugt: ",Dateiname)
    
# Frequenzberechnungen

if(frAtight==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str("fr"+"-"+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+"-tight"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K-tight.chk"+"\n"))
    Datout.write(str("%OldChk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight.chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=SaveNM Temperature="+str(Temperatur)+" "+str(Modelzeile)+" Geom=AllChk Guess=Read NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    if(Anhangszeile==" "):
      Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+"   "+str(Info)+"\n"))
      Datout.write( str("    "+"\n") )
    else:
     Datout.write( str(str(Anhangszeile)+"\n"))
     Datout.write( str("    "+"\n") )
     Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+"   "+str(Info)+"\n"))
     Datout.write( str("    "+"\n") )

    Datout.close()    
    print("Folgende Datei wurde erzeugt: ",Dateiname)

if(frAverytight==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str("fr"+"-"+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+"-verytight"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K-verytight.chk"+"\n"))
    Datout.write(str("%OldChk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight.chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=SaveNM Temperature="+str(Temperatur)+" "+str(Modelzeile)+" Geom=AllChk Guess=Read Density=Current NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    if(Anhangszeile==" "):
      Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+"   "+str(Info)+"\n"))
      Datout.write( str("    "+"\n") )
    else:
     Datout.write( str(str(Anhangszeile)+"\n"))
     Datout.write( str("    "+"\n") )
     Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+"   "+str(Info)+"\n"))
     Datout.write( str("    "+"\n") )
 
    Datout.close()    
    print("Folgende Datei wurde erzeugt: ",Dateiname)




 
if(frB==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str("fr"+"-"+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K.chk"+"\n"))
    Datout.write(str("%OldChk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+".chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=SaveNM Temperature="+str(Temperatur)+" "+str(Modelzeile)+ " " + str(S1Modus)  + " Geom=AllChk Guess=Read Density=Current  NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    if(Anhangszeile==" "):
      Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"   "+str(Info)+"\n"))
      Datout.write( str("    "+"\n") )
    else:
     Datout.write( str(str(Anhangszeile)+"\n"))
     Datout.write( str("    "+"\n") )
     Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"   "+str(Info)+"\n"))
     Datout.write( str("    "+"\n") )
     Datout.close()   
    print("Folgende Datei wurde erzeugt: ",Dateiname)


if(frBtight==1):
    # Ausgabe der Matrix in die Datei
    Dateiname=str("fr"+"-"+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+"-tight"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K-tight.chk"+"\n"))
    Datout.write(str("%OldChk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight.chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=SaveNM Temperature="+str(Temperatur)+" "+str(Modelzeile)+ " " + str(S1Modus)  + " Geom=AllChk Guess=Read  Density=Current NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    if(Anhangszeile==" "):
      Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+"   "+str(Info)+"\n"))
      Datout.write( str("    "+"\n") )
    else:
     Datout.write( str(str(Anhangszeile)+"\n"))
     Datout.write( str("    "+"\n") )
     Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-tight"+"   "+str(Info)+"\n"))
     Datout.write( str("    "+"\n") )
    Datout.close()    
    print("Folgende Datei wurde erzeugt: ",Dateiname)


if(frBverytight==1):
    # Ausgabe in die Datei
    Dateiname=str("fr"+"-"+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+"-verytight"+".inp")  
    Datout = open (Dateiname, "w+")    # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K-verytight.chk"+"\n"))
    Datout.write(str("%OldChk="+str(Molekuelname)+"-"+"Opt"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight.chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=SaveNM Temperature="+str(Temperatur)+" "+str(Modelzeile)+" "+str(S1Modus)+" Geom=AllChk Guess=Read  Density=Current  NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))
    
    Datout.write( str("    "+"\n") )
    if(Anhangszeile==" "):
      Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+"   "+str(Info)+"\n"))
      Datout.write( str("    "+"\n") )
    else:
     Datout.write( str(str(Anhangszeile)+"\n"))
     Datout.write( str("    "+"\n") )
     Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-verytight"+"   "+str(Info)+"\n"))
     Datout.write( str("    "+"\n") )
     Datout.close()   
    print("Folgende Datei wurde erzeugt: ",Dateiname)

   
if(EmSpec==1):
    # Emissionsspektrum Erzeugung
    Dateiname=str("EmSpec"+"-"+str(Molekuelname)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+".inp")  
    Datout = open (Dateiname, "w+")                                             # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K.chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=(ReadFC,ReadFCHT,FC,SaveNM,Emission) Temperature="+str(Temperatur)+" "+str(Modelzeile)+" Geom=AllChk Guess=Read NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))    
    Datout.write( str("    "+"\n") )
    Datout.write(str("DoTemp"+"\n") )
    Datout.write(str("MaxC1=100"+"\n") )
    Datout.write(str("MaxC2=65"+"\n") )
    Datout.write(str("MaxInt=150"+"\n") )
    Datout.write(str("MaxBands=20"+"\n") )
    Datout.write(str("SpecMin=10000"+"\n") )
    Datout.write(str("SpecMax=30000"+"\n") )
    Datout.write(str("NoRelI00"+"\n") )
    Datout.write(str("ForcePrtSpectrum"+"\n") )
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K.chk"+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.close() 
    print("Folgende Datei wurde erzeugt: ",Dateiname)    

if(EmSpectight==1):
    # Emissionsspektrum Erzeugung
    Dateiname=str("EmSpec"+"-"+str(Molekuelname)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+"-tight.inp")  
    Datout = open (Dateiname, "w+")                                             # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K-tight.chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=(ReadFC,ReadFCHT,FC,SaveNM,Emission) Temperature="+str(Temperatur)+" "+str(Modelzeile)+" Geom=AllChk Guess=Read NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))    
    Datout.write( str("    "+"\n") )
    Datout.write(str("DoTemp"+"\n") )
    Datout.write(str("MaxC1=100"+"\n") )
    Datout.write(str("MaxC2=65"+"\n") )
    Datout.write(str("MaxInt=150"+"\n") )
    Datout.write(str("MaxBands=20"+"\n") )
    Datout.write(str("SpecMin=10000"+"\n") )
    Datout.write(str("SpecMax=30000"+"\n") )
    Datout.write(str("NoRelI00"+"\n") )
    Datout.write(str("ForcePrtSpectrum"+"\n") )
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K-tight.chk"+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.close() 
    print("Folgende Datei wurde erzeugt: ",Dateiname)    

if(AbsSpec==1):
    # Absorptionssionsspektrum Erzeugung--Anforderung durch Weglassen des "EmissionsKeywords"!
    Dateiname=str("AbsSpec"+"-"+str(Molekuelname)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(Temperatur)+"K"+".inp")  
    Datout = open (Dateiname, "w+")                                             # oeffnet beschreibbare Datei, die wenn sie nicht schon existiert neu erzeugt wird.    
    Datout.write(str("$ RunGauss"+"\n"))
    Datout.write(str("%NProcShared="+str(Nproc)+"\n"))
    Datout.write(str("%Mem="+str(Memory)+"\n"))
    Datout.write(str("%Chk="+str(Molekuelname)+"-"+"fr"+"-"+str(ZustandA)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K.chk"+"\n"))
    Datout.write(str("#P "+str(Funktionalzeile)+"/"+str(Basissatz)+" Freq=(ReadFC,ReadFCHT,FC,SaveNM) Temperature="+str(Temperatur)+" "+str(Modelzeile)+" Geom=AllChk Guess=Read NoSymm "+"\n"))
    Datout.write(str("# GFINPUT IOP(6/7=3) 6D 10F"+"\n"))    
    Datout.write( str("    "+"\n") )
    Datout.write(str("DoTemp"+"\n") )
    Datout.write(str("MaxC1=200"+"\n") )
    Datout.write(str("MaxC2=130"+"\n") )
    Datout.write(str("MaxInt=300"+"\n") )
    Datout.write(str("MaxBands=40"+"\n") )
    Datout.write(str("SpecMin=2000"+"\n") )
    Datout.write(str("SpecMax=35000"+"\n") )
    Datout.write(str("NoRelI00"+"\n") )
    #Datout.write(str("ForcePrtSpectrum"+"\n") )
    Datout.write( str("    "+"\n") )
    Datout.write(str(str(Molekuelname)+"-"+"fr"+"-"+str(ZustandB)+"-"+str(Funktional)+"-"+str(Modelname)+"-"+str(int(Temperatur))+"K.chk"+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.write( str(str(Anhangszeile)+"\n"))
    Datout.write( str("    "+"\n") )
    Datout.close() 
    print("Folgende Datei wurde erzeugt: ",Dateiname)    

   
print( "Fertig:","MA=",MA, "MB=",MB )  
