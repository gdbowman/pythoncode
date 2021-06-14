import sys

'''
This python script generates a pymol script that colors a PDB file according to the Zp parameter values of duplex DNA
calculated by the program X3DNA. To run this script, the pdb file itself is not needed, just the X3DNA parameter file output.

To obtain the X3DNA parameter file, you can go here: http://web.x3dna.org/

This script should be run at the command line by typing something like the following:
python3 makepymolscript-color-DNAchain-by-x3DNAparam-Zp-[version].py X3DNAparameterfile.out

The user should edit the parameters below ("EDIT THESE") so that the numbers and variable names are appropriate.
The parameters include variable names used in the pymol script that is generated. 

Be sure that the X3DNA outputfile has the expected keywords. For the X3DNA version I used, the block 
containing the Zp parameter values starts like this:

****************************************************************************
Classification of each dinucleotide step in a right-handed nucleic acid
structure: A-like; B-like; TA-like, or other cases.

    step       Xp      Yp      Zp     XpH     YpH     ZpH    Form
   1 CC/GG   -3.09    8.85   -0.58   -2.51    8.87   -0.16     B
   2 CG/CG   -3.10    8.86   -0.58   -2.53    8.87   -0.15     B
   3 GC/GC   -2.46    8.92   -0.93   -0.96    8.62   -2.49     B
   4 CC/GG   -2.05    8.27   -0.92   -0.64    8.29   -0.70     B
[snip...]

If the block containing the Zp parameters is different than above, then I would recommend instead using a different script 
that creates the same pymol coloring output script, but using a 2-column CSV file as input. Please see documentation for 
the python script "makepymolscript-color-DNAchainduplex-by-2columnCSV-[version].py" for more info.

Notes:
1. This script needs to be run with python3. Running this script should produce a text file that is a coloring script for pymol.

2. Once the output pymol script is created, that script is run within pymol using the "@" symbol to run text files, 
such as "@colorDNA-Chd1apo-chainIJ-ZpParameter.pml".

3. Parameters that the user needs to set are the following:
- molecule name in pymol
    Note that this script assumes that the molecule being colored (with the same molname) has already been loaded up in pymol. Then,
    executing the script created (e.g. "@colorDNA-Chd1apo-chainIJ-ZpParameter.pml") will color the molecule in pymol (here called 
    "Chd1apo").
- chain IDs - for the DNA chains to be colored.
    There are two expected chains, for double-stranded DNA (or RNA). Only numbering for one chain (the first encountered)
    is given in the first column of the X3DNA output file. Assuming that the two chains are antiparallel, the numbering for
    the first chain is reversed for the second chain, so that the first nucleotide of chain 1 base pairs with the last of chain 2.
    Again, the assumption with this coloring script is that the two chains are antiparallel, so if something else is going on (e.g.
    an RNA hairpin structure), then this coloring script will not give the correct output. 
- resiOFFSET_chain1, resiOFFSET_chain2
    Since here we are looking at base step parameters, it may be necessary to adjust the numbering for the second chain (by adding 1) 
    so that the two base-paired nucleotides have the same color. For most structures, where the two DNA strands are the same length, 
    the default values of 0 and 1 should be fine. These offsets are also useful if the PDB file has a different residue numbering 
    other than 1 for the first residue (by default, X3DNA considers the first nucleotide it encounters to be residue 1, and numbering
    continues from there). So if the first DNA residue is, say 11, then resiOFFSET_chain1 should be set to 10 (since default of 1
    by X3DNA plus 10 will give 11 expected in pymol). Note that Pymol does not appear to work well with negative residue numbering
    (common in nucleosome structures), so I generally renumber the residues so that the first residue is 1 anyway.
- MIN, MID and MAX values of the parameter of interest - these should be related to the Zp values used for coloring
    For this version, the MAX color is red [1,0,0], MID color is white [1,1,1], and MIN color is blue [0,0,1].
    The program extrapolates colors between MID and MID for shades of red and MID and MIN for shades of blue.
- stopearlyFLAG,stopearlyRESI 
    In some cases, it is useful to only color a segment of the DNA. I therefore set up the option to only color DNA up to a particular
    residue (which is stopearlyRESI-1). The program only does this when stopearlyFLAG=1; when stopearlyFLAG is set to zero, all of the 
    DNA residues and parameters will be used.
- scriptversion
    this is a flag that is useful to remind the user what parameter is being used for coloring; alternatively, if trying out coloring
    using different MAX, MID, and MIN values, it is useful to not overwrite previous versions by changing the scriptversion flag. 
    This variable is included in the name of the output filename.

    GDBowman, June 2021
'''
####################################################################################
####################################################################################
### set variables - EDIT THESE
####################################################################################
####################################################################################
molname="Chd1apo"
chainIDfirst="I"
chainIDsecond="J"
resiOFFSET_chain1 = 0 ### this is to adjust numbering so that first column gives resi numbering for  DNA first chain 
resiOFFSET_chain2 = 1 ### this is to adjust numbering so that first column (plus this offset) gives resi numbering for second DNA chain 
parametervalueMAX=2.0
parametervalueMID=0.5
parametervalueMIN=-1.5
colorMAX=[1,0,0] ### red
colorMID=[1,1,1] ### white
colorMIN=[0,0,1] ### blue
colorBLACK=[0,0,0] ### this color is used for the unique cases where the parameter value is far out of range (like -1000),
#which is used as an indication that the nucleotide should be skipped
stopearlyFLAG=0
stopearlyRESI=132
scriptversion="ZpParameter" ### this is useful to remind what parameter coloring by

####################################################################################
####################################################################################
#  The code below should not need editing for standard use of this script
####################################################################################
####################################################################################

print("#"*60)
print()
print("Usage: ")
print("python3 thisscript.py x3dnaoutputfile.out ")
print()
print("#"*60)

### read in X3DNA output file given on commandline 
try:
    parameterfileNAME = sys.argv[1]
    file = open(parameterfileNAME,'r')
    parameterfileINPUT=file.readlines()
    file.close()
    isthiscorrectfiletype = parameterfileINPUT[1].split()
    if isthiscorrectfiletype [0] != "3DNA":
        print("This does not seem to be in the expected X3DNA format")
        sys.exit("Exiting...\n\n")
except IOError:  
    print("\nCannot read ", sys.argv[1])
    print("please give the filename of an X3DNA output file\n")
    sys.exit("EXITING...")
except IndexError:
    print("please give the filename of an X3DNA output file\n")
    sys.exit("EXITING...")

### open the output file that will be created
pymolscriptname="colorDNA-"+molname+"-chain"+chainIDfirst+chainIDsecond+"-v"+str(scriptversion)+".pml"
pymoloutput = open(pymolscriptname,"w")

### header for pymol script
pymoloutput.write("### this script was created using the python script {}.".format(sys.argv[0]))
pymoloutput.write("### this script should color DNA residues according to parameters given in file\n")
pymoloutput.write("# {}\n\n".format(parameterfileNAME))
pymoloutput.write("### this script was run with the following parameters: \n")
pymoloutput.write("# molname = {}\n".format(molname))
pymoloutput.write("# chainIDfirst = {}\n".format(chainIDfirst))
pymoloutput.write("# chainIDsecond = {}\n".format(chainIDsecond))
pymoloutput.write("# resiOFFSET_chain1 = {}\n".format(resiOFFSET_chain1))
pymoloutput.write("# resiOFFSET_chain2 = {}\n".format(resiOFFSET_chain2))
pymoloutput.write("# parametervalueMAX = {}\n".format(parametervalueMAX))
pymoloutput.write("# parametervalueMIN = {}\n".format(parametervalueMIN))
pymoloutput.write("# parametervalueMID = {}\n".format(parametervalueMID))
pymoloutput.write("# scriptversion = {}\n".format(scriptversion))

### First, collect the data block desired from the x3dna output file

founddatablock=0 ### flag to turn on when hit keyword
stopearly=0
parameterZp=[]


for i in range(len(parameterfileINPUT)):
    linetemplist=parameterfileINPUT[i].split()
    #print(i,x3dnaoutput_chd1apo[i])
    if len(linetemplist)==10:
        ### 
        #print(linetemplist)
        if (linetemplist[0]=="Classification") and (linetemplist[1]=="of") and (linetemplist[2]=="each"):
            founddatablock=i
            continue
    elif (founddatablock > 0) and (linetemplist):
        #print(linetemplist)
        if ((linetemplist[0]=="step") or (linetemplist[0]=="structure:")):
            continue ### second line of block, skipping
        if "*****" in parameterfileINPUT[i]:
            break ### end of datablock; escape from for loop
        if "~~~~~~~" in parameterfileINPUT[i]:
            break ### end of datablock; escape from for loop
        if stopearlyFLAG==1 and linetemplist[0]==str(stopearlyRESI):
            break ### only look at the bp on the nucleosome prior to this residue
        ### fill up lists defined above
        parameterZp.append([int(linetemplist[0]),float(linetemplist[4])])

### Here we are going to color two DNA chains at once, which we assume are antiparallel; 
#   assuming that the first column is the residue numbering for first chain, second chain should
#   have the reverse numbering; so we will make two lists, where one is key and other is value in 
#   a dictionary (that maps one residue onto the other)

firstchain = []
for item in parameterZp:
    firstchain.append(item[0])

secondchain=list(firstchain)
secondchain.reverse()

secondchainnumbering = {} ### here, the resi for first chain is key and the value is resi for second chain
for i in range(len(firstchain)):
    secondchainnumbering[firstchain[i]] = secondchain[i] + resiOFFSET_chain2

### for each line in the parameterZp list created above, where first and second item are numbers, set the resi to first
#    and use second to generate the color
#    if the second parameter is above MAX, then make color MAX
#    if the second parameter is bewteen MID and MAX, then divide param by MAX and use that fraction
#    to scale all values of colorMAX
#    if the second parameter is below MIN, then make color MIN
#    if the second parameter is between MIN and MID, then divide param by min and use that fraction
#
#    we will use the ranges of MAX-MID and MID-MIN for the division steps above (cannot assume that MID is zero)

for item in parameterZp:
    resiCURRENT=int(item[0]) + resiOFFSET_chain1
    valueRAW=item[1]
    #print("resi = {}; value = {}".format(resiCURRENT, valueRAW),end="")
    ### here scaling the input parameter (second column) relative to MAX, MID, MIN values
    #   and creating 3-number list for rgb color
    if valueRAW > parametervalueMAX: 
        tempCOLOR = colorMAX
    elif valueRAW < -1000: ### this is the unique case where bpstep is missing, and set parameters to 2000 for reference
        tempCOLOR = colorBLACK
    elif valueRAW < parametervalueMIN:
        tempCOLOR = colorMIN
    elif valueRAW == parametervalueMID:
        tempCOLOR = colorMID
    elif valueRAW > parametervalueMID: # Between MID and MAX
        colorFRACTION = (parametervalueMAX - valueRAW)/(parametervalueMAX - parametervalueMID)
        ### here, only want to change rgb values that are otherwise zero, so making mask of [0,1,1] for red
        tempCOLOR = [colorMAX[0]+(colorMID[0]-colorMAX[0])*colorFRACTION,colorMAX[1]+(colorMID[1]-colorMAX[1])*colorFRACTION,colorMAX[2]+(colorMID[2]-colorMAX[2])*colorFRACTION]
    elif valueRAW < parametervalueMID: #between MID and MIN
        colorFRACTION = (valueRAW - parametervalueMIN)/(parametervalueMID - parametervalueMIN)
        tempCOLOR = [colorMIN[0]+(1-colorMIN[0])*colorFRACTION,colorMIN[1]+(colorMID[1]-colorMIN[1])*colorFRACTION,colorMIN[2]+(colorMID[2]-colorMIN[2])*colorFRACTION]
    else: ### this should never happen
        tempCOLOR = [0,0,0]
    #print("tempCOLOR = {}".format(tempCOLOR))
    ### here need to have a unique color for each residue, chain, mol
    pymoloutput.write("set_color tempcolor{}{}{},[{},{},{}]\n".format(molname, chainIDfirst,resiCURRENT,tempCOLOR[0],tempCOLOR[1],tempCOLOR[2]))
    pymoloutput.write("color tempcolor{}{}{},{} and chain {} and resi {}\n".format(molname, chainIDfirst,resiCURRENT,molname, chainIDfirst,resiCURRENT))
    pymoloutput.write("color tempcolor{}{}{},{} and chain {} and resi {}\n".format(molname, chainIDfirst,resiCURRENT,molname, chainIDsecond,secondchainnumbering[int(item[0])]))

pymoloutput.close()

print("#"*60)
print()
print("pymol script written:")
print()
print("\t{}".format(pymolscriptname))
print()
print("#"*60)
print()
