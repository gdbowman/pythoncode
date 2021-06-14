import sys

'''
This python script generates a pymol script that colors one chain of a PDB file according to values given in a simple 
2-column text file, where the two numbers are separated by a space. The first column should be integers that correspond 
to the residue numbers for one nucleic acid or protein chain, and the second column gives the values for each residue that 
will be converted to a color (to be mapped onto a structure in pymol). 

An excellent program that calculates a number of standard nucleic acid parameters is X3DNA, which can be found here:  
http://web.x3dna.org (note that for base-pair and base-step parameters, X3DNA by default first lists the residue numbering
of the first nucleic acid chain it encounters. Therefore, changing the order of the chains in the pdb file and rerunning
X3DNA will list parameters based on the first->last residues for the other (now first) chain).

To run this script, the pdb file itself is not needed, just a 2-column text file, with residue numbers in first column and
values (floats) corresponding to each residue number in the second column. Only the first two columns are used in this version, 
and any rows with characters/words (that is, non-float, non-integer) in columns 1 or 2 are ignored.

This script should be run at the command line by typing something like the following:
python3 makepymolscript-color-singlechain-by-2columntext-[version].py parameterfile.txt

The user should edit the parameters below ("EDIT THESE") so that the numbers and variable names are appropriate.
The parameters include variable names used in the pymol script that is generated. 


Notes:
1. This script needs to be run with python3. Running this script should produce a text file that is a coloring script for pymol.

2. Once the output pymol script is created, that script is run within pymol using the "@" symbol to run text files, 
such as "@colorscript-Chd1apo-chainI-ZpParameter.pml".

3. Parameters that the user needs to set are the following:
- molecule name in pymol
    Note that this script assumes that the molecule being colored (with the same molname) has already been loaded up in pymol. Then,
    executing the script created (e.g. "@colorscript-Chd1apo-chainI-SlideParameter.pml") will color that molecule (here called "Chd1apo").
- chain ID - what chain to be colored.
- resiOFFSET_chain1
    Ideally, the first column of the input gives the residue numbering. In that case, the resiOFFSET should be zero. In some
    cases (using X3DNA output), the numbering for the parameters does not necessarily match with that of the PDB file. Also,
    note that Pymol does not appear to work well with negative residue numbering (common for DNA in nucleosome structures), 
    so I generally renumber the residues so that the first residue is 1 anyway. But alternatively, the resiOFFSET_chain1 can be
    set so that the values (second column) are properly mapped to the corresponding residues.
- MIN, MID and MAX values of the parameter of interest - these should be related to the Zp values used for coloring
    For this version, the MAX color is red [1,0,0], MID color is white [1,1,1], and MIN color is blue [0,0,1].
    The program extrapolates colors between MID and MID for shades of red and MID and MIN for shades of blue.
- stopearlyFLAG,stopearlyRESI 
    In some cases, it is useful to only color a segment of the chain. I therefore set up the option to only color up to a particular
    residue (which is stopearlyRESI-1). The program only does this when stopearlyFLAG=1; when stopearlyFLAG is set to zero, all of 
    the residues and parameters will be used.
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
chainID="I"
resiOFFSET_chain1 = 0 ### this is to adjust residue numbering (first column) so that the second column is mapped to the correct residue in PDB file
parametervalueMAX=0.4
parametervalueMIN=-2.0
parametervalueMID=-0.8
colorMAX=[1,0,0] ### red
colorMID=[1,1,1] ### white
colorMIN=[0,0,1] ### blue
colorBLACK=[0,0,0] ### this color is used for the unique cases where the parameter value is far out of range (like -1000),
# which is used as an indication that the nucleotide should be skipped
scriptversion="SlideParameter"

####################################################################################
####################################################################################
#  The code below should not need editing for standard use of this script
####################################################################################
####################################################################################

print("#"*60)
print()
print("Usage: ")
print("python3 thisscript.py DNAparameterfile.txt ")
print()
print("#"*60)

### read in 2-column file given on commandline 
try:
    parameterfileNAME = sys.argv[1]
    file = open(parameterfileNAME,'r')
    parameterfileINPUT=file.readlines()
    file.close()
except IOError:  
    print("\nCannot read ", sys.argv[1])
    print("please supply a file with parameters.\n")
    sys.exit("EXITING...")
except IndexError:
    print("need to supply filename where first column is resi numbering (integer) and second is parameter to color by (float) ")
    sys.exit("EXITING...")

### for each line of the inputfile, write a line of pymol script
pymolscriptname="colorscript-"+molname+"-chain"+chainID+"-"+str(scriptversion)+".pml"
pymoloutput = open(pymolscriptname,"w")

### header for pymol script
pymoloutput.write("### this script was created using the python script {}.".format(sys.argv[0]))
pymoloutput.write("### this script should color DNA residues according to parameters given in file\n")
pymoloutput.write("# {}\n\n".format(parameterfileNAME))
pymoloutput.write("### this script was run with the following parameters: \n")
pymoloutput.write("# molname = {}\n".format(molname))
pymoloutput.write("# chainID = {}\n".format(chainID))
pymoloutput.write("# resiOFFSET_chain1 = {}\n".format(resiOFFSET_chain1))
pymoloutput.write("# parametervalueMAX = {}\n".format(parametervalueMAX))
pymoloutput.write("# parametervalueMIN = {}\n".format(parametervalueMIN))
pymoloutput.write("# parametervalueMID = {}\n".format(parametervalueMID))
pymoloutput.write("# scriptversion = {}\n".format(scriptversion))

### for each line in file, where first and second item are numbers, set the resi to first
#    and use second to generate the color
#    if the second parameter is above MAX, then make color MAX
#    if the second parameter is bewteen MID and MAX, then divide param by MAX and use that fraction
#    to scale all values of colorMAX
#    if the second parameter is below MIN, then make color MIN
#    if the second parameter is between MIN and MID, then divide param by min and use that fraction
#
#    we will use the ranges of MAX-MID and MID-MIN for the division steps above (cannot assume that MID is zero)   

for line in parameterfileINPUT:
    try:  ### only using lines that have two numbers; skipping those with words
        inputDATA = line.split()
        resiCURRENT=int(inputDATA[0]) + resiOFFSET_chain1
        valueRAW=float(inputDATA[1])
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
        pymoloutput.write("set_color tempcolor{}{}{},[{},{},{}]\n".format(molname, chainID,resiCURRENT,tempCOLOR[0],tempCOLOR[1],tempCOLOR[2]))
        pymoloutput.write("color tempcolor{}{}{},{} and chain {} and resi {}\n".format(molname, chainID,resiCURRENT,molname, chainID,resiCURRENT))
    except:
        print("skipping line {}".format(inputDATA))
        continue

pymoloutput.close()

print("#"*60)
print()
print("pymol script written:")
print()
print("\t{}".format(pymolscriptname))
print()
print("#"*60)
print()
