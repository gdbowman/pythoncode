# pythoncode
Here are some python2 and python3 scripts that I have found useful.

Brief descriptions are given below.
Detailed instructions for each should be in a commented block at the top of each file.

Scripts (python3) to generate a Pymol script that colors a PDB file based on a text file input. 
---
Based on values in the text input, these scripts color in a ramp from blue to white to red.
###
• makepymolscript-color-DNAchainduplex-by-2columnCSV-2021june08a.py

Colors both chains in a DNA duplex according to text input that gives [residue numbering, value] for one chain.

Input expected in CSV format.
###

• makepymolscript-color-DNAchainduplex-by-2columntext-2021june09a.py

Colors both chains in a DNA duplex according to text input that gives [residue numbering, value] for one chain.

Input expected is text file with two columns separated by a space.
###

• makepymolscript-color-DNAchainduplex-by-X3DNAparam-Zp-2021june09a.py

Colors both chains in a DNA duplex according to Zp parameter as calculated by X3DNA program [https://x3dna.org/]

Input expected is X3DNA output file.
###

• makepymolscript-color-singlechain-by-2columnCSV-2021june09a.py

Colors a single chain of a PDB file according to text input that gives [residue numbering, value] in two columns.

Input expected in CSV format.
###

• makepymolscript-color-singlechain-by-2columntext-2021june09a.py
Colors a single chain of a PDB file according to text input that gives [residue numbering, value] in two columns.

Input expected is text file with two columns separated by a space.

###
MISC
---


---
GDBowman
