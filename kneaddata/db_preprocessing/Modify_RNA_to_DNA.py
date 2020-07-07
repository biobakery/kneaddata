import sys
import pdb
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import datetime
from Bio.Alphabet import IUPAC




###############################################################################################
#                                                                                             #
#    Modify_RNA_to_DNA.py                                                                     #
#    --------------------                                                                     #
#    Program objective:      Back transcribe RNA sequences to DNA (Replace "U"s with "T"s)    #
#                                                                                             #
#    Usage:                                                                                   #
#    python -u Modify_RNA_to_DNA.py  Input_RNA.fa      Output_DNA.fa                          #
#                                                                                             # 
#    Written by George Weingart george.weingart@gmail.com   07/01/2020                        #
#                                                                                             #
#                                                                                             #
#                                                                                             #    
###############################################################################################


#####################################################################################
#Copyright (C) <2020>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# This file is a component of the KneadData 
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact george.weingart@gmail.com)
#####################################################################################


###############################################################################################
#   Program Start                                                                             #
###############################################################################################
CurrentTime = datetime.datetime.now() 
print("Program Started ", CurrentTime.strftime("%Y-%m-%d %H:%M:%S")  )


cntOut = 0
lSeqs = list()

iDisplayAfter = 1000
InputFasta = sys.argv[1]
Output_Fasta_FileName = sys.argv[2]

print("Input  filename : ", InputFasta)
print("Output filename : ", Output_Fasta_FileName)
print("The program will back transcribe the input fasta sequences into the output fasta sequences")
print("   ")
fasta_sequences = SeqIO.parse(open(InputFasta),'fasta')
cntSeq = 0
for seq in fasta_sequences:
        cntSeq+=1
        if cntSeq % iDisplayAfter ==0:
            CurrentTime = datetime.datetime.now()
            print("Read ", str(cntSeq), " Input Sequences ",CurrentTime.strftime("%Y-%m-%d %H:%M:%S")   )
        SeqID = seq.name
        Seq = (str(seq.seq))
        Updated_Seq = str(seq.seq.back_transcribe())
        seq.seq._data  = Updated_Seq
        lSeqs.append(seq)

CurrentTime = datetime.datetime.now()   
print("   ")
print("Starting to write the output Sequences", CurrentTime.strftime("%Y-%m-%d %H:%M:%S"))
with open(Output_Fasta_FileName, "w") as Output_Fasta:
    for seq in lSeqs:
        cntOut+=1
        if cntOut % iDisplayAfter == 0:
            CurrentTime = datetime.datetime.now()
            print("Wrote ", str(cntOut), " OutputSequences ",CurrentTime.strftime("%Y-%m-%d %H:%M:%S"))
        SeqIO.write(seq, Output_Fasta, "fasta")


print("   ")
print("Total input sequences processed = ", str(cntSeq))
print("Total output sequences written  = ", str(cntOut))
CurrentTime = datetime.datetime.now() 
print("Program completed successfully ",CurrentTime.strftime("%Y-%m-%d %H:%M:%S")  )
sys.exit(0)   
