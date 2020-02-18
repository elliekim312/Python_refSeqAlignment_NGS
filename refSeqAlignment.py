#!/cm/shared/Tools/python/anaconda/bin/python

## RefFlat Mission
## Coding by Eunjung Kim
## Date: 2017. 08. 28

from __future__ import print_function
from sys import argv
from scipy import stats
from decimal import*
import os, operator


## Open hg19 data
cwd = os.getcwd()       #start path from the current location
os.chdir('../Mission1/'+'hg19') #find hg19 directory

chrDic = {}

for i in list(range(1,23)) + ['X','Y']:
        print(i)
        sSeq    = open('chr%s' %i + '.fa','r')
        row     = sSeq.readlines()
        sequence = ""   #initialize the sequence
        for sNuc in row:
                if sNuc.startswith('>'):
                        continue
                else:
                        sequence += sNuc.rstrip('\n').upper()
        chrDic[str(i)] = sequence       #store sequence
os.chdir(cwd)   #back to the

## Open refFlat.txt through the below path
cwd = os.getcwd()
os.chdir('../../data/')
InFile  = open("refFlat.txt", "r")


## Open Mission5_Dataset1.txt
Dataset = open('Mission5_Dataset3.txt', 'r')
lData = []
#print self.sGeneSym    #test
for line in Dataset.readlines():
        Dataline = line.strip().split('\t')
        lData.append(Dataline)

def main():
        Aligned_NM(cRefSeqList) #Mission1 - Mission3
        Call_RegulMotif(MotifList) ##Mission4
        Call_miRNA(mData) ##Mission5
        Sort_Enrichment(lUniqueMotif) ##Mission5

class cRefSeq:
        ## Initialize RefFlat variables
        def __init__(self):
                self.data       = []            #Store all data into list
                self.sGeneSym   = 'NULL'
                self.sNMID      = 'NULL'
                self.sChrID     = 'NULL'
                self.sChrStrand = 'NULL'
                self.ntxStart   = 0
                self.ntxEnd     = 0
                self.ncodStart  = 0
                self.ncodEnd    = 0
                self.nNumExon   = 0
                self.ExonStart  = []
                self.ExonEnd    = []

                self.ExonSeq    = 'NULL'
                self.UTR5       = 'NULL'
                self.UTR3       = 'NULL'
                self.ORF        = 'NULL'

                self.lRegulation= []
                self.Regulated  = 'NULL'
                self.fRegulationValue = 0.0

        ##Parse the information of RefFlat.txt
        def Parse_RefFlat(self, sReadLine):
                self.data       = sReadLine.strip().split('\t') #Readline by tab
                self.sGeneSym   = self.data[0].upper()
                self.sNMID      = self.data[1]
                self.sChrID     = self.data[2]
                self.sChrStrand = self.data[3]
                self.ntxStart   = int(self.data[4])
                self.ntxEnd     = int(self.data[5])
                self.ncodStart  = int(self.data[6])
                self.ncodEnd    = int(self.data[7])
                self.nNumExon   = int(self.data[8])
                self.ExonStart  = map(int, self.data[9].strip(',').split(','))
                self.ExonEnd    = map(int, self.data[10].strip(',').split(','))


        def GetCodon(self, ExonStart, ExonEnd, ncodStart, ncodEnd, strand):
                cds             = 0     #codon start location of exon
                cde             = 0     #codon end location of exon
                Exon5Start      = []
                Exon5End        = []
                ORFStart        = []
                ORFEnd          = []
                Exon3Start      = []
                Exon3End        = []


                ## Get location of start/end codon
                for i in range(len(ExonStart)):
                        if ExonStart[i] <= ncodStart <= ExonEnd[i]:
                                cds = i
                        if ExonStart[i] <= ncodEnd <= ExonEnd[i]:
                                cde = i

                ## Get 5'UTR coordinates
                for i in range(cds+1):
                        Exon5Start.append(ExonStart[i])
                        Exon5End.append(ExonEnd[i])
                Exon5End[-1]    = ncodStart     #change the last position to startcodon


                ## Get ORF coordinates
                for i in range(cds,cde+1):
                        ORFStart.append(ExonStart[i])
                        ORFEnd.append(ExonEnd[i])
                ORFStart[0] = ncodStart
                ORFEnd[-1] = ncodEnd


                ## Get 3'UTR c  rdinates
                for i in range(cde,self.nNumExon):
                        Exon3Start.append(ExonStart[i])
                        Exon3End.append(ExonEnd[i])
                Exon3Start[0] = ncodEnd

                sequence = chrDic[self.sChrID[3:]]


                self.ExonSeq    = Get_Seq(ExonStart, ExonEnd, sequence)                 #Get the Exon Sequence
                self.UTR5       = Get_Seq(Exon5Start, Exon5End, sequence)               #UTR5 Seq
                self.ORF        = Get_Seq(ORFStart, ORFEnd, sequence)                   #ORF Seq
                self.UTR3       = Get_Seq(Exon3Start, Exon3End, sequence)               #UTR3 Seq

                if strand == '-':
                        s3UTR = self.UTR5
                        s5UTR = self.UTR3
                        self.ExonSeq    = Reverse_Complement(self.ExonSeq)               #Get the Exon Sequence
                        self.UTR5       = Reverse_Complement(s5UTR)                      #UTR5 Seq
                        self.ORF        = Reverse_Complement(self.ORF)                   #ORF Seq
                        self.UTR3       = Reverse_Complement(s3UTR)


                ## Find Stop Codon
                self.stopcodon  = 0
                for i in range(0, len(self.ORF)-3, 3):
                        if self.ORF[i:i+3] in ["TAG", "TGA", "TAA"]:
                                self.stopcodon += 1

        ## throuhg Mission4_Dataset, find the reguation
        def Get_Regulation(self, sGeneSym):
                for i in lData:
                        if sGeneSym == i[0]:
                                self.fRegulationValue = float(i[1])
                                if Decimal(str(self.fRegulationValue)) < Decimal(str(-0.50)):
                                        self.Regulated = 'HighlyDown Regulated'
                                else:
                                        self.Regulated = 'NotHighlyDown Regulated'

def Get_Seq(lExonStart, lExonEnd, sequence):
        sConcatSeq      = ""
        for i in range(0, len(lExonStart)):
                sConcatSeq += sequence[lExonStart[i]:lExonEnd[i]]
        return sConcatSeq


def Reverse_Complement(dna):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U':'A'}
        return ''.join([complement[base] for base in reversed(dna)])

##################################################################################### Mission1 and 2
cRefSeqList = []                                #Store the refFlat.txt into a list of classes
for sReadLine in InFile.readlines():            #Open refFlat.txt and read by line
        cReadRefSeq = cRefSeq()                 #Creating a new instance object named cReadRefSeq
        cReadRefSeq.Parse_RefFlat(sReadLine)    #Calling function "Parse_RefFlat()"
        cRefSeqList.append(cReadRefSeq)         #Store them into a list of class
#       cReadRefSeq.GetRegualtion(sReadLine)
print ('Answer1 : ', str(len(cRefSeqList)))     #Number of RefFlat List


cRefNM = []                                     #List for aligned NM sequence
cRefMultiNM = []        #NM  multipul seq list
cRefNewNM = []          #NM  w/o multipul list
lcRef = []
cRefORF = []            #ORF seq list
cUnique = []
cRefUnique = []         #unique sGeneSym
Isoform = []


## Removed non-NM and Keep the ones aligned only on chr 1-22, X, and Y
def Aligned_NM(cRefSeqList):
        for cRefSeq in cRefSeqList:
                if cRefSeq.sNMID[:2] == "NM" and "_" not in cRefSeq.sChrID:     #If NMID has no "NM" and chromosomeID has "-":
                        cRefNM.append(cRefSeq)
        print('Answer2 : ', str(len(cRefNM)))

        nDic = {}
        for cRefSeq in cRefNM:  #in cRefNM
                try:
                        nDic[cRefSeq.sNMID] += 1
                        cRefMultiNM.append(cRefSeq.sNMID)
                except KeyError:
                        nDic[cRefSeq.sNMID] = 1

        for cRefSeq in cRefNM:
                if cRefSeq.sNMID not in  cRefMultiNM:
                        cRefNewNM.append(cRefSeq)
        print ('Answer3 : ', str(len(cRefNewNM)))

        cRefNewNM.sort(key = lambda x: int(x.sNMID[3:]))
#       for cRefSeq in cRefNewNM[:10]:
#               print (cRefSeq.sNMID, cRefSeq.sGeneSym)
########################################################################################################################################################## Mission 3
        for cRefSeq in cRefNewNM:
                cRefSeq.GetCodon(cRefSeq.ExonStart, cRefSeq.ExonEnd, cRefSeq.ncodStart, cRefSeq.ncodEnd, cRefSeq.sChrStrand)
                lcRef.append(cRefSeq)

        for cRefSeq in lcRef:   #remove NM SEq that have wrong ORFs
                if len(cRefSeq.ORF) % 3 == 0 and cRefSeq.ORF[:3] == "ATG" and cRefSeq.ORF[-3:] in ["TAG", "TGA", "TAA"] and cRefSeq.stopcodon == 0:
                                cRefORF.append(cRefSeq)
        print ('Answer4 : ', str(len(cRefORF)))

        for cRefSeq in cRefORF: #remove the multiple same sGeneSym
                if cRefSeq.sGeneSym not in Isoform:
                        Isoform.append(cRefSeq.sGeneSym)
                        cRefUnique.append(cRefSeq)
        print ('Answer5 : ', str(len(cRefUnique)))


        cRefUnique.sort(key = lambda x: int(x.sNMID[3:]))
#       for i in cRefUnique[:10]:
#               print (str(i.sNMID), '\t', str(i.sGeneSym), '\t',  str(len(i.UTR5)), '\t', str(len(i.ORF)), '\t', str(len(i.UTR3)))

        for unique in cRefUnique:
                p = unique.sGeneSym
                unique.Get_Regulation(p)        #Calling function "Get_Regulation()"
                cUnique.append(unique)
###################################################################################################################################################################### Mission 4 & 5
def permutations(items, n):
    if n == 0:
        yield ''
    else:
        for i  in range(len(items)):            #permute the nucleotides
            for base in permutations(items, int(n-1)):  #call permutations starting ith the initial n passed until we reach 0.
                yield str(items[i] + str(base))

MotifList       = []    #Get 7mer Motif DNA List
Nucleotides     = ['A','C','G','T']
for motif in permutations(Nucleotides, 7):      #make 7mer seq
        MotifList.append(motif)


class OverRepresented_Motif():
        def __init__(self):
                self.Motif = 'NULL'
                self.n1 = 0
                self.n2 = 0
                self.n3 = 0
                self.n4 = 0
                self.fFractioA  = 0.0
                self.fFractionB = 0.0
                self.Enrichment = 0.0
                self.fPValue = 0.0
                self.cPValue = 0.0

        ## Find Over-represented Motifs
        def Get_RegulatedMotif(self, sReadMotif):
                #for i in MotifList:
                self.Motif = sReadMotif
                for u in cUnique:
                        if u.Regulated != 'NULL':               #Taking off Null
                                if sReadMotif in u.ORF:         #If there are motif:
                                        if u.Regulated == 'HighlyDown Regulated':
                                                self.n1 += 1
                                        else:
                                                self.n3 += 1
                                else:                           #If no motifs exist:
                                        if u.Regulated == 'HighlyDown Regulated':
                                                self.n2 += 1
                                        else:
                                                self.n4 += 1

                PSEUDO_COUNT = 0.0000001
                self.fFractionA = float((self.n1+PSEUDO_COUNT)/(self.n1+self.n2+PSEUDO_COUNT))
                self.fFractionB = float((self.n3+PSEUDO_COUNT)/(self.n3+self.n4+PSEUDO_COUNT))
                self.Enrichment = float(self.fFractionA/self.fFractionB)
                OddsRatio, PValue = stats.fisher_exact([[self.n1,self.n2],[self.n3,self.n4]])
                self.fPValue = PValue


lUniqueMotif = []
mData = []
def Call_RegulMotif(MotifList):
        for unique in MotifList:
                OverMotif = OverRepresented_Motif()
                OverMotif.Get_RegulatedMotif(unique)
                lUniqueMotif.append(OverMotif)

        maturefa = open('mature.fa', 'r')
        miRNAs = {}
        miRNASeq = maturefa.readlines()
        mName = ''
        mSeq = ''
        for i in miRNASeq:
                if i.startswith('>'):
                        mName = (i.lstrip('>').split(' '))[0]
                else:
                        mSeq = i.rstrip('\n')
                miRNAs[mName] = mSeq

        temp = []
        for key, value in miRNAs.items():
                temp = [key, value]
                mData.append(temp)

class miRNA_investigator:
        def __init__(self):
                self.miRname    = 'Null'
                self.miRseq     = 'Null'
                self.miRNAs_6mer = 'Null'
                self.miRNA_7mer_m8 = 'Null'
                self.miRNA_7mer_a1 = 'Null'
                self.miRNA_8mer = 'Null'
        def SetMiRNAs(self, sReadmiR):
                self.miRname    = sReadmiR[0]
                self.miRseq    = sReadmiR[1]
                self.miRNAs_6mer        = Reverse_Complement(self.miRseq[0:6])
                self.miRNAs_7mer_m8     = Reverse_Complement(self.miRseq[1:8])
                self.miRNAs_7mer_a1     = Reverse_Complement(self.miRseq[0:7])
                self.miRNAs_8mer        = Reverse_Complement(self.miRseq[0:8])

lmiRNA = []
def Call_miRNA(mData):
        for i in mData:
                if i[0].startswith('hsa'):
                        cmiR = miRNA_investigator()
                        cmiR.SetMiRNAs(i)
                        lmiRNA.append(cmiR)

def Sort_Enrichment(lUniqueMotif):
        lUniqueMotif.sort(key = lambda x: Decimal(str(x.fPValue)))
        print ("########################## Mission5 Dataset_1 - UTR3 #################################")
        numberofEnrich = 0
        for i in lUniqueMotif:
                if Decimal(str(i.Enrichment)) > Decimal('1.0'):
                        numberofEnrich += 1

        print (numberofEnrich)
        for i in lUniqueMotif[:5]:
                if Decimal(str(i.Enrichment)) > Decimal('1.0'):
                        for m in lmiRNA:
                                if i.Motif[6] == "A":
                                        if i.Motif == m.miRNAs_6mer:
                                                print (m.miRname, " : ", "6mer")
                                        if i.Motif[:6] == m.miRNAs_7mer_a1[:6]:
                                                print (m.miRname, " : ", "7mer_a1")
                                        if i.Motif == m.miRNAs_7mer_a1:
                                                print (m.miRname, " : ", "7mer_a1")
                                        if i.Motif == m.miRNAs_7mer_m8:
                                                print (m.miRname, " : ", "7mer_m8")
                                else:
                #                        if i.Motif == m.miRNAs_7mer_a1:
                 #                               print (m.miRname, " : ", "7mer_a1")
                                        if i.Motif == m.miRNAs_7mer_m8:
                                                print (m.miRname, " : ", "7mer_m8")
                        print (i.Motif, str(i.fPValue*numberofEnrich), str(i.n1), str(i.n2), str(i.n3), str(i.n4), str(i.Enrichment))
                        print ("   ")
        print ("######################################################################################")

main()
os.chdir(cwd)
