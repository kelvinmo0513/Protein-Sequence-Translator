import sys
import os
import math

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TGC':'C', 'TGT':'C', 
    'TGG':'W', 'CGN':'R', 'GGN':'G', 'ACN':'T',
    'CCN':'P', 'CTN':'L', 'TCN':'S', 'GCN': 'A',
    'GTN':'V'}

def getSeq(inputFile):
    seqFile = os.path.join(os.getcwd(), inputFile)
    fileReader = open(seqFile, 'r')
    seq = ""
    for line in fileReader.readlines():
        if '>' in line:
            continue
        else:
            seq = seq + line.rstrip('\n')
    return seq

def findFrame(position):
    finalpos = position + 1
    if finalpos%3 == 0:
        return 3
    else:
        if finalpos%2 == 0:
            return 2
        else:
            return 1

# Method to get reverse complement of strand
def reverseComp(sequence):
    revseq = sequence[::-1]
    revcomp = ""
    for i in range(len(revseq)):
        curBase = revseq[i]
        if curBase == 'A':
            revcomp = revcomp + 'T'
        if curBase == 'C':
            revcomp = revcomp + 'G'
        if curBase == 'G':
            revcomp = revcomp + 'C'
        if curBase == 'T':
            revcomp = revcomp + 'A'
    return revcomp

def findORFs(sequence):
    reversecomp = reverseComp(sequence)
    print("Reverse Complement Found!")
    print(len(sequence))
    orfs = []
    for i in range(3):
        thisFrame = i + 1
        cur = i
        while(cur < len(sequence)):
            if sequence[cur:cur+3] == "ATG":
                start = cur
                cur = cur + 3
                while(cur < len(sequence)):
                    if sequence[cur:cur+3] == "TAA" or sequence[cur:cur+3] == "TGA" or sequence[cur:cur+3] == "TAG":
                        #print("orf found")
                        
                        orf = [start + 1, cur + 2, thisFrame]
                        #print(orfs)
                        orfs.append(orf)
                        cur = cur + 3
                        break
                    cur = cur + 3
                continue
            else:
                cur = cur + 3
    for i in range(3):
        thisFrame = 0
        if i == 0 or i == 2:
            if i == 0:
                thisFrame = 3
            else:
                thisFrame = 1
        else:
            thisFrame = 2
        cur = i
        while(cur < len(reversecomp)):
            if reversecomp[cur:cur+3] == "ATG":
                start = cur
                cur = cur + 3
                while(cur < len(reversecomp)):
                    if reversecomp[cur:cur+3] == "TAA" or reversecomp[cur:cur+3] == "TGA" or reversecomp[cur:cur+3] == "TAG":
                        orf = [len(reversecomp) - (cur) + 3, len(reversecomp) - start, -1*(thisFrame)]
                        orfs.append(orf)
                        cur = cur + 3
                        break
                    cur = cur + 3
                continue
            else:
                cur = cur + 3
    print(len(orfs))
    return orfs


def filterORFs(orfs):
    filteredORFs = []
    for orf in orfs:
        if (orf[1] - orf[0] + 1) >= 90:
            filteredORFs.append(orf)
    return filteredORFs


def writeORFSeq(filteredORFs, sequence, filename):
    outputFile = os.path.join(os.getcwd(), filename)
    fileWriter = open(outputFile, 'w')
    cur = 1
    for orf in filteredORFs:
        fileWriter.write(">Found ORF #" + str(cur) + '\n')
        if orf[2] < 0:
            fileWriter.write(getProt(reverseComp(sequence[orf[0]-1:orf[1]])) + '\n')
        else:
            fileWriter.write(getProt(sequence[orf[0]-1:orf[1] + 1]) + '\n')
        cur = cur + 1
    fileWriter.close()

def getProt(orf):
    orflen = len(orf)
    protein = ""
    cur = 0
    while cur < orflen:
        if orf[cur:cur+3] == 'TAG' or orf[cur:cur+3] == 'TAA' or orf[cur:cur+3] == 'TGA' or len(orf[cur:cur+3]) < 3:
            return protein
        protein = protein + gencode[orf[cur:cur+3]]
        cur = cur + 3

def getFile(sequenceFile):
    outputFile = ''
    for i in range(len(sequenceFile)):
        if sequenceFile[i:i+3] == 'fna':
            return outputFile + "txt"
        outputFile = outputFile + sequenceFile[i]

# Main method
def main():
    seqfile = sys.argv[1]
    #print(seqfile)
    sequence = getSeq(seqfile)
    prelimORFs = findORFs(sequence)
    outputName = getFile(seqfile)
    print("Preliminary ORFs Found!")
    filteredORFs = filterORFs(prelimORFs)
    print("ORFs filtered!")
    print(outputName)
    writeORFSeq(filteredORFs, sequence, outputName)
    print("ORFs written to file!")

if __name__ == '__main__':
    main()
