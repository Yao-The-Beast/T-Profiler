import csv
import sys
import NameConversion

# Stroe a list of genes that affected by a certain TF
# TF: the TF we are looking at
# targetGenes: [genes] (in commonName) that are affected by a certain TF
class TargetedGenes:
    def __init__(self, TF, targetGenes):
        self.TF = TF
        self.targetGenes = targetGenes

    def display(self):
        print(self.TF)
        print(self.targetGenes)

    def addATargetGene(self, thisGene):
        self.targetGenes.append(thisGene)



# we return TF (commonName), [genes (commonName)]
def processTargetGene(filename, nameDict):
    # read file
    with open(filename) as dataFile:
        content = dataFile.readlines()
    content = [x.strip() for x in content]

    targetGeneLists = []

    for line in content:
        entries = line.split(",")
        tg = TargetedGenes("",[])
        for (index, col) in enumerate(entries):
            value = ""
            # convert to commonName if possible
            if col != "":
                value = NameConversion.systematic2common(col, nameDict)
            # the first one is TF
            if index == 0:
                setattr(tg, "TF", value)
            # others are genes
            else:
                if col != "":
                    tg.addATargetGene(value)
        targetGeneLists.append(tg)

    return targetGeneLists
