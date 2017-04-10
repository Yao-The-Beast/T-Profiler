import csv
import sys
import re
import NameConversion
from math import log

# Experiment object
# TF                        Transcription Factors [string]
# genes_expressionValue     dict: GeneName (commonName): ExpressionValue (string : float)
# condition                 describe the experiment (String)
# id                        the column number in the database (csv)
# description               the original x label
class Experiment:
    def __init__(self, TF, genes_expressionValue, condition, description,id):
        self.TF = TF
        self.genes_expressionValue = genes_expressionValue
        self.condition = condition
        self.id = id
        self.description = description

    def displayVerbose(self, randomGeneName, dictName):
        print("TF:")
        print(self.TF)
        print("Condition:")
        print(self.condition)
        print("ID:")
        print(self.id)
        print("Description:")
        print(self.description)
        print(randomGeneName + " Gene Expression Value")
        print(self.genes_expressionValue[NameConversion.systematic2common(randomGeneName, dictName)])
        print("----------------------------")

    def display(self):
        print("ID:")
        print(str(self.id))
        print("TF:")
        print(self.TF)
        print("Condition:")
        print(self.condition)
        print("Description:")
        print(self.description)
        print("-------------------------------")

    def addAGeneExpressionValue(self, thisGene, thisExpressionValue):
        self.genes_expressionValue[thisGene] = thisExpressionValue

    def getThisGeneExpressionValue(self, thisGene):
        return self.genes_expressionValue[thisGene]

    def addATF(self, thisTF):
        self.TF.append(thisTF)





# helper of processExperiments
# return a list of experiments [Experiment]
# TF are not extracted from experiment.condition
# experiment.genes_expressionValue is complete
def processExperiments_1(filename, nameDict):
    # read file
    with open(filename) as dataFile:
        content = dataFile.readlines()
    content = [x.strip() for x in content]

    # the experiments we interested in
    experiments = []
    # process first row
    firstRow = content[0]
    entries = firstRow.split(",")
    for (index, col) in enumerate(entries):
        # we omit CONTROL-CONTROL and the first entry
        # if "CONTROL-CONTROL" not in col and index != 0:
        #     thisExperiment = Experiment([], dict(), col, col, index)
        #     experiments.append(thisExperiment)
        if index != 0:
            thisExperiment = Experiment([], dict(), col, col, index)
            experiments.append(thisExperiment)

    # process the whole matrix
    for line in content:
        entries = line.split(",")
        # first we check the y column is the gene we want
        # the gene we want start withs Y
        thisGene = ""
        match = re.search('(^Y\S+)', entries[0])
        if match:
            thisGene = entries[0]
            thisGene = NameConversion.systematic2common(thisGene, nameDict)
        else:
            continue
        # then we loop through the valid experiments
        # add the {gene, geneExpressionValue} pair to each experiment.genes_expressionValue
        for (index, thisExperiment) in enumerate(experiments):
            thisCol = thisExperiment.id
            thisGeneExpressionValue = entries[thisCol]
            if (float(thisGeneExpressionValue) != 0):
                thisGeneExpressionValue = float(thisGeneExpressionValue)
            else:
                thisGeneExpressionValue = 0
            experiments[index].addAGeneExpressionValue(thisGene, thisGeneExpressionValue)
    return experiments


# helper of the processExperiments
# input:  a list of experiments [experiment] from processExperiments_1
# return complete experiments
# TF and conditions are extracted from the X label
# the output [Experiment] are ready to use
def processExperiments_2(experiments, nameDict, targetGeneList):

    # get a list of TFs
    TF = []
    for thisTargetGeneObj in targetGeneList:
        TF.append(thisTargetGeneObj.TF)

    # iterate through the experiments
    for thisExperiment in experiments:
        thisDescription = thisExperiment.description
        # condition delimiter is "-"
        entries = thisDescription.split("-")
        for block in entries:
            # MB0479_EDS1_RGT1-plusLys-2-03.04.16-04.14.16-13
            # split each block further more using _
            for col in block.split("_"):
                # first we check if the col is a systematicName (eg: YBR033W)
                # if it is, we convert col to commonName
                # the commonName should be one of the entry in TF[]
                if NameConversion.isInDict(col, nameDict):
                    thisTFCommonName = NameConversion.systematic2common(col, nameDict)
                    if thisTFCommonName not in TF:
                        print("TF: " + thisTFCommonName + " is not defined in TargetGene TF list")
                    thisExperiment.addATF(thisTFCommonName)
                # if not, col might already be a commonName (eg: EDS1)
                # so we check whether col is already contained in the TF[]
                elif col in TF:
                    thisExperiment.addATF(col)
                # if it is neither a commonName or a systematicName,
                # we check if it describes the environment (eg: plusLys, minusLys)
                else:
                    plusMatch = re.search('(^plus\S+)', col)
                    minusMatch = re.search('(^minus\S+)', col)
                    if plusMatch or minusMatch:
                        thisExperiment.condition = col
    return experiments


# USE THIS FUNCTION !
# Output a txt file that contains the description of all the experiments
# it is done by calling display(experiments)
def summarizeExperiments(experiments, outputFileName):
    f = open(outputFileName, 'w')
    for thisExp in experiments:
        outputLine = ""
        outputLine += ("ID: " + str(thisExp.id) + " \n")
        outputLine += ("Description: " + str(thisExp.description) + " \n")
        outputLine += "---------------------------- \n"
        f.write(outputLine)
    f.close()

# USE THIS FUNCTION !!
# Arguments: dataFileName, systematic&common name dictionary, [targetGene Object]
# return [Experiment Object]
def processExperiments(dataFileName, nameDict, targetGeneList):
    experiments = processExperiments_1(dataFileName, nameDict)
    experiments = processExperiments_2(experiments, nameDict, targetGeneList)
    return experiments

# USE THIS FUNCTION !!
# Arguments: allExperiments, [baseCaseExperimentsID]
# We extract the base case experiments from all the experiments based on the ID
# Then we average the value among these experiments to generate an averaged base experiment
# Return two parts [averagedBaseExperiment] [targetExperiment]
def divideExperiments_base_target(experiments, baseCaseExpIDs):
    baseExperiments = []
    targetExperiments = []
    # find the base case experiments based on IDs
    for thisExperiment in experiments:
         if thisExperiment.id in baseCaseExpIDs:
             baseExperiments.append(thisExperiment)
         else:
             targetExperiments.append(thisExperiment)

    # average the base experiments
    averagedBaseExperiment = Experiment([], {}, "", "", -1)
    for (index, thisBaseExperiment) in enumerate(baseExperiments):
        if index == 0:
            averagedBaseExperiment.genes_expressionValue = thisBaseExperiment.genes_expressionValue
            averagedBaseExperiment.description = thisBaseExperiment.description  + ";\n"
        else:
            for (key, value) in thisBaseExperiment.genes_expressionValue.items():
                averagedBaseExperiment.genes_expressionValue[key] += value
            averagedBaseExperiment.description += (thisBaseExperiment.description +  ";\n")

    for (key, value) in averagedBaseExperiment.genes_expressionValue.items():
        averagedBaseExperiment.genes_expressionValue[key] = value / len(baseExperiments)

    return [averagedBaseExperiment], targetExperiments
