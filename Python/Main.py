import NameConversion
import TargetGene
import GeneExpressionData
import T_Profiler

# return a dict
# key: value = systematicName: commonName
nameDict = NameConversion.processDict("../Calling_Cards_data/YeastCommonAndSystematicGeneNames.csv")

# return a list of targetGene Object
# TF                string
# targetGenes       [string]
# all the name are commonName
targetGeneList = TargetGene.processTargetGene("../Calling_Cards_data/targetGenes.csv", nameDict)
# for entry in targetGeneList:
#     entry.display()

# read data from Gene
experiments = GeneExpressionData.processExperiments("../Calling_Cards_data/callingCardsPlates.csv", nameDict, targetGeneList)
GeneExpressionData.summarizeExperiments(experiments, "output/Experiments_Summary.txt")

baseCaseExpIDs = [1,2,3]

baseExperiments, targetExperiments = GeneExpressionData.divideExperiments_base_target(experiments, baseCaseExpIDs)
baseExperiment = baseExperiments[0]

results = []
# generate T values for experiments
for (index, thisExperiment) in enumerate(targetExperiments):
    thisTProfileCalculator = T_Profiler.TProfilerCalculator(thisExperiment, baseExperiment)
    thisTProfileCalculator.calculate_t_values(targetGeneList, nameDict)

    # extract the information from TProfilerCalculator Object to TProfilerResult Object
    for (key, value) in thisTProfileCalculator.t.items():
        thisT = thisTProfileCalculator.t[key]
        thisP = thisTProfileCalculator.p[key]
        thisTProfilerResult = T_Profiler.TProfilerResult(thisTProfileCalculator.targetExperiment.id, key, thisTProfileCalculator.targetExperiment.description, thisT, thisP)
        results.append(thisTProfilerResult)

# sort by absolute t value
# And print out the TProfilerResult
results.sort(key=lambda x: float('-inf') if x.experiment_t == None else abs(x.experiment_t) , reverse=True)
for thisResult in results:
    thisResult.display()
