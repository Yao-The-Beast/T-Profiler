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
# summary the data
GeneExpressionData.summarizeExperiments(experiments, "output/Experiments_Summary.txt")
# normalize data
# experiments = GeneExpressionData.normalize_data(experiments, 0, 10)
# experiments = GeneExpressionData.normalize_data_2(experiments, [0])

#
# # here to specify base and target experiments
targetExpIDs = [4,5,6]
baseCaseExpIDs = [21,22,23]
testTFList = ["EDS1","RGT1","LYS14"]

baseExperiments, targetExperiments = GeneExpressionData.divideExperiments_base_target(experiments, baseCaseExpIDs, targetExpIDs)
baseExperiment = baseExperiments[0]
targetExperiment = targetExperiments[0]

# run the T_Profiler
results = []
thisTProfileCalculator = T_Profiler.TProfilerCalculator(targetExperiment, baseExperiment)
thisTProfileCalculator.calculate_t_values_paired_t_test(targetGeneList, nameDict, testTFList)

# extract the information from TProfilerCalculator Object to TProfilerResult Object
for (key, value) in thisTProfileCalculator.t.items():
    thisT = thisTProfileCalculator.t[key]
    thisP = thisTProfileCalculator.p[key]
    thisTProfilerResult = T_Profiler.TProfilerResult(thisTProfileCalculator.targetExperiment.id, key, thisTProfileCalculator.targetExperiment.description, thisT, thisP)
    results.append(thisTProfilerResult)

# sort by absolute t value
# And print out the TProfilerResult
results.sort(key=lambda x: float('-inf') if x.experiment_t == None else abs(x.experiment_t) , reverse=True)

# display base case and t-profiler output
print("Base Case is: ")
print(baseExperiment.description.split()[0])
print("TF KO:" + str(baseExperiment.TF))
print("=====================================================" + "\n")
print("T Profiler Output is: \n")
for thisResult in results:
    thisResult.display()
