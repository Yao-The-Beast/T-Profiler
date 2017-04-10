import NameConversion
import GeneExpressionData
import TargetGene
import statistics
import math
from scipy import stats
import numpy as np

# t profiler caculation class (objectd)
# Eg: for experiment: MB0479_EDS1_RGT1-plusLys-2-03.04.16-04.14.16-13
# We need to take EDS1 and RGT1 apart to view them separately.
# We need to look at EDS1 first, than RGT1
# U_G = Mean [ log( [the expression value of the genes influenced by EDS1 in this targetExperiment]
#           / [the expression value of those Genes in base case] ) ]
# U_G' = Mean [ log( [the expression value of the genes not influenced by EDS1 in this targetExperiment]
#           / [the expression value of those Genes in base case] ) ]
# Base case are the average of experiments that NO TFs are knocked out

class TProfilerCalculator:
    def __init__(self, targetExperiment, baseExperiment):
        self.baseExperiment = baseExperiment
        self.targetExperiment = targetExperiment
        self.t = dict()
        self.p = dict()

    def display(self):
        print("Target Experiment ID: " + str(self.targetExperiment.id))
        print("Target Experiment Description: ")
        print(self.targetExperiment.description)
        # print("Base Case Experiment:")
        # display(baseExperiment)
        print("T & P Value is: ")
        for thisKey in self.t.keys():
            print("t: " + str(self.t[thisKey]) + "; " + "p: " + str(self.p[thisKey]))
        print("---------------------------------")


    # calculate mean
    def calculate_mean(self, genesDict):
        values = [v for v in genesDict.values()]
        mean = statistics.mean(values)
        return mean

    # calculate_estimated standard deviation
    def calculate_estimated_variance(self, genesDict):
        values = [v for v in genesDict.values()]
        variance = statistics.pvariance(values)
        return variance

    # calculate pooled standard deviation s
    def calculate_pooled_standard_deviation(self, targetGenesDict, unTargetGenesDict):
        variance_target = self.calculate_estimated_variance(targetGenesDict)
        variance_untarget = self.calculate_estimated_variance(unTargetGenesDict)
        num_target = len(targetGenesDict)
        num_untarget = len(unTargetGenesDict)
        # s = sqrt( [(N - 1) * variance + (N' - 1) * variance'] / (N + N' - 2) )
        numerator = (num_target - 1.0) * variance_target + (num_untarget - 1.0) * variance_untarget
        denumerator = (num_target + num_untarget - 2.0)
        s = math.sqrt(float(numerator) / float(denumerator))
        return s

    def isBaseCase(self):
        if len(self.targetExperiment.TF) == 0:
            return True
        else:
            return False

    # yao's implementation of T_Profiler
    def calculate_t_values_yao_version(self, targetGeneList, nameDict):
        # if there is no TF existed in the experiment
        # we return empty list (as we are not able to generate t values)
        if self.isBaseCase():
            # print("Experiment ID " + str(self.targetExperiment.id) + " has no target TFs")
            return

        t_values = []
        allGenesDict = self.targetExperiment.genes_expressionValue
        allGenesBaseCaseDict = self.baseExperiment.genes_expressionValue
        for thisTF in self.targetExperiment.TF:
            # retrieve the genes affected by this TF from targetGeneLists
            # append it to the targetGenes
            targetGenes = []
            for thisTargetGeneObj in targetGeneList:
                if thisTargetGeneObj.TF == thisTF:
                    targetGenes = thisTargetGeneObj.targetGenes
                    break
            # separate the overall genes into two set
            # targetGenes vs unTargetGenes
            # unTargetGenes = allGenesDict.keys - targetGenes
            unTargetGenes = list(allGenesDict.keys())
            for thisTargetGene in targetGenes:
                if thisTargetGene in unTargetGenes:
                    unTargetGenes.remove(thisTargetGene)
            # generate targetGenesDict & unTargetGenesDict
            # {key : Value = Gene: ExpVal_in_target_exp / ExpVal_in_base_exp}
            targetGenesDict = dict()
            unTargetGenesDict = dict()
            # calculate log expression ratio
            for thisGene in targetGenes:
                if allGenesBaseCaseDict[thisGene] != 0 and allGenesDict[thisGene] != 0:
                    targetGenesDict[thisGene] = math.log(float(allGenesDict[thisGene]) / allGenesBaseCaseDict[thisGene])
            for thisGene in unTargetGenes:
                if allGenesBaseCaseDict[thisGene] != 0 and allGenesDict[thisGene] != 0:
                    unTargetGenesDict[thisGene] =  math.log(float(allGenesDict[thisGene]) / allGenesBaseCaseDict[thisGene])

            # calculate mean log expression value ratio
            targetGeneMean = self.calculate_mean(targetGenesDict)
            untargetGeneMean = self.calculate_mean(unTargetGenesDict)
            # calculate pooled standard deviation s
            s = self.calculate_pooled_standard_deviation(targetGenesDict, unTargetGenesDict)
            num_target = len(targetGenesDict)
            num_untarget = len(unTargetGenesDict)
            # calculate t
            # t = (mean_target - mean_untarget) / [s * sqrt(1/num_target + 1/num_untarget)]
            t = (targetGeneMean - untargetGeneMean) / float(s * math.sqrt(float(1.0/num_target) + float(1.0/num_untarget)))

            # append the result to t_values
            self.t[thisTF] = t

    # USE THIS FUNCTION
    # Calculate T values
    # Outputs are stored in the self.t dictionary
    def calculate_t_values(self, targetGeneList, nameDict):
        # if there is no TF existed in the experiment
        # we return empty list (as we are not able to generate t values)
        if self.isBaseCase():
            # print("Experiment ID " + str(self.targetExperiment.id) + " has no target TFs")
            return

        t_values = []
        allGenesDict = self.targetExperiment.genes_expressionValue
        allGenesBaseCaseDict = self.baseExperiment.genes_expressionValue
        for thisTF in self.targetExperiment.TF:
            # retrieve the genes affected by this TF from targetGeneLists
            # append it to the targetGenes
            targetGenes = []
            for thisTargetGeneObj in targetGeneList:
                if thisTargetGeneObj.TF == thisTF:
                    targetGenes = thisTargetGeneObj.targetGenes
                    break
            # separate the overall genes into two set
            # targetGenes vs unTargetGenes
            # unTargetGenes = allGenesDict.keys - targetGenes
            unTargetGenes = list(allGenesDict.keys())
            for thisTargetGene in targetGenes:
                if thisTargetGene in unTargetGenes:
                    unTargetGenes.remove(thisTargetGene)
            # generate targetGenesDict & unTargetGenesDict
            # {key : Value = Gene: ExpVal_in_target_exp / ExpVal_in_base_exp}
            targetGenesDict = dict()
            unTargetGenesDict = dict()
            # calculate log expression ratio
            for thisGene in targetGenes:
                if allGenesBaseCaseDict[thisGene] != 0 and allGenesDict[thisGene] != 0:
                    targetGenesDict[thisGene] = math.log(float(allGenesDict[thisGene]) / allGenesBaseCaseDict[thisGene])
            for thisGene in unTargetGenes:
                if allGenesBaseCaseDict[thisGene] != 0 and allGenesDict[thisGene] != 0:
                    unTargetGenesDict[thisGene] =  math.log(float(allGenesDict[thisGene]) / allGenesBaseCaseDict[thisGene])

            # scipy two sample t test
            numpyArray1 = np.fromiter(iter(targetGenesDict.values()), dtype=float)
            numpyArray2 = np.fromiter(iter(unTargetGenesDict.values()), dtype=float)
            (t, p) = stats.ttest_ind(numpyArray1, numpyArray2, equal_var = True)

            # append the result to t_values
            self.t[thisTF] = t
            self.p[thisTF] = p


class TProfilerResult:
    def __init__(self, experiment_id, experiment_TF, experiment_description, experiment_t, experiment_p):
        self.experiment_id = experiment_id
        self.experiment_TF = experiment_TF
        self.experiment_description= experiment_description
        self.experiment_t = experiment_t
        self.experiment_p = experiment_p

    def display(self):
        print("Experiment ID: " + str(self.experiment_id))
        print("Target TF: " + str(self.experiment_TF))
        print("Experiment Description: " + str(self.experiment_description))
        print("T Value: " + str(self.experiment_t))
        print("P Value: " + str(self.experiment_p))
        print("----------------------------------")
