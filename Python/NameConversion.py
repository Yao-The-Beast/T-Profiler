import csv
import sys


# key, value = systematicName, commonName
def processDict(filename):
    # read file
    with open(filename) as dataFile:
        content = dataFile.readlines()
    content = [x.strip() for x in content]
    content.pop(0)
    nameDict = dict()
    for line in content:
        entries = line.split(",")
        commonName = entries[0]
        systematicName = entries[1]
        nameDict[systematicName] = commonName
    return nameDict

def systematic2common(systematicName, dictName):
    if systematicName in dictName:
        return dictName[systematicName]
    else:
        return systematicName

def isInDict(systematicName, dictName):
    if systematicName in dictName:
        return True
    else:
        return False
