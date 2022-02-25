#!/usr/bin/env python
# merge 10x outs web info txt extracted from web.html

import os


def mergeInfo(filePathListFilePath, outFilePath):

    totalInfoDict = dict()
    totalKeys = set()

    with open(filePathListFilePath) as filePathListFile:
        for filePath in filePathListFile:
            filePath = filePath.strip()
            if filePath == "": continue

            sampleName = filePath[2:filePath.index("_analysis")]

            totalInfoDict[sampleName] = dict()

            with open(filePath) as txtFile:
                for line in txtFile:
                    line = line.strip()
                    if line == "": continue

                    listLine = line.split("\t")
                    totalInfoDict[sampleName][listLine[0]] = listLine[1]
                    totalKeys.add(listLine[0])

    with open(outFilePath, 'w') as outFile:
        keyList = list(totalKeys)

        outFile.write("{0}\t{1}\n".format("SampleID", "\t".join(keyList)))

        for sample in totalInfoDict.keys():
            outFile.write("{}".format(sample))
            for samplekey in keyList:
                if samplekey in totalInfoDict[sample].keys():
                    outFile.write("\t{}".format(totalInfoDict[sample][samplekey]))
                else:
                    outFile.write("\tNA")
            outFile.write("\n")

def main():
    #  import argparse
    #  parser = argparse.ArgumentParser(description = 'software info')
    #  parse.add_argument()

    #  # parser.add_argument('--fqFilePathList', nargs='+', help="fqFilePathList")
    #  # parser.add_argument('--bowtie2outPath', dest='bowtie2outPath', help='bowtie2outPath')
    #  # parser.add_argument('--metaDataPath', dest='metaDataPath', help='metaDataPath')
    #  # parser.add_argument('--outFilePath', dest='outFilePath', help='outFilePath')

    #  args = parser.parse_args()

    filePathListFilePath = "./file.path.list"
    outFilePath = "info.merged.txt"

    mergeInfo(filePathListFilePath, outFilePath)

if __name__ == '__main__':
    main()
