#!/usr/bin/env python
# stat markers

def makeMarker2CellTypeDict(markersDatabaseFilePath):
    # {marker:{set([celltype])}}
    marker2CellTypeDict = dict()

    with open(markersDatabaseFilePath) as inFile:
        recordIdx = 0
        for line in inFile:
            line = line.strip("\n")
            line = line.strip(" ")
            line = line.strip()
            if line == "":
                continue
            recordIdx = recordIdx + 1

            if recordIdx == 1:
                listLine = line.lower().split("\t")
                markerIdx = listLine.index("marker")
                cellTypeIdx = listLine.index("celltype")
                print("index:")
                print(markerIdx)
                print(cellTypeIdx)
                continue

            listLine = line.split("\t")
            marker = listLine[markerIdx].lower()
            print(listLine)
            cellType = listLine[cellTypeIdx]

            if marker not in marker2CellTypeDict.keys():
                marker2CellTypeDict[marker] = set([cellType])
            else:
                marker2CellTypeDict[marker].add(cellType)

    return marker2CellTypeDict

def run(markersTop, markersDatabase, outputFilePath):
    m2ctd = makeMarker2CellTypeDict(markersDatabase)
    with open(markersTop)as meFile:
        with open(outputFilePath, "w") as outputFile:
            for line in meFile:
                line = line.strip()
                if line == "": continue
                if "gene" in line:
                    outputFile.write(line + "\tcelltype\n")
                    continue
                listLine = line.split("\t")
                gene = listLine[-1].lower()
                if gene not in m2ctd.keys():
                    outputFile.write(line + "\t-\n")
                else:
                    outputFile.write(line + "\t{0}\n".format("|".join(m2ctd[gene])))

def main():
    import argparse
    parser = argparse.ArgumentParser(description = 'software info')
    parser.add_argument('--markersTop', dest='markersTop', help="markersTop")
    parser.add_argument('--markersDatabase', dest='markersDatabase', help='markersDatabase')
    parser.add_argument('--out', dest='out', help='out')
    args = parser.parse_args()
    run(args.markersTop, args.markersDatabase, args.out)

if __name__ == '__main__':
    main()
