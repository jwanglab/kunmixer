import argparse
import sys, glob, json, os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import pandas as pd
sns.set(style="white", font_scale=1.2)


inputPath = "./kunmixer"

def plotMatchRateDistribution(matchRate, *threshold_rate):
    plt.figure(figsize=(12, 6))
    plotdf = pd.DataFrame({"matchRate":matchRate})
    if threshold_rate:
        hue = []
        for rate in matchRate:
            if rate >=threshold_rate[0]:
                hue.append(">" + str(threshold_rate[0]))
            else:
                hue.append("<" + str(threshold_rate[0]))
        plotdf["type"] = hue
        fig, ax = plt.subplots(1, figsize=(12,10) , tight_layout=True)
        fig.suptitle('Similarity Rate Distribution')

        sns.histplot(ax=ax, data = plotdf, x="matchRate", bins=100, color="dodgerblue")
        ax.set_title("All Similarity Rates")
        ax.set_ylabel("Number of Similarity Rates")

        plt.xlabel("Similarity Rate")
        plt.ylabel("Number of Similarity Rates")
        outputName = inputPath + "_matchRate_" + str(threshold_rate[0]) + "_distribution.png"
        plt.savefig(outputName, dpi=400, bbox_inches='tight')
    else:
        sns.histplot(data = plotdf, x="matchRate", rug=True, kde=False, bins=100)
        outputName = inputPath + "_matchRate_distribution.png"
        plt.savefig(outputName, dpi=400, bbox_inches='tight')

def plotCompareMatrix(matchMatrix, n):
    matplotlib.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=(15,15), dpi=100)
    ax = plt.gca()
    im = plt.imshow(matchMatrix, interpolation="none", vmin=0, vmax=1)
    plt.yticks(range(matchMatrix.shape[0]), names[:n])
    plt.xticks(range(matchMatrix.shape[1]), names[:n], rotation="vertical")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    plt.tight_layout()
    plt.savefig("compare_matrix.png")

def buildFileFingerPrint(files):
    genotypes = []

    for i in range(len(files)):
        print("Processing {}".format(files[i].split("/")[-1]))

        data = [line.strip().split('\t') for line in open(files[i])]
        data = data[1:] # first line is header
        genotypeArray = np.zeros(len(data), dtype="u1")

        two_als = 0
        toomany_als = 0
        one_als = 0

        for d in range(len(data)):
            oneLine = data[d]
            alleles = [int(a) for a in oneLine[4:8]]
            maxCount = max(alleles)
            nalleles = 0
            for bit in range(len(alleles)):
                if (alleles[bit] > maxCount/4 and maxCount > 16) or (alleles[bit] > 6 and maxCount <=16):
                    genotypeArray[d] = genotypeArray[d] | (1<<bit)
                    nalleles += 1
            if nalleles == 2:
                two_als += 1
            elif nalleles > 2:
                toomany_als += 1
            elif nalleles == 1:
                one_als += 1
        print("=1 alleles: {}".format(one_als))
        print("=2 alleles: {}".format(two_als))
        print(">2 alleles: {}".format(toomany_als))
        genotypes.append(genotypeArray)
    return genotypes

def buildMatrix(genotypes):
    matchRate = []
    n = len(genotypes)
    matchMatrix = np.zeros((n, n))
    totalMatrix = np.zeros((n, n))
    rateMatrix = np.zeros((n, n))

    for i in range(n):
        rateMatrix[i,i] = 1
        if genotypes[i].shape[0] == 0 or np.count_nonzero(genotypes[i]) < 100:
            continue
        for j in range(i+1, n):
            if genotypes[j].shape[0] == 0 or np.count_nonzero(genotypes[j]) < 100:
                continue
            if genotypes[i].shape != genotypes[j].shape:
                print("{} ({}) and {} ({}) not the same shape".format(i, genotypes[i].shape[0], j, genotypes[j].shape[0]))
                continue
            # oneone means num of alleles are equal and both > 0
            oneone = np.count_nonzero(np.logical_and(np.equal(genotypes[i], genotypes[j]), np.logical_or(genotypes[i], genotypes[j])))
            total = np.count_nonzero(np.logical_and(genotypes[i], genotypes[j]))
            if total == 0:
                continue
            rate = oneone/total
            matchRate.append(rate*100)

            matchMatrix[i, j] = oneone
            matchMatrix[j, i] = matchMatrix[i, j]

            totalMatrix[i, j] = total
            totalMatrix[j, i] = totalMatrix[i, j]

            rateMatrix[i, j] = rate
            rateMatrix[j, i] = rateMatrix[i, j]

    return matchRate, matchMatrix, totalMatrix, rateMatrix

def keepdata(name, match, total):
    datakeeper[name] = {}
    datakeeper[name]['Same'] = match.tolist()
    datakeeper[name]['SharedSNPSites'] = total.tolist()


def generateCompareResult(genotypes, matchMatrix, totalMatrix, rateMatrix, names, minRate=0.6, goodRate=0.75):
    outputFilename = inputPath + "_SNP_check_result.csv"
    file = open(outputFilename, "w")
    file.write("Similarity,Status,Target,NumSNPs1,BestMatched,NumSNPs2,Same,SharedSNPSites,Notes\n")

    for i in range(len(names)):

        rateMatrix[i,i] = 0
        j = np.argmax(rateMatrix[i,:])
        maxTotal = totalMatrix[i, j]
        HighRate = rateMatrix[i, j]

        secHighIndex = np.argsort(-rateMatrix[i, :])[1]
        #secTotal = totalMatrix[i, secHighIndex]
        secHigh = rateMatrix[i, secHighIndex]
        rateMatrix[i,i] = 1

        if HighRate > minRate:
            if secHigh > goodRate:
                keepdata(names[i].split(".ktype")[0], matchMatrix[i,:], totalMatrix[i,:])
                addition = 'Second Highest: ' + names[secHighIndex].split(".ktype")[0] + ' ' + "{:.2f}".format(secHigh*100) + '; NumOfSameSNPsites: 2nd 3rd 4th 5th : ' + ' '.join("{:.0f}".format(a) for a in sorted(matchMatrix[i,:])[-2:-6:-1])
                file.write("{:.2f},MUlTI.GOOD,{},{},{},{},{},{},{}\n".format(HighRate * 100,names[i].split(".ktype")[0], np.count_nonzero(genotypes[i]), names[j].split(".ktype")[0], np.count_nonzero(genotypes[j]), int(matchMatrix[i,j]), int(maxTotal), addition))
            elif (HighRate - secHigh > 0.1 or HighRate > goodRate):
                keepdata(names[i].split(".ktype")[0], matchMatrix[i,:], totalMatrix[i,:])
                addition = 'Second Highest: ' + names[secHighIndex].split(".ktype")[0] + ' ' + "{:.2f}".format(secHigh*100) + '; NumOfSameSNPsites: 2nd 3rd 4th 5th : ' + ' '.join("{:.0f}".format(a) for a in sorted(matchMatrix[i,:])[-2:-6:-1])
                file.write("{:.2f},GOOD,{},{},{},{},{},{},{}\n".format(HighRate * 100,names[i].split(".ktype")[0], np.count_nonzero(genotypes[i]), names[j].split(".ktype")[0], np.count_nonzero(genotypes[j]), int(matchMatrix[i,j]), int(maxTotal), addition))
            else:
                keepdata(names[i].split(".ktype")[0], matchMatrix[i,:], totalMatrix[i,:])
                addition = 'Second Highest: ' + names[secHighIndex].split(".ktype")[0] + ' ' + "{:.2f}".format(secHigh*100) + '; NumOfSameSNPsites: 2nd 3rd 4th 5th : ' + ' '.join("{:.0f}".format(a) for a in sorted(matchMatrix[i,:])[-2:-6:-1])
                file.write("{:.2f},UNCERTAIN,{},{},{},{},{},{},{}\n".format(HighRate * 100,names[i].split(".ktype")[0], np.count_nonzero(genotypes[i]), names[j].split(".ktype")[0], np.count_nonzero(genotypes[j]), int(matchMatrix[i,j]), int(maxTotal), addition))

    file.close()

def generateNonDuplicateRecord():
	dt = pd.read_csv(inputPath + "_SNP_check_result.csv")
	manifest = []
	indexlist = []

	for index, row in dt.iterrows():
	    capsule = sorted([row['Target'], row['BestMatched']])
	    if capsule not in manifest:
	        indexlist.append(index)
	        manifest.append(capsule)
	unique = dt.iloc[indexlist,:]
	unique.to_csv(inputPath + "_SNP_check_result.csv", index=False)

def generateSimpleSampleFingerPrintCsv(names, genotypes):
    names = [x.split(".ktype")[0] for x in names]
    dt = pd.DataFrame(columns=names)
    for index in range(len(names)):
        dt[names[index]] = genotypes[index]
    dt.to_csv("SimpleSampleFingerPrint.csv", index=False)

def calculateAverageMatchRate():
    df = pd.read_csv(inputPath + "_SNP_check_result.csv")
    print("Avg similarity rate is: {:.2f}".format(sum(df['Similarity'].tolist())/df.shape[0]))
    print("The 5 smallest similarity rates are: {}".format(str(sorted(df['Similarity'].tolist())[:5])))
    print("The 5 smallest match numbers are: {}".format(str(sorted(df['Same'].tolist())[:5])))
    print("The similarity output files are: " + inputPath + "_SNP_check_result.csv and " + inputPath + "_matchRate_55_distribution.png")
    print("The fingerprint output of all samples is in current path: SimpleSampleFingerPrint.csv")

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Compare kunmixer profiles across samples")
    parser.add_argument("profiles", help="Kunmixer (ktype or fastttype) profile files", nargs="+")
    args = parser.parse_args()
    files = args.profiles
    files = [f for f in files if os.stat(f).st_size > 0]
    print(f"Using {len(files)} of {len(args.profiles)} files, the rest are empty...")
    names = [f.split('/')[-1][:-6] for f in files]

    # RNAseq(FF, FFPE)
    # WXS, WGS
    threshold_rate = 0.55
    match_rate = 0.75

    # build snp sites count array for each sample
    genotypes = buildFileFingerPrint(files)
    matchRate, matchMatrix, totalMatrix, rateMatrix = buildMatrix(genotypes)

    # plot data distribution and compare matrix
    plotMatchRateDistribution(matchRate, int(threshold_rate * 100))
    plotCompareMatrix(rateMatrix, len(names))

    # find and generate match pairs
    datakeeper = {}
    generateCompareResult(genotypes, matchMatrix, totalMatrix, rateMatrix, names, threshold_rate, match_rate)

    # remove the duplicated records from newly created file
    generateNonDuplicateRecord()

    # generate simple sample fingerprint for later snp screening
    generateSimpleSampleFingerPrintCsv(names, genotypes)

    calculateAverageMatchRate()
