#######VIRGO2
#######Authors: Michael France
#######Contact: mfrance@som.umaryland.edu

"""
Copyright 2025

Licensed under the MIT License, (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at:

    https://mit-license.org/

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

#import required packages
#importing packages to be used with error handling
try:
    import pandas as pd
except:
    print("Required package pandas not available")
    exit()
try:
    import numpy as np
except:
    print("Required package numpy not available")
    exit()
try:
    import argparse
except:
    print("Required package argparse not available")
    exit()
try:
    import sys
except:
    print("Required package sys not available")
    exit()
try:
    import json
except:
    print("Required package json not available")
    exit()

try:
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.Seq import Seq
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except:
    print("Required package Bio not available")
    exit()

try:
    from collections import Counter
except:
    print("Required package collections not available")
    exit()

try:
    import gzip
except:
    print("Required package gzip not available")
    exit()

try:
    import gc
except:
    print("Required package gc not available")
    exit()

########CORE FUNCTIONS
def parseFasta(refdir,fastaFile,ending):
    with open("{refdir}/{fastaFile}.{ending}".format(refdir=refdir,fastaFile=fastaFile,ending=ending)) as fasta_file:
        headers = []
        seqs = []
        for rec, sequence in SimpleFastaParser(fasta_file):
            headers.append(rec)
            seqs.append(str(sequence))
    fastaDF = pd.DataFrame(dict(Genome=headers, seq=seqs))
    return fastaDF

def ComputeRefSig(refdir,fastaFile,ending,k):
    genome = parseFasta(refdir,fastaFile,ending)
    GenomeSeq = Seq(genome['seq'][0])
    revGenomeSeq = GenomeSeq.reverse_complement()
    GenomeSeq = str(GenomeSeq)
    revGenomeSeq = str(revGenomeSeq)
    GenomeSeq = GenomeSeq[-k:] + GenomeSeq 
    revGenomeSeq = revGenomeSeq[-k:] + revGenomeSeq
    forwardKmerList = [GenomeSeq[i:i+k] for i in range(len(GenomeSeq)-(k-1))]
    reverseKmerList = [revGenomeSeq[i:i+k] for i in range(len(revGenomeSeq)-(k-1))]
    KmerList = forwardKmerList+reverseKmerList
    KmerDict = Counter(KmerList)
    return KmerDict

def generateOutput(outputPrefix,defaultOut):
    LBPout=pd.DataFrame.from_dict(defaultOut,orient='index')
    LBPout.columns = [outputPrefix]
    LBPout=LBPout.div(LBPout.sum(axis=0),axis=1)
    LBPout=LBPout.T
    LBPout.to_csv("{out}_strainRelAbund.csv".format(out=outputPrefix))

parser = argparse.ArgumentParser(description="kSanity is a tool to detect and quantify bacterial strains in shotgun metagenomic data")
subparsers = parser.add_subparsers(dest='command')

###Viable commands
#build-ref -> builds the kMer reference from a set of conspecific whole genome sequences
#reads -> builds the kMer database from a fastq file containing sequence reads
#DnQ -> performs strain detection and quantification
#compile -> compiles the relative abundance estimates across all the samples


subparsers.required = True
#  subparser for mapping subroutine
parser_build = subparsers.add_parser('build-ref')

# adding arguments for build subroutine
parser_build.add_argument('-r','--refdir',help='Path to directory containing reference genomes',required=True)
parser_build.add_argument('-m','--manifest',help='File containing list of strains, 1 entry per line',required=True)
parser_build.add_argument('-e','--ending',choices=['fa','fna','fasta'],help='File ending of reference genomes default: fa',default='fa',required=False)
parser_build.add_argument('-k','--kmer',help='kMer length used in analysis, odd numbers recommended between 33-77, default:55',default=55,required=False)
parser_build.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the mapping output',required=True)

#subparser for reads subroutine
parser_build = subparsers.add_parser('reads')

# adding arguments for build subroutine
parser_build.add_argument('-r','--reads',help='Path to the reads file, expected to be a gzipped fastq file',required=True)
parser_build.add_argument('-k','--kmer',help='kMer length used in analysis, odd numbers recommended between 33-77, default:55, must match that used in building the reference',default=55,required=False)
parser_build.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the mapping output',required=True)

#subparser for the DnQ subroutine
parser_build = subparsers.add_parser('DnQ')

# adding arguments for build routine
parser_build.add_argument('-r','--ref',help='Full path to the reference files generated by the build command with prefix but without ending e.g. /path/to/PREFIX where /path/to/ contains PREFIX.',required=True)
parser_build.add_argument('-p','--profile',help='Full path to the sample kMer profile json file produced by the reads command',required=True)
parser_build.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the mapping output',required=True)
parser_build.add_argument('-t','--threshold',help="Percent of unique kMers needed to be identified for a strain to be considered detected, default=0.5",default=0.5,required=False)

#subparser for the DnQ subroutine
parser_build = subparsers.add_parser('compile')

# adding arguments for build routine
parser_build.add_argument('-i','--inputDir',help='Full path to the directory containing the per sample kSanity output files',required=True)
parser_build.add_argument('-o','--outputPrefix',help='Prefix used in the filename for the compiled output',required=True)

args = parser.parse_args()

#
if args.command == 'build-ref':
    try:
        with open(args.manifest) as file:
            strains = [line.rstrip() for line in file]
    except:
        print("Manifest file not read")

    k = int(args.kmer)

    genomeSignatures = {}

    for strain in strains:
        print(strain)
        genomeSignatures[strain]=ComputeRefSig(args.refdir,strain,args.ending,k)

    genomeSignatures_revised = {}
    for strain in strains:
        for kMer in genomeSignatures[strain].keys():
            try:
                genomeSignatures_revised[kMer].append(strain)
            except KeyError as e:
                genomeSignatures_revised[kMer]=[]
                genomeSignatures_revised[kMer].append(strain)

    uniqueKmers = {}
    for strain in strains:

        strainSig = list(genomeSignatures[strain].keys())
        otherSigs = [value for key, value in genomeSignatures.items() if key != strain]
        otherSig = list()

        for x in otherSigs:
            otherSig = otherSig + list(x.keys())
        strainUnique = list(set(strainSig).difference(set(otherSig)))
        
        uniqueKmers[strain] = strainUnique

    with open('{out}.genomeSignatures.json'.format(out=args.outputPrefix), 'w', encoding='utf-8') as aK:
        json.dump(genomeSignatures_revised, aK, ensure_ascii=False, indent=4)

    with open('{out}.uniqueKmers.json'.format(out=args.outputPrefix), 'w', encoding='utf-8') as uK:
        json.dump(uniqueKmers, uK, ensure_ascii=False, indent=4)

if args.command == 'reads':
    
    k = int(args.kmer)
    count=0
    KmerDict = {}
    with gzip.open('{readsFile}'.format(readsFile=args.reads),"rt") as reads:
        for read, seq, qual in FastqGeneralIterator(reads):       
            read_len = len(seq)
            count+=1
            if count % 1000000 == 0:
                print(count)
            if read_len > k:        
                readKmers = [seq[i:i+k] for i in range(read_len-(k-1))]
                for readKmer in readKmers:
                    try:
                        KmerDict[readKmer]+=1
                    except KeyError as e:
                        KmerDict[readKmer]=1

    with open('{}.json'.format(args.outputPrefix), 'w', encoding='utf-8') as rK:
        json.dump(KmerDict, rK, ensure_ascii=False, indent=4)

if args.command== 'DnQ':

    ######## preparing uniqueKmer data
    uniqueKmers = json.load(open("{ref}.uniqueKmers.json".format(ref=args.ref)))
    inverseSignatures = {}
    for key,val in uniqueKmers.items():
        for x in val:
            inverseSignatures[x]=key
    del uniqueKmers
    gc.collect()
    signatureKmers = pd.DataFrame.from_dict(inverseSignatures,orient='index')
    del inverseSignatures
    gc.collect()
    signatureKmers.columns = ['Strain']
    uniquekMerCounts = signatureKmers.groupby(['Strain']).size().to_dict()

    defaultOut={}
    for key, val in uniquekMerCounts.items():
        defaultOut[key]=0

    #######readKmers handling
    readKmers = json.load(open("{}".format(args.profile)))
    readKmers = pd.DataFrame.from_dict(readKmers,orient='index')
    readKmers.columns = ['Count']

    detectReplicons = pd.merge(left=signatureKmers,right=readKmers,left_index=True,right_index=True,how="left")
    del signatureKmers
    gc.collect()

    detectReplicons["DetectedkMer"]=np.where((detectReplicons["Count"]>0),1,0)


    detectReplicons = detectReplicons.groupby('Strain').agg({'Count':'median','DetectedkMer':'mean'})
    detectReplicons=detectReplicons.reset_index()
    detectReplicons["uniqueKmerCount"] = detectReplicons['Strain'].map(uniquekMerCounts)
    
    detectReplicons.columns=['Strain','MedianCov','PercentDetect','NumUniqueKmer']
    detectReplicons=detectReplicons[['Strain','NumUniqueKmer','PercentDetect','MedianCov']]
    #median version of line
    detectReplicons.to_csv("{prefix}_strainDetect.csv".format(prefix=args.outputPrefix),sep=",",index=None)

    detectedReplicons = detectReplicons[detectReplicons['PercentDetect']>=float(args.threshold)]
    listDetectedStrains=list(detectedReplicons['Strain'])

    listDetectedStrains=list(set(listDetectedStrains) & set(defaultOut.keys()))

    if len(listDetectedStrains)==0:
        generateOutput(args.outputPrefix,defaultOut)
        print('No LBP strains detected')
        sys.exit()

    detectedReplicons=detectedReplicons.set_index("Strain")
    detectedReplicons = detectedReplicons.drop(['NumUniqueKmer','PercentDetect'],axis=1)

    allRefKmers = json.load(open("{}.genomeSignatures.json".format(args.ref)))

    coreKmers={}
    for key,val in allRefKmers.items():
        if val == list(defaultOut.keys()):
            coreKmers[key]='Core'

    coreReadKmers= readKmers.loc[readKmers.index.isin(list(coreKmers.keys()))]
    del readKmers
    gc.collect()
    popGenomeCopies=coreReadKmers['Count'].median()

    detectedReplicons['PopPercent'] = detectedReplicons['MedianCov'].div(popGenomeCopies)

    for index,row in detectedReplicons.iterrows():
        defaultOut[index]=row['PopPercent']

    generateOutput(args.outputPrefix,defaultOut)

if args.command== 'compile':

    CompiledDFs = []

    for resultFile in sorted(os.listdir(args.inputDir)):

        if "strainRelAbund" in resultFile:

            DF=pd.read_csv(resultFile,sep=",",index_col=0)
            CompiledDFs.append(DF)

    CompiledDFs=pd.concat(CompiledDFs,axis=0)

    CompiledDFs.to_csv("{}.compiled.relAbund.csv".format(args.outputPrefix),sep=",")