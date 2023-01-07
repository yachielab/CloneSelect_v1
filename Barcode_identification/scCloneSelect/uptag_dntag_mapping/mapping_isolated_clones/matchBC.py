"""
To match dnBC of sorted samples based on list of barcode whitelist
05.12.22

Run for each sample by:
python matchBC.py dnBCreadfile libraryreadfile --output --distance 

Directory: /home/ha5666/DryLab/Projects/scCloneSelect/scCSExp17.3-3/CloneIsolation/barcodematch
Test run by: python matchBC.py ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/CloneIsolation/barcodeextraction/67-67/BClibsplit.csv ./scCSlibv3.csv --o matchedBC.csv --distance 3
"""
import argparse
from symspellpy import SymSpell, Verbosity
import subprocess
import pandas as pd

def matchBCs(queryBCs, scCSlib,output):
    #read the queryBC list
    dfqueryBC = pd.read_csv(queryBCs, header = None, names = ['upBC','dnBC','count'])

    #prepare dictionary
    retreiveupBC = {}
    retreiveupBC['Query dnBC'] = dfqueryBC['dnBC'].to_list()
    retreiveupBC['Read count'] = dfqueryBC['count'].to_list()
    retreiveupBC['Best matching dnBC'] = []
    retreiveupBC['Edit distance'] = []
    retreiveupBC['library read count'] = []

    #build up whitelist upBC dict for symspell
    p1 = subprocess.Popen(["awk -F, 'FNR > 1 {print $3, $4}' %s > frequency_dictionary.txt" %scCSlib], shell = True)
    p1.wait()

    sym_spell = SymSpell(max_dictionary_edit_distance = maxdistance, prefix_length=10)
    sym_spell.load_dictionary('frequency_dictionary.txt', 0, 1)

    ## find matching upBC within queryBC
    input_terms = retreiveupBC['Query dnBC']
    suggest_list = [sym_spell.lookup(i, Verbosity.TOP, include_unknown=True) for i in input_terms]
    matchlist = []

    for suggestions in suggest_list:
        for suggestion in suggestions:
            matchlist.append (str (suggestion).split(','))

    #input information for matching BC into retreiveupBCdict
    retreiveupBC['Best matching dnBC'] = [item[0] for item in matchlist]
    retreiveupBC['Edit distance'] = [item[1] for item in matchlist]
    retreiveupBC['library read count'] = [item[2] for item in matchlist]

    #Processing for plotting 
    retreiveupBCdf = pd.DataFrame(retreiveupBC)
    retreiveupBCdf['library read count'] = pd.to_numeric(retreiveupBCdf['library read count'])
    #remove reads that are not in the whitelist
    retreiveupBCdf = retreiveupBCdf[retreiveupBCdf['library read count'] != 0]
    #compute frequency of each matching dnBCs
    retreiveupBCdf['Read count freq'] = retreiveupBCdf['Read count']/retreiveupBCdf['Read count'].sum()

    #aggregate reads matching to same dnnBC
    retreiveupBCdf = retreiveupBCdf.groupby(['Best matching dnBC'],as_index = False).agg({'Read count freq': 'sum'})
    retreiveupBCdf = pd.DataFrame(retreiveupBCdf)

    #if no matching dnBC, change corresponding values to NA
    retreiveupBCdf.to_csv(output)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('queryBCfile', help='csv file of querybarcodes (output of cellBC_dnBCretreival)')
    parser.add_argument('scCSlib', help='csv file of scCS whitelist (output of scCSlibrary3)')
    parser.add_argument("--distance", help="max Damerau-Levenshtein distance between query and library BC", type=int, default = 5) #optional
    parser.add_argument("--output", help="name of csv file", default = 'matchedBC.csv') #optional
    args = parser.parse_args()

    queryBCs = args.queryBCfile
    scCSlib = args.scCSlib
    maxdistance = args.distance
    output = args.output

    print ('\n... matching reads with whitelistBCs ...')
    matchBCs(queryBCs, scCSlib,output)
    print ('### results saved as %s \n'% output)
