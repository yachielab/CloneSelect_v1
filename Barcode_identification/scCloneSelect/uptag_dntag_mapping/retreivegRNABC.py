"""
A program to retreive upBC based on list of dnBC and barcode whitelist

Run by:
python retreivegRNABCv2.py queryBC.csv scCSlib.csv --distance 

changes from original:
keep only "TOP" Top suggestion with the highest term frequency of the suggestions of smallest edit distance found.

"""
import argparse
from symspellpy import SymSpell, Verbosity
import subprocess
import pandas as pd

def retreiveupBC(queryBCs, scCSlib):
    #read the queryBC list
    dfqueryBC = pd.read_csv(queryBCs)

    #prepare dictionary
    retreiveupBC = {}
    retreiveupBC['Query dnBC'] = dfqueryBC['dnBC'].to_list()
    retreiveupBC['Cell count'] = dfqueryBC['count'].to_list()
    retreiveupBC['Best matching dnBC'] = []
    retreiveupBC['Edit distance'] = []
    retreiveupBC['library read count'] = []
    retreiveupBC['Corresponding upBC'] = []

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

    #match the corresponding dnBC
    scCSlibdf = pd.read_csv(scCSlib)
    for item in retreiveupBC['Best matching dnBC']:
        try:
            upBC = scCSlibdf.loc[scCSlibdf['dnBC'] == item, 'upBC'].iloc[0]
        except IndexError:
            upBC = 'N/A'
        retreiveupBC['Corresponding upBC'].append(upBC)

    retreiveupBCdf = pd.DataFrame(retreiveupBC)

    #if no matching dnBC, change corresponding values to NA
    retreiveupBCdf.loc[retreiveupBCdf['Corresponding upBC'] == 'N/A', 'Best matching dnBC'] = 'N/A'
    retreiveupBCdf.loc[retreiveupBCdf['Corresponding upBC'] == 'N/A', 'Edit distance'] = 'N/A'
    retreiveupBCdf.loc[retreiveupBCdf['Corresponding upBC'] == 'N/A', 'library read count'] = 'N/A'

    retreiveupBCdf.to_csv('retreiveupBC.csv')

    return retreiveupBCdf


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('queryBCfile', help='csv file of querybarcodes (output of cellBC_dnBCretreival)')
    parser.add_argument('scCSlib', help='csv file of scCS whitelist (output of scCSlibrary3)')
    parser.add_argument("--distance", help="max Damerau-Levenshtein distance between query and library BC", type=int, default = 5) #optional
    args = parser.parse_args()

    queryBCs = args.queryBCfile
    scCSlib = args.scCSlib
    maxdistance = args.distance

    print ('\n... retreiving gRNAbarcodes from querydnBCs ...')
    retreiveupBCdf = retreiveupBC(queryBCs,scCSlib)
    queryBCcount = len(retreiveupBCdf.index)
    nonmatchcount = len(retreiveupBCdf.loc[retreiveupBCdf['Best matching dnBC'] == 'N/A'])
    print ('## Number of query barcodes: ', queryBCcount)
    print ('## Percentage of query barcodes with matched barcodes in scCS lib:', round((queryBCcount-nonmatchcount)/queryBCcount*100,2), '%')
    print ('### retreived gRNABC saved as "retreiveupBC.csv" \n')
