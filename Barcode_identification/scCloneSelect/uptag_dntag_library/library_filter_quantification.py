"""
proceeds library_extraction_trimming.sh for filtering barcode library based on abundance
or run by python library_filter_quantification.py output of library_extraction_trimming.sh
"""

import argparse
import pandas as pd
from kneed import DataGenerator, KneeLocator
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from Bio.Seq import Seq


#plot readout split barcode table
def BClib (bartenderoutput):
    #open csvfile as dataframe
    updnBCdf = pd.read_csv(bartenderoutput, names=['upBC', 'dnBC', 'count'])
    #remove BC2, BC4 that seemed to be mixed in sample
    updnBCdf = updnBCdf.drop(updnBCdf[(updnBCdf['upBC'] == 'AGCGTGTCAGGGTGACC')].index)
    updnBCdf = updnBCdf.drop(updnBCdf[(updnBCdf['upBC'] == 'AGTCTGTCTCTCACAGC')].index)
    updnBCdf = updnBCdf.drop(updnBCdf[(updnBCdf['dnBC'] == 'TATGAAGCGGGCTCCCTGTTCGGCTTGCCT')].index)
    updnBCdf = updnBCdf.drop(updnBCdf[(updnBCdf['dnBC'] == 'TGTTTGGGTCTCCGTTTGTCATGTTGTGGC')].index)
    #sort inorder of count, reset index
    updnBCdf = updnBCdf.sort_values('count', ascending=False).reset_index(drop=True)
    updnBCdf.index = updnBCdf.index + 1
    #save as csv
    updnBCdf.to_csv('scCSlibv3.csv')
    return updnBCdf

#plot cumulative frequency
def kneebclib ():
    #BCdf = pd.read_csv('BClibsplit.csv' , names=['upBC', 'dnBC', 'count'])
    BCdf = BClib ('BClibsplit.csv')
    knBCdf = BCdf.drop(BCdf.columns[[0]], axis=1).sort_values('count', ascending=False).reset_index(drop=True)
    knBCdf['Cumulative Frequency'] = knBCdf['count'].cumsum().reset_index(drop=True)
    knBCdf.index = knBCdf.index + 1
    index = np.array(knBCdf.index + 1)
    cumfrq = np.array(knBCdf['Cumulative Frequency'])
    cumfrqnor = np.array(knBCdf['Cumulative Frequency']/knBCdf['count'].sum())
    
    xaxis = np.insert(index,0, 0)
    yaxis = np.insert(cumfrqnor,0, 0)
    yaxis2 = np.insert(cumfrq,0, 0)
    
    kneedle = KneeLocator(xaxis, yaxis, S=1, curve='concave', direction='increasing')
    kneepoint = round(kneedle.knee, 3)
                      
    plt.rcParams.update({'figure.autolayout': True})
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax.plot(xaxis,yaxis, color='#F7792B', linewidth=2.0)
    ax2.plot(xaxis, yaxis2, color='#F7792B', linewidth=2.0)
 
    print (' === number of BC based on knee point : ', kneepoint, '\n')
    ax.set_xlabel('barcode number sorted by number of reads')  # Add an x-label to the axes.
    ax.set_ylabel('cumulative frequency of reads')  # Add a y-label to the axes.
    ax2.set_ylabel('cumulative number of reads')
    ax.vlines(x=578, ymin=0, ymax=1,linestyle="dashed", color = 'grey')

    plt.savefig('cumulativekneee.png', format='png', dpi = 200)
    knBClibdf = knBCdf.loc[knBCdf.index[0:kneepoint]].drop(columns ='Cumulative Frequency')
    return knBCdf


## select only most abundant barcode
#since the number of barcodes is not so high, instead of filtering based on knee point for the whitelist,
#filtered for unique barcodes pairs based on abundance.

def refineBClib (updnBCdf):
    #choose best combination per upBC
    upBCdf = updnBCdf.loc[updnBCdf.groupby("upBC")["count"].idxmax()]
    #choose best combination per dnBC
    BCdf = upBCdf.loc[upBCdf.groupby("dnBC")["count"].idxmax()]
    #sort inorder of count, reset index
    BCdf = BCdf.sort_values('count', ascending=False).reset_index(drop=True)
    BCdf.index = BCdf.index + 1
    #save as csv
    BCdf.to_csv('scCSlibv3.2.csv')
    
    return BCdf

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bartenderoutput', help='csv file, output of cellBC_dnBCretreival)')
    args = parser.parse_args()

    splitBCs = args.bartenderoutput
    BCdf = BClib (splitBCs)
    refineBClib (BCdf)
    print ('### library saved as "scCSlibv3.2.csv"')