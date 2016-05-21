#!/ccnc_bin/mini_env/bin/python
from __future__ import division
__author__ = 'kcho'
import re
import os
import sys
import pandas as pd
#import matplotlib
#matplotlib.use('GTK')
import matplotlib.pyplot as plt
import argparse
import valueSwap
import textwrap
import ccncpy.ccncpy as ccncpy

def openStatsTable(freesurfer_dir):
    statsROI = os.path.join(freesurfer_dir,'stats')
    filesToRead = 'aparc.stats'

    dfList = []
    for side in ['lh','rh']:
        statsFile = os.path.join(statsROI,side+'.'+filesToRead)
        with open(statsFile,'r') as f:
            lines = f.readlines()

        linesEdited = [re.search('(\w+)\s+(\d+)\s+(\d+)\s+(\d+)',x).group(1,4) for x in lines if re.search('^\w',x)]
        df = pd.DataFrame(linesEdited)
        df.columns = ['ROI','Volume']
        df['Volume'] = df['Volume'].astype(int)
        df['side']=side
        dfList.append(df)

    return pd.concat(dfList)


