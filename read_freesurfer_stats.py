#!/ccnc_bin/mini_env/bin/python
from __future__ import division
__author__ = 'kcho'
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import textwrap
plt.style.use('ggplot')

def openStatsTable(freesurfer_dir):
    '''
    Reads freesurfer aparc.stats and
    returns pandas dataframe of all variables
    '''
    statsROI = os.path.join(freesurfer_dir,'stats')
    filesToRead = 'aparc.stats'

    dfList = []
    for side in ['lh','rh']:
        statsFile = os.path.join(statsROI,side+'.'+filesToRead)
        with open(statsFile,'r') as f:
            lines = f.readlines()

        ROIlines = [x.strip() for x in lines if re.search('^\w',x)]

        re_roi_lines= [re.search('(?P<ROI>\w+)\s+(?P<numvert>\S+)\s+(?P<surfarea>\S+)\s+(?P<grayvol>\S+)\s+(?P<thickavg>\S+)\s+(?P<thickstd>\S+)\s+(?P<meancurv>\S+)\s+(?P<gauscurv>\S+)\s+(?P<foldind>\S+)\s+(?P<curvind>\S+)', x) for x in ROIlines]
        df = pd.DataFrame([x.groups() for x in re_roi_lines])

        df.columns = ['roi','numvert','surfarea','grayvol','thickavg','thickstd','meancurv','gauscurv','foldind','curvind']
        df['side']=side

        for i in ['numvert','surfarea','grayvol','thickavg','thickstd','meancurv','gauscurv','foldind','curvind']:
            df[i] = df[i].astype(float)

        dfList.append(df)

    return pd.concat(dfList)

def openStatsTable_big(freesurfer_dir):
    '''
    Using openStatsTable, it first reads freesurfer aparc.stats.
    Then it averages each index within the bigger regions,
    namely
    LPFC, OFC, MPFC, LTC, MTC, SMC, PC, OCC
    '''
    df = openStatsTable(freesurfer_dir)
    roi_region_map = get_roi_region_map()
    df['region'] = df['roi'].map(roi_region_map)
    df.loc[df['region'].isnull(), 'region'] = 'other'

    region_gb = df.groupby(['region','side'])
    meanDf = region_gb.mean().reset_index()

    gb = meanDf.groupby('region')
    meanDf = pd.concat([gb.get_group('LPFC'),
                        gb.get_group('OFC'),
                        gb.get_group('MPFC'),
                        gb.get_group('LTC'),
                        gb.get_group('MTC'),
                        gb.get_group('SMC'),
                        gb.get_group('PC'),
                        gb.get_group('OCC')])

    return meanDf


def graph_ind(ind_df, index):
    '''
    ind_df : pandas dataframe from openStatsTable_big
    index : string, eg) 'numvert', 'surfarea', 'grayvol', 'thickavg'
    '''

    fig = plt.figure(figsize=(16,8))
    fig.suptitle("{index} in eight regions".format(index = index), fontsize=20)

    lh_g = plt.subplot2grid((2,2),(0, 0), colspan=2)
    rh_g = plt.subplot2grid((2,2),(1, 0), colspan=2)

    gb = ind_df.groupby('side')
    lh_g.plot(gb.get_group('lh')[index], 'r', label='Left')
    rh_g.plot(gb.get_group('rh')[index], 'r', label='Right')
    label = ind_df.region.unique()

    min_y = ind_df[index].min()
    max_y = ind_df[index].max()
    std_y = ind_df[index].std()

    for i, s in zip([lh_g, rh_g], ['Left', 'Right']):
        i.set_xlabel(s, fontsize=16)
        i.set_ylabel(index, fontsize=16)

        i.set_xticks(range(8))
        i.set_xticklabels(label)

        i.set_xlim(-.5, 7.5)
        i.set_ylim(min_y - std_y, max_y + std_y)
        i.legend()
        legend = i.legend(frameon = 1)
        frame = legend.get_frame()
        frame.set_facecolor('white')

    return fig

def getICV(freesurfer_dir):
    '''
    Returns ICV in string
    '''
    aseg_stats_file = os.path.join(freesurfer_dir, 'stats/aseg.stats')
    with open(aseg_stats_file, 'r') as f:
        lines = [x for x in f.readlines() if 'Intracranial' in x]
    ICV = re.search('\d+\.\d+', lines[0]).group(0)
    return ICV

def get_roi_region_map():
    a = get_cortical_rois()
    newDict = {}
    for region, roi in a.iteritems():
        for i in roi:
            newDict[i] = region
    return newDict

def get_cortical_rois():
    '''
    returns more detailed cortical divisions
    '''
    roiDict = {'OFC' : ['parsorbitalis', 'medialorbitofrontal', 'lateralorbitofrontal'],
           'MPFC' : ['caudalanteriorcingulate', 'rostralanteriorcingulate', 'superiorfrontal'],
           'LPFC' : [ 'parstriangularis', 'rostralmiddlefrontal', 'frontalpole', 'parsopercularis'],
           'SMC' : [ 'precentral', 'caudalmiddlefrontal', 'postcentral', 'paracentral'],
           'PC' : ['inferiorparietal', 'supramarginal', 'precuneus', 'posteriorcingulate', 'isthmuscingulate', 'superiorparietal'],
           'MTC' : ['entorhinal', 'parahippocampal', 'fusiform'],
           'LTC' : ['transversetemporal', 'superiortemporal', 'bankssts', 'inferiortemporal', 'middletemporal', 'temporalpole'],
           'OCC' : ['pericalcarine', 'lingual', 'lateraloccipital', 'cuneus']}


    return roiDict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            {codeName} :
            ========================================
            eg) {codeName} --input {in_put} --output {output}
            '''.format(codeName=os.path.basename(__file__),
                       in_put = 'subjectLoc',
                       output = 'outLoc')))

    parser.add_argument(
        '-i', '--freesurferDir',
        help='Subject location',
        default=os.getcwd())

    parser.add_argument(
        '-k', '--index',
        help='information needed, ie) grayvol, thickavg, surfarea',
        default='grayvol')

    parser.add_argument(
        '-d', '--dataframe',
        help='Print dataframe',
        action='store_true',
        default=False)

    parser.add_argument(
        '-g', '--graph',
        help='Draw graph',
        action='store_true',
        default=False)

    parser.add_argument(
        '-v', '--ICV',
        help='print ICV',
        action='store_true',
        default=True)

        #default=[x for x in os.listdir(os.getcwd()) if os.path.isdir(x)])

    args = parser.parse_args()

    if args.ICV:
        print getICV(args.freesurferDir)

    df = openStatsTable_big(args.freesurferDir)

    if args.dataframe:
        print df[['roi','region',args.index]]

    if args.graph:
        fig = graph_ind(df, args.index)
        fig.savefig('fig.png')
    #plt.show()



