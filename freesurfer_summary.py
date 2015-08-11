__author__ = 'kcho'

import re
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main(subject_loc = '/Users/kcho/T1', locations=['/Users/kcho/T1','/Users/kcho/T1'], roi_list = ['ctx_lh_G_cuneus']):
    freesurfer_dir = get_freesurfer_loc(subject_loc)
    print freesurfer_dir


    os.environ["FREESURFER_HOME"] = '/Applications/freesurfer'
    os.environ["SUBJECTS_DIR"] = '{0}'.format(subject_loc)

    #freesurfer_table = openTable('/Applications/freesurfer/FreeSurferColorLUT.txt')
    #print freesurfer_table

    #for roi in roi_list:
    #    print freesurfer_table[roi]

    volumeDf = openStatsTable(freesurfer_dir)
    volumeDf['name'] = volumeDf.side + '_' + volumeDf.ROI

    volumeDf = getSummary(volumeDf)
    print volumeDf


    # graph
    draw_graph(volumeDf)

    # collect stats
    if len(locations) > 1:
        meanDf = collectStats(locations)


def collectStats(locations):
    subjectDict = {}
    for location in locations:
        print location
        freesurfer_dir = get_freesurfer_loc(location)
        volumeDf = openStatsTable(freesurfer_dir)
        volumeDf['name'] = volumeDf.side + '_' + volumeDf.ROI
        volumeDf = getSummary(volumeDf)
        subjectDict[location] = volumeDf

    print len(subjectDict.values())
    #print pd.concat(subjectDict.values())


def draw_graph(volumeDf):
    gb = volumeDf.groupby(['side','cortex'])
    lh_volume_sums={}
    rh_volume_sums={}
    cortexList = []
    for (side,cortex), grp in gb:
        print cortex
        if side == 'lh':
            lh_volume_sums[side+'_'+cortex] = grp['Volume'].sum()
            cortexList.append(cortex)
        else:
            rh_volume_sums[side+'_'+cortex] = grp['Volume'].sum()

    plt.plot(lh_volume_sums.values(),'r')
    plt.plot(rh_volume_sums.values(),'b')
    plt.xticks(range(len(cortexList)), cortexList)
    plt.show()
    # plt.plot(gb.get_group('lh')['Volume'],'r-')
    # plt.plot(gb.get_group('rh')['Volume'],'b-')
    # plt.xticks(range(len(gb.get_group('lh')['ROI'])),gb.get_group('lh')['ROI'])
    # plt.show()

def getSummary(volumeDf):
    # thalamus : lh, 10, rh, 49
    # lh_OFC : 1019 1014 1012
    roiDict = {'OFC' : ['parsorbitalis', 'medialorbitofrontal', 'lateralorbitofrontal'],
               'MPFC' : ['caudalanteriorcingulate', 'rostralanteriorcingulate', 'superiorfrontal'],
               'LPFC' : [ 'parstriangularis', 'rostralmiddlefrontal', 'frontalpole', 'parsopercularis'],
               'SMC' : [ 'precentral', 'caudalmiddlefrontal', 'postcentral', 'paracentral'],
               'PC' : ['inferiorparietal', 'supramarginal', 'precuneus', 'posteriorcingulate', 'isthmuscingulate', 'superiorparietal'],
               'MTC' : ['entorhinal', 'parahippocampal', 'fusiform'],
               'LTC' : ['transversetemporal', 'superiortemporal', 'bankssts', 'inferiortemporal', 'middletemporal', 'temporalpole'],
               'OCC' : ['pericalcarine', 'lingual', 'lateraloccipital', 'cuneus']}

    columnMake = pd.DataFrame.from_dict(roiDict,orient='index').T.stack().reset_index()
    columnMake.columns = ['order','cortex','subroi']

    volumeDf =  pd.merge(volumeDf,
             columnMake[['cortex','subroi']],
             left_on='ROI',
             right_on='subroi',
             how='inner').drop('subroi',axis=1)

    return volumeDf



#    ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
#    ts = ts.cumsum()
#    plt.figure(); ts.plot();plt.show()


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


def openTable(f_file):
    with open(f_file,'r') as f:
        lines = f.readlines()
    lines_edited = [re.search('^(\d+)\s+(\S+)',x).group(2,1) for x in lines if re.search('^\d',x)]
    lines_edited_dict = dict(lines_edited)

    return lines_edited_dict

def get_freesurfer_loc(location):
    freesurferDIR = re.search('freesurfer',
            ' '.join(os.listdir(location)),
            re.IGNORECASE).group(0)
    return os.path.join(location, freesurferDIR)

def roi_extraction(subjectDir, roiName, roiNumber=False, outputDir=False):
    '''
    Extracts ROI from freesurfer output
    roiName : string
    roiNumber : list of number
    '''

    # FREESURFER vs freesurfer
    inputName =  '{subjectDir}/{freesurferDir}/mri/aparc+aseg.mgz'.format(
                subjectDir=subjectDir,
                freesurferDir=freesurferDIR)

    # If the output directory is specified
    if outputDir:
        outputName = '{outputDir}/{roiName}.nii.gz'.format(
                outputDir = outputDir,
                roiName = roiName)
    else:
        outputName = '{subjectDir}/ROI/{roiName}.nii.gz'.format(
                subjectDir=subjectDir,
                roiName=roiName)


    # If the number of the ROI is secificed
    if roiNumber:
        command = 'mri_binarize --i {inputName} \
                                --match {roiNumber} \
                                --o {outputName}'.format(
                                    inputName=inputName,
                                    roiNumber=roiNumber,
                                    outputName=outputName)

    else:


        command = 'mri_binarize --i {inputName} \
                                --match {roiNumber} \
                                --o {outputName}'.format(
                                    inputName=inputName,
                                    roiNumber=roiNumber,
                                    outputName=outputName)

    command = re.sub('\s+',' ',command)
    output = os.popen(command).read()
    print output



if __name__ == '__main__':
    main()



