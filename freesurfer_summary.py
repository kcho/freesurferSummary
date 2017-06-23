#!/ccnc/anaconda2/bin/python
from __future__ import division
import os
import re
from os.path import join, basename, dirname
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import itertools
import textwrap
from progressbar import ProgressBar
import time
#import tksurferCapture

__author__ = 'kcho'
plt.style.use('ggplot')

def freesurferSummary(args):
    '''
    Output freesurfer summary
    '''

    # User input for subject name (to be on the graph)
    subjectInitials = raw_input('Subject initial :')

    # Collect statistics using freesurfer functions
    #     mri_annotation2label
    #     mri_mergelabels
    #     mris_anatomical_stats --> takes time
    # Returns infoDf
    #     ,subroi,numvert,surfarea,grayvol,thickavg,thickstd,meancurv,gauscurv,foldind,curvind,side
    #     0,lh_bankssts,1803.0,1207.0,2743.0,2.301,0.439,0.107,0.019,13.0,1.5,lh
    infoDf = collectStats(args.fsDir)

    # if secondary dir is given
    if args.fsDir2:
        infoDf2 = collectStats(args.fsDir2)
        subjectInitials_2 = raw_input('Second subject initial :')
    else:
        # Read CCNC healthy control average information
        # meanDf eg
        #     ,subroi,thickavg,side,thickstd,subject
        #    0,lh_bankssts,2.577,lh,0.439,NOR01_XXX
        #meanDfLoc = '/ccnc_bin/meanThickness/detailed_mean_2015_12_28.csv'
        #     ,subroi,      thickavg, side, thickstd, subject,   roi
        #    0,lh_bankssts, 2.577,    lh,   0.439,    NOR01_XXX, bankssts
        meanDfLoc = '/ccnc_bin/meanThickness/detailed_mean_2017_06_16.csv'
        meanDf = pd.read_csv(meanDfLoc, index_col=0)
        if 'roi' not in meanDf.columns or 'region' not in meanDf.columns:
            meanDf['roi'] = meanDf.subroi.str[3:]
            meanDf['region'] = meanDf.roi.apply(getRegion)
            meanDf.to_csv(meanDfLoc)

        # meanDf to have averaged values for all groups
        infoDf2 = meanDf.groupby(['roi','side','region']).mean().reset_index()
        infoDf2.columns = ['roi','side','region', 'thickavg','thickstd']
        subjectInitials_2 = 'CCNC_mean'

    ## Amend information Dfs for previously created csvs
    #for df_tmp in infoDf, meanDf: 
        #print df_tmp
        #if 'roi' not in df_tmp.columns:
            #df_tmp['roi'] = df_tmp.subroi.str[3:]
            #print("The subject {subjName}'s csv requires editing, to have 'roi' column".
                  #format(subjName=subjName))
        #if 'region' not in df_tmp.columns:
            #df_tmp['roi'] = df_tmp.subroi.str[3:]
            #df_tmp['region'] = df_tmp.roi.apply(getRegion)
            #print("The subject {subjName}'s csv requires editing, to have 'region' column".
                  #format(subjName=subjName))

    # Make line plots of cortical thickness for each hemisphere
    draw_thickness(args.fsDir,
                   infoDf, infoDf2,
                   subjectInitials, subjectInitials_2)

    # tksurferCapture.main(args.fsDir, join(args.fsDir,
                                          # 'tmp',
                                          # 'thick_kev_detailed_new.csv'))


def getRegion(roi):
    '''
    find region of the detailed freesurfer label
    eg) regions : 'LPFC', 'OFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC' or 'OCC'
    '''
    roiDict = get_cortical_rois()
    for region, roiList in roiDict.iteritems():
        if roi in roiList:
            return region


def reorder_df(df, colName, orderList):
    gb = df.groupby(colName)
    newDf = pd.concat([gb.get_group(x) for x in orderList])
    newDf = newDf.reset_index()
    return newDf

def draw_thickness(fsDir, infoDf, meanDf, subjName, meanDfName):
    # graph order from left
    roiOrder = ['LPFC', 'OFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC', 'OCC']


    # side group by
    infoDf_gb = infoDf.groupby('side')
    meanDf_gb = meanDf.groupby('side')


    fig, axes = plt.subplots(nrows=2, figsize=(22,12), facecolor='white')
    fig.suptitle("Cortical thickness in all regions", fontsize=20)

    for snum, side in enumerate(['lh', 'rh']):
        infodf = infoDf_gb.get_group(side)
        meandf = meanDf_gb.get_group(side)

        # Reorder Dfs
        infodf = infodf.sort_values(['roi','side'])
        meandf = meandf.sort_values(['roi','side'])
        infodf = reorder_df(infodf, 'region', roiOrder)
        meandf  = reorder_df(meandf, 'region', roiOrder)


        ax = axes[snum]

        ax.plot(infodf.thickavg, 'r', label=subjName)
        ax.plot(meandf.thickavg, 'b--', label='CCNC HCs')

        # errorbar(x, y, std)
        roiList = infodf.roi.unique()
        eb1 = ax.errorbar(range(len(roiList)),
                          meandf.thickavg,
                          meandf.thickstd,
                          marker='^',
                          label='_nolegend_', 
                          color='b',
                          ecolor='b')
        eb1[-1][0].set_linestyle('--')

        # Graph settings
        ax.set_ylim(1.0, 5)
        ax.set_xlabel(side, fontsize=16)
        ax.set_xticks(range(len(roiList)))
        ax.set_xticklabels(['' for x in roiList])
        ax.set_xlim(-.5, 32.5)
        ax.grid(False)

        ax.legend()
        legend = ax.legend(frameon = 1)
        frame = legend.get_frame()
        frame.set_facecolor('white')

        # Background fill (group discrimination)
        roiDict = get_cortical_rois()
        roiOrder_full = [[x]*len(roiDict[x]) for x in roiOrder]
        roiOrder_one_list = list(itertools.chain.from_iterable(roiOrder_full))
        roiOrder_array = np.array(roiOrder_one_list)

        regionToHighlight = roiOrder[0::2]
        regionToNotHighlight = roiOrder[0::2]
        xCoords = [np.where(roiOrder_array==x)[0] for x in regionToHighlight]
        for x in xCoords:
            col='green'
            ax.axvspan(x[0], x[-1], alpha=0.5)

        #boundArr = np.array([])

        #for regionNum, region in enumerate(roiOrder):
            #boundArr = np.append(boundArr, np.tile(regionNum, len(roiDict[region])))

        #print boundArr

        startNum = 0
        switch = 0
        x_starting_point = 0
        for region in roiOrder:
            if switch == 0:
                switch = 1 
                pass
                alpha = 1
                col='white'
                #p = ax.axvspan(x_starting_point-.5, startNum-.5, facecolor=col, alpha=alpha)
                x_starting_point = startNum
                startNum = startNum + len(roiDict[region])

                ax.text((x_starting_point-.5 + startNum-.5)/2, 1.2,
                        region,
                        horizontalalignment='center',
                        alpha=.4,
                        fontsize=15)
            else:
                alpha = 0.2
                col='yellow'
                #p = ax.axvspan(x_starting_point-.5, startNum-.5, facecolor=col, alpha=alpha)
                x_starting_point = startNum
                startNum = startNum + len(roiDict[region])
                switch = 0 

                ax.text((x_starting_point-.5 + startNum-.5)/2, 1.2,
                        region,
                        horizontalalignment='center',
                        alpha=.4,
                        fontsize=15)

        ## annotation
        mergedDf = pd.merge(meandf,
                            infodf,
                            on=['roi','side','region'],
                            how='inner')

        mergedDf['mean_sub_indv'] = mergedDf.thickavg_x - mergedDf.thickavg_y
        for row in mergedDf[abs(mergedDf['mean_sub_indv']) > 0.5].iterrows():
            if row[1].mean_sub_indv < 0:
                col = 'green'
            else:
                col = 'red'

            ax.annotate(row[1].roi,
                    xy=(row[0], row[1].thickavg_y),
                    xytext=(row[0], 1.5-row[1].mean_sub_indv/3),
                    arrowprops=dict(facecolor=col, shrink=0.05),
                    horizontalalignment=side,
                    fontsize=20)

    ax.set_xticklabels(infodf.roi)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=30)
    plt.tight_layout(pad=7, w_pad=3, h_pad=0.2)
    #fig.show()
    fig.savefig(join(fsDir, basename(subjName)+'_thickness.png'))
    print(join(fsDir, basename(subjName)+'_thickness.png'))



def dictWithTuple2df(infoDict):
    df = pd.DataFrame.from_dict(infoDict)
    df = df.T.reset_index()
    df.columns = ['subroi', 'numvert', 'surfarea', 'grayvol',
                  'thickavg', 'thickstd',
                  'meancurv', 'gauscurv', 'foldind', 'curvind']
    df['side'] = df['subroi'].str[:2]
    df['roi'] = df.subroi.str[3:]
    df['region'] = df.roi.apply(getRegion)

    return df



def getInfoFromLabel(fsDir,roiDict):
    '''
    Change this function to use pandas and numpy
    '''

    infoDict={}

    print 'fsDir',fsDir
    print 'roiDict', roiDict

    pbar = ProgressBar().start()
    totalNum = 2 * len(roiDict.keys())
    num = 1
    for side in ['lh','rh']:
        for cortex, rois in roiDict.iteritems():
            if len(rois) > 1:
                command = 'mris_anatomical_stats \
                -l {loc}/{side}_{cortex} {name} {side} 2>/dev/null'.format(
                    loc=join(fsDir,'tmp'),
                    side=side,
                    cortex=cortex,
                    name=basename(fsDir)
                )
            else:
                command = 'mris_anatomical_stats \
                -l {loc}/{side}.{cortex}.label {name} {side} 2>/dev/null'.format(
                    loc=join(fsDir,'tmp'),
                    side=side,
                    cortex=cortex,
                    name=basename(fsDir)
                )
            output=os.popen(re.sub('\s+',' ',command)).read()
            pbar.update((num/totalNum) * 100)
            num+=1

            print output

            thickness = re.search('thickness\s+=\s+(\S+)\s+mm\s+\S+\s+(\S+)', output).group(1,2)
            numvert = re.search('number of vertices\s+=\s+(\S+)', output).group(1)
            surfarea = re.search('total surface area\s+=\s+(\S+)', output).group(1)
            grayvol = re.search('total gray matter volume\s+=\s+(\S+)', output).group(1)
            meancurv = re.search('average integrated rectified mean curvature\s+=\s+(\S+)', output).group(1)
            gauscurv = re.search('average integrated rectified Gaussian curvature\s+=\s+(\S+)', output).group(1)
            foldind = re.search('folding index\s+=\s+(\S+)', output).group(1)
            curvind = re.search('intrinsic curvature index\s+=\s+(\S+)', output).group(1)

            thickness = tuple([float(x) for x in thickness])
            infoDict[side+'_'+cortex] = [float(x) for x in [numvert, surfarea, grayvol,
                                              thickness[0], thickness[1],
                                              meancurv, gauscurv, foldind, curvind]]
    pbar.finish()
    return infoDict


def mergeLabel(fsDir, roiDict):
    '''
    Merge labels into one label
    merge the labels in the roiDict.values --> label  roiDict.keys
    '''

    for side in ['lh','rh']:
        for cortex, rois in roiDict.iteritems():
            inLabelLocs = [join(fsDir,'tmp',side+'.'+x+'.label') for x in rois]
            inLabelForms = ' '.join(['-i '+x for x in inLabelLocs])
            outLabel = join(fsDir,'tmp',side+'_'+cortex)

            command = 'mri_mergelabels \
                    {inLabel} \
                    -o {outLabel} \
                    2>/dev/null'.format(inLabel = inLabelForms,
                                        outLabel = outLabel)
            os.popen(command).read()


def makeLabel(fsDir):
    '''
    Run mri_annotation2label for lh and rh hemisphere.
    Creates labels in $fsDir/tmp
    '''

    fsDirName = basename(fsDir)
    labelOutDir = join(fsDir, 'tmp')

    # below must be defined here, in order for the FS to run
    os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
    os.environ["SUBJECTS_DIR"] = dirname(fsDir)

    for side in ['lh','rh']:
        command = 'mri_annotation2label \
            --subject {basename} \
            --hemi {side} \
            --outdir {outDir} \
            --ctab {outDir}/{side}_ctab.txt \
                2>/dev/null'.format(basename=fsDirName, 
                                    side=side, 
                                    outDir=labelOutDir)

        print re.sub('\s+',' ',command)
        os.popen(re.sub('\s+',' ',command)).read()

def collectStats(fsDir):
    '''
    CollectStats
    Summarise cortical thickness in more than one subjects.
    'background_subject_locs' should given in python list format.
    For each background,
    1. Creates labels using makeLabel
    2. Merges labels according to the roiDict using mergeLabel
    3. Estimates cortical thickness in each merged labels (dict)
    4. Converts dict to pandas Dataframe using dictWithTuple2df
    5. save df to fsDir/tmp/thick_kev_detailed.csv
    Than mean df is returned.
    '''

    # cortical regions as a dictionary
    roiDict = get_cortical_rois_detailed()

    if not os.path.isfile(join(fsDir,'tmp','thick_kev_detailed_new.csv')):
        makeLabel(fsDir)
        mergeLabel(fsDir, roiDict)
        infoDict = getInfoFromLabel(fsDir, roiDict)
        infoDf = dictWithTuple2df(infoDict)
        infoDf.to_csv(join(fsDir,'tmp','thick_kev_detailed_new.csv'))
    else:
        infoDf = pd.read_csv(join(fsDir,'tmp','thick_kev_detailed_new.csv'),
                             index_col=0)
        if 'roi' not in infoDf.columns or 'region' not in infoDf.columns:
            infoDf['roi'] = infoDf.subroi.str[3:]
            infoDf['region'] = infoDf.roi.apply(getRegion)
            infoDf.to_csv(join(fsDir,'tmp','thick_kev_detailed_new.csv'))

    # return mean
    return infoDf

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

def getSummary(volumeDf,roiDict):
    # thalamus : lh, 10, rh, 49
    # lh_OFC : 1019 1014 1012
    columnMake = pd.DataFrame.from_dict(roiDict,orient='index').T.stack().reset_index()
    columnMake.columns = ['order','cortex','subroi']

    volumeDf =  pd.merge(volumeDf,
             columnMake[['cortex','subroi']],
             left_on='ROI',
             right_on='subroi',
             how='inner').drop('subroi',axis=1)

    return volumeDf

def openStatsTable(fsDir):
    statsROI = join(fsDir,'stats')
    filesToRead = 'aparc.stats'

    dfList = []
    for side in ['lh','rh']:
        statsFile = join(statsROI,side+'.'+filesToRead)
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


def concatFsDf(freesurferDirList):
    dfList = []

    for freesurferDir in freesurferDirList:
        dfList.append(collectStats_v2(freesurferDir))
    dfMerged = pd.concat(dfList)

    return dfMerged

def concatDf_to_meanDf(concatDf):
    meanDf = concatDf.groupby(concatDf.index).mean()
    infoCols = concatDf.ix[concatDf.index, ['subroi','side']].drop_duplicates()
    print infoCols
    return pd.concat([infoCols, meanDf], axis=1)

def subjDirs_to_fsDirs(subjectDirList):
    freesurferDirList = []

    for subjDir in subjectDirList:
        freesurferDirList.append(get_freesurferDir(subjDir))

    return freesurferDirList

def get_cortical_rois():
    '''
    returns 8 cortical divisions in dictionary
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

def get_cortical_rois_detailed():
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

    detailed = []
    for key, roi_list in roiDict.iteritems():
        detailed = detailed + roi_list

    detailed_ROIs = {}
    for roi in detailed:
        detailed_ROIs[roi] = [roi]

    return detailed_ROIs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=textwrap.dedent('''\
            {codeName} : 
            ========================================
            eg) {codeName} --input {in_put} --output {output}
            '''.format(codeName=basename(__file__),
                       in_put = 'subjectLoc',
                       output = 'outLoc')))

    parser.add_argument(
        '-d', '--fsDir',
        help='Freesurfer directory location',
        default=os.getcwd())

    parser.add_argument(
        '-d2', '--fsDir2',
        help='Secondary freesurfer directory location',
        default=False)

    args = parser.parse_args()

    os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
    os.environ["SUBJECTS_DIR"] = dirname(dirname(args.fsDir))

    freesurferSummary(args)
