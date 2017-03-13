#!/ccnc_bin/mini_env/bin/python
from __future__ import division

__author__ = 'kcho'
import os
import re
from os.path import join, basename, dirname
import sys
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import textwrap
from progressbar import ProgressBar
import time

plt.style.use('ggplot')

def freesurferSummary(args):
    infoDf = collectStats(args.fsDir)

    meanDfLoc = '/ccnc_bin/meanThickness/detailed_mean_2015_12_28.csv'
    meanDf = pd.read_csv(meanDfLoc, index_col=0)

    subjectInitials = raw_input('Subject initial :')

    # Graph
    draw_thickness_detailed(args.fsDir,
                            infoDf,
                            meanDf,
                            subjectInitials,
                            'CCNC_mean')

    # valueSwap.main(main_fsDir, 
                   # join(main_fsDir,
                    # 'tmp/thick_kev_detailed.csv'))


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

def draw_thickness_detailed(fsDir, infoDf, meanDf, subjName, meanDfName):
    # graph order from left
    roiOrder = ['LPFC', 'OFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC', 'OCC']

    # Amend information Dfs
    infoDf['roi'] = infoDf.subroi.str[3:]
    infoDf['region'] = infoDf.roi.apply(getRegion)
    meanDf['roi'] = meanDf.subroi.str[3:]
    meanDf = meanDf.groupby(['roi','side']).mean().reset_index()
    meanDf['region'] = meanDf.roi.apply(getRegion)
    meanDf.columns = ['roi','side','thickavg','thickstd','region']

    # Reorder Dfs
    infoDf = reorder_df(infoDf, 'region', roiOrder)
    meanDf  = reorder_df(meanDf, 'region', roiOrder)

    infoDf_side_gb = infoDf.groupby('side')
    label = infoDf.subroi.str[3:].unique()

    fig, (lh_g, rh_g) = plt.subplots(nrows=2, figsize=(22,12))
    fig.suptitle("Cortical thickness in all regions", fontsize=20)

    # lh_g = plt.subplot2grid((2,2),(0, 0), colspan=2)
    # rh_g = plt.subplot2grid((2,2),(1, 0), colspan=2)

    lh_g.plot(infoDf.groupby('side').get_group('lh').index,
             infoDf.groupby('side').get_group('lh').thickavg,
            label=subjName)
    # lh_g.plot(infoDf.groupby('side').get_group('lh')['thickavg'],'r',label=subjName)
    # lh_g.plot(meanDf.groupby('side').get_group('lh')['thickavg'],'r--',label=meanDfName)

    # error bar

    # eb1 = lh_g.errorbar(range(len(meanDf.roi.unique())),
                        # meanDf.groupby('side').get_group('lh')['thickavg'],
                        # meanDf.groupby('side').get_group('lh')['thickstd'],
                        # linestyle='None',
                        # marker='^')
    # eb1[-1][0].set_linestyle('--')

    lh_g.set_xlabel('Left', fontsize=16)

    lh_g.set_ylim(1.0, 5)

    lh_g.set_xticks(range(len(label)))
    lh_g.set_xticklabels(['' for x in label])
    lh_g.set_xlim(-.5, 32.5)
    lh_g.legend()
    legend = lh_g.legend(frameon = 1)
    frame = legend.get_frame()
    frame.set_facecolor('white')

    ######fill
    roiDict = get_cortical_rois()
    startNum = 0
    switch = 0
    x_starting_point = 0
    for region in roiOrder:
        if switch == 0:
            switch = 1 
            pass
            x_starting_point = startNum
            startNum = startNum + len(roiDict[region])

            lh_g.text((x_starting_point-.5 + startNum-.5)/2, 1.2,
                    region,
                    horizontalalignment='center',
                    alpha=.4,
                    fontsize=15)
        else:
            alpha = 0.2
            col='green'
            p = lh_g.axvspan(x_starting_point-.5, startNum-.5, facecolor=col, alpha=alpha)
            x_starting_point = startNum
            startNum = startNum + len(roiDict[region])
            switch = 0 

            lh_g.text((x_starting_point-.5 + startNum-.5)/2, 1.2,
                    region,
                    horizontalalignment='center',
                    alpha=.4,
                    fontsize=15)

    ## annotation
    mergedDf = pd.merge(meanDf.groupby('side').get_group('lh'),
                        infoDf_side_gb.get_group('lh'),
                        on=['roi','side','region'],
                        how='inner')

    mergedDf['mean_sub_indv'] = mergedDf.thickavg_x - mergedDf.thickavg_y
    for row in mergedDf[abs(mergedDf['mean_sub_indv']) > 0.5].iterrows():
        lh_g.annotate(row[1].roi,
                xy=(row[0], row[1].thickavg_y),
                xytext=(row[0], 1.5-row[1].mean_sub_indv/3),
                arrowprops=dict(facecolor='green', shrink=0.05),
                horizontalalignment='left',
                fontsize=20)


    rh_g.plot(infoDf_side_gb.get_group('rh')['thickavg'],'b',label=subjName)
    rh_g.plot(meanDf.groupby('side').get_group('rh')['thickavg'],'b--',label=meanDfName)

    #error bar
    eb2 = rh_g.errorbar(range(len(meanDf.roi.unique())),
                        meanDf.groupby('side').get_group('rh')['thickavg'],
                        meanDf.groupby('side').get_group('rh')['thickstd'],
                        linestyle='None',
                        marker='^',
                        color='b')
    eb2[-1][0].set_linestyle('--')

    #label = ['LPFC' 'OFC' 'MPFC' 'LTC' 'MTC' 'SMC' 'PC' 'OCC','b']
    xticksNum = range(8)

    #rh_g.set_xticklabels(label)
    rh_g.set_xlabel('Right', fontsize=16)
    rh_g.set_ylim(1, 5)
    rh_g.set_xlim(-.5, 32.5)
    rh_g.set_xticks(range(len(label)))
    rh_g.set_xticklabels(label)


    rh_g.legend()
    legend = rh_g.legend(frameon = 1)
    frame = legend.get_frame()
    frame.set_facecolor('white')

    ######fill
    startNum = 0
    switch = 0
    x_starting_point = 0
    for region in roiOrder:
        if switch == 0:
            switch = 1 
            pass
            x_starting_point = startNum
            startNum = startNum + len(roiDict[region])

            rh_g.text((x_starting_point-.5 + startNum-.5)/2, 1.2,
                    region,
                    horizontalalignment='center',
                    alpha=.4,
                    fontsize=15)
        else:
            alpha = 0.2
            col='green'
            p = rh_g.axvspan(x_starting_point-.5, startNum-.5, facecolor=col, alpha=alpha)
            x_starting_point = startNum
            startNum = startNum + len(roiDict[region])
            switch = 0 

            rh_g.text((x_starting_point-.5 + startNum-.5)/2, 1.2,
                    region,
                    horizontalalignment='center',
                    alpha=.4,
                    fontsize=15)

    ## annotation
    mergedDf = pd.merge(meanDf.groupby('side').get_group('rh'),
                        infoDf_side_gb.get_group('rh'),
                        on=['roi','side','region'],
                        how='inner')

    mergedDf['mean_sub_indv'] = mergedDf.thickavg_x - mergedDf.thickavg_y
    for row in mergedDf[abs(mergedDf['mean_sub_indv']) > 0.5].iterrows():
        rh_g.annotate(row[1].roi,
                xy=(row[0], row[1].thickavg_y),
                xytext=(row[0], 1.5-row[1].mean_sub_indv/3),
                arrowprops=dict(facecolor='green', shrink=0.05),
                horizontalalignment='left',
                fontsize=20)

    plt.tight_layout(pad=7, w_pad=3, h_pad=0.2)


    labels = rh_g.get_xticklabels()
    plt.setp(labels, rotation=30)
    labels = lh_g.get_xticklabels()
    plt.setp(labels, rotation=30)
    fig.savefig(join(fsDir, basename(subjName)+'_thickness'))


# def main(subject_loc, background_subject_locs, graph, meanDfLoc,verbose, brain):


    # ##########################################################
    # # Freesurfer setting
    # ##########################################################

    # ##########################################################
    # # Summarize cortical thickness
    # ##########################################################
    # # if no mean table is given but background list is given
    # if meanDfLoc == False and background_subject_locs != None:
        # print 'Calculating cortical thickness in {0}'.format(background_subject_locs)
        # # make mean table
        # meanDfName = raw_input('Name of the background subject : ')

        # if verbose:
            # meanDfList = []
            # for subjectDir in background_subject_locs:
                # meanDfList.append(collectStats_v2(subjectDir))
            # meanDf = pd.concat(meanDfList)
        # else:
            # meanDf = collectStats(background_subject_locs)

    # # if no mean table is given, and backround list is empty
    # # use subject_loc as the background
    # elif meanDfLoc == False and background_subject_locs == None:
        # print 'No background subjects are given, now running only {0}'.format(subject_loc)
        # meanDfName = ''

        # if verbose:
            # meanDf = collectStats_v2([subject_loc])
        # else:
            # meanDf = collectStats([subject_loc])

    # # if meanDfLoc is given
    # else:
        # print 'Now comparing with mean_thickness.csv in /ccnc_bin/meanThickness/'
        # meanDfName = 'NOR'
        # #meanDf = pd.read_csv('/ccnc_bin/meanThickness/new_mean_OCT15.csv')
        # if verbose:
            # meanDf = pd.read_csv('/ccnc_bin/meanThickness/detailed_mean_2015_12_28.csv', index_col=0)
        # else:
            # meanDf = pd.read_csv('/ccnc_bin/meanThickness/new_mean_OCT15.csv', index_col=0)

    # ##########################################################
    # # Get roi dict : 8 cortex each hemisphere
    # ##########################################################
    # roiDict = get_cortical_rois()

    # ##########################################################
    # # annotation2label --> merge labels --> freesurfer/tmp
    # ##########################################################
    
    # if graph:
        # if verbose:
            # thicknessDf = collectStats_v2([dirname(fsDir)])#background_subject_locs)
            # draw_thickness_detailed(thicknessDf,meanDf,basename(subject_loc), meanDfName, subject_loc)
        # else:
            # thicknessDf = collectStats([dirname(fsDir)])#background_subject_locs)
            # draw_thickness(thicknessDf,meanDf,basename(subject_loc), meanDfName, subject_loc)


    # if brain:
        # if verbose:
            # valueSwap.main(fsDir, 
                    # join(fsDir,
                        # 'tmp/thick_kev_detailed.csv'))
        # else:
            # print 'Please select verbose mode'
            # pass

    # print "*"*80
    # print 'Outputs are saved in /ccnc/mri_team'
    # print "*"*80



def draw_thickness(infoDf,meanDf, subjName, meanDfName, subject_loc):
    infoDf['roi'] = infoDf.subroi.str[3:]
    meanDf['roi'] = meanDf.subroi.str[3:]

    gb = infoDf.groupby('roi')
    infoDf = pd.concat([gb.get_group('LPFC'),
                        gb.get_group('OFC'),
                        gb.get_group('MPFC'),
                        gb.get_group('LTC'),
                        gb.get_group('MTC'),
                        gb.get_group('SMC'),
                        gb.get_group('PC'),
                        gb.get_group('OCC')])


    gbmean = meanDf.groupby('roi')
    meanDf = pd.concat([gbmean.get_group('LPFC'),
                        gbmean.get_group('OFC'),
                        gbmean.get_group('MPFC'),
                        gbmean.get_group('LTC'),
                        gbmean.get_group('MTC'),
                        gbmean.get_group('SMC'),
                        gbmean.get_group('PC'),
                        gbmean.get_group('OCC')])
    print meanDf
    gb = infoDf.groupby('side')
    label = infoDf.subroi.str[3:].unique()

    fig = plt.figure(figsize=(12,8))
    fig.suptitle("Cortical thickness in eight regions", fontsize=20)

    lh_g = plt.subplot2grid((2,2),(0, 0), rowspan=2)
    rh_g = plt.subplot2grid((2,2),(0, 1), rowspan=2)
    lh_g.plot(gb.get_group('lh')['thickavg'],'r',label=subjName)

    lh_g.plot(meanDf.groupby('side').get_group('lh')['thickavg'],'r--',label=meanDfName)

    # error bar

    eb1 = lh_g.errorbar(range(len(meanDf.roi.unique())),
                        meanDf.groupby('side').get_group('lh')['thickavg'],
                        meanDf.groupby('side').get_group('lh')['std'],
                        linestyle='None',
                        marker='^')
    eb1[-1][0].set_linestyle('--')

    lh_g.set_xlabel('Left', fontsize=16)
    lh_g.set_ylabel('Cortical thickness in mm', fontsize=16)
    lh_g.set_ylim(1.0, 4)

    lh_g.set_xticks(range(8))
    lh_g.set_xticklabels(label)
    lh_g.set_xlim(-.5, 7.5)
    lh_g.legend()
    legend = lh_g.legend(frameon = 1)
    frame = legend.get_frame()
    frame.set_facecolor('white')




    rh_g.plot(gb.get_group('rh')['thickavg'],'b',label=subjName)
    rh_g.plot(meanDf.groupby('side').get_group('rh')['thickavg'],'b--',label=meanDfName)

    # error bar
    eb2 = rh_g.errorbar(range(len(meanDf.roi.unique())),
                        meanDf.groupby('side').get_group('rh')['thickavg'],
                        meanDf.groupby('side').get_group('rh')['std'],
                        linestyle='None',
                        marker='^',
                        color='b')
    eb2[-1][0].set_linestyle('--')

    #label = ['LPFC' 'OFC' 'MPFC' 'LTC' 'MTC' 'SMC' 'PC' 'OCC','b']
    xticksNum = range(8)

    #rh_g.set_xticklabels(label)
    rh_g.set_xlabel('Right', fontsize=16)
    rh_g.set_ylim(1, 4)
    rh_g.set_xlim(-.5, 7.5)
    rh_g.set_xticks(range(8))
    rh_g.set_xticklabels(label)


    rh_g.legend()
    legend = rh_g.legend(frameon = 1)
    frame = legend.get_frame()
    frame.set_facecolor('white')

    plt.savefig('/ccnc/mri_team/'+subject_loc)


def dictWithTuple2df(infoDict):
    df = pd.DataFrame.from_dict(infoDict)
    df = df.T.reset_index()
    df.columns = ['subroi', 'numvert', 'surfarea', 'grayvol',
                  'thickavg', 'thickstd',
                  'meancurv', 'gauscurv', 'foldind', 'curvind']
    df['side'] = df['subroi'].str[:2]
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
            print command
            output=os.popen(re.sub('\s+',' ',command)).read()
            pbar.update((num/totalNum) * 100)
            num+=1

            print 'output', output
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
    '''

    for side in ['lh','rh']:
        for cortex, rois in roiDict.iteritems():
            inLabelLocs = [join(fsDir,'tmp',side+'.'+x+'.label') for x in rois]
            inLabelForms = ' '.join(['-i '+x for x in inLabel])
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

    for side in ['lh','rh']:
        command = 'mri_annotation2label \
            --subject {basename} \
            --hemi {side} \
            --outdir {outDir} \
                2>/dev/null'.format(basename=fsDirName, 
                                    side=side, 
                                    outDir=labelOutDir)

        os.popen(re.sub('\s+',' ',command)).read()

        command = 'mri_annotation2label \
            --subject {basename} \
            --hemi {side} \
            --outdir {outDir} \
            --ctab {outDir}/{side}_ctab.txt \
                2>/dev/null'.format(basename=fsDirName, 
                                    side=side, 
                                    outDir=labelOutDir)

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

    ## FREESURFER settings
    #os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
    #os.environ["SUBJECTS_DIR"] = dirname(fsDir)

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
        '-f', '--fsDir',
        help='Freesurfer directory location',
        default=os.getcwd())


    args = parser.parse_args()

    os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
    os.environ["SUBJECTS_DIR"] = dirname(dirname(args.fsDir))

    freesurferSummary(args)

