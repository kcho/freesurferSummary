#!/ccnc_bin/mini_env/bin/python
from __future__ import division
__author__ = 'kcho'
import re
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import valueSwap
import textwrap
sys.path.append('/ccnc_bin')
import ccncpy.ccncpy as ccncpy
from progressbar import ProgressBar
import time


plt.style.use('ggplot')

def get_freesurferDir(dir):
    FS_description= ['bem','mri','scripts',
                     'src','stats','surf','tmp']

    freesurfer_dir = ccncpy.subDirSearch(FS_description,
                                         dir)

    if len(freesurfer_dir) > 1:
        print 'Please choose one data'
        print '======================'

        for num, i in enumerate(freesurfer_dir):
            print '{number} : {location}'.format(number=num,
                                                 location = i)
            choice = raw_input('\t: ')

        return ''.join(freesurfer_dir[choice])
    else:
        return ''.join(freesurfer_dir)


def main(subject_loc, background_subject_locs, graph, meanDfLoc,verbose, brain):


    ##########################################################
    # Freesurfer setting
    ##########################################################

    ##########################################################
    # Summarize cortical thickness
    ##########################################################
    # if no mean table is given but background list is given
    if meanDfLoc == False and background_subject_locs != None:
        print 'Calculating cortical thickness in {0}'.format(background_subject_locs)
        # make mean table
        meanDfName = raw_input('Name of the background subject : ')

        if verbose:
            meanDfList = []
            for subjectDir in background_subject_locs:
                meanDfList.append(collectStats_v2(subjectDir))
            meanDf = pd.concat(meanDfList)
        else:
            meanDf = collectStats(background_subject_locs)

    # if no mean table is given, and backround list is empty
    # use subject_loc as the background
    elif meanDfLoc == False and background_subject_locs == None:
        print 'No background subjects are given, now running only {0}'.format(subject_loc)
        meanDfName = ''

        if verbose:
            meanDf = collectStats_v2([subject_loc])
        else:
            meanDf = collectStats([subject_loc])

    # if meanDfLoc is given
    else:
        print 'Now comparing with mean_thickness.csv in /ccnc_bin/meanThickness/'
        meanDfName = 'NOR'
        #meanDf = pd.read_csv('/ccnc_bin/meanThickness/new_mean_OCT15.csv')
        if verbose:
            meanDf = pd.read_csv('/ccnc_bin/meanThickness/detailed_mean_2015_12_28.csv', index_col=0)
        else:
            meanDf = pd.read_csv('/ccnc_bin/meanThickness/new_mean_OCT15.csv', index_col=0)

    ##########################################################
    # Get roi dict : 8 cortex each hemisphere
    ##########################################################
    roiDict = get_cortical_rois()

    ##########################################################
    # annotation2label --> merge labels --> freesurfer/tmp
    ##########################################################
    
    if graph:
        if verbose:
            thicknessDf = collectStats_v2([os.path.dirname(freesurfer_dir)])#background_subject_locs)
            draw_thickness_detailed(thicknessDf,meanDf,os.path.basename(subject_loc), meanDfName, subject_loc)
        else:
            thicknessDf = collectStats([os.path.dirname(freesurfer_dir)])#background_subject_locs)
            draw_thickness(thicknessDf,meanDf,os.path.basename(subject_loc), meanDfName, subject_loc)


    if brain:
        if verbose:
            valueSwap.main(freesurfer_dir, 
                    os.path.join(freesurfer_dir,
                        'tmp/thick_kev_detailed.csv'))
        else:
            print 'Please select verbose mode'
            pass

    print "*"*80
    print 'Outputs are saved in /ccnc/mri_team'
    print "*"*80


def draw_thickness_detailed(infoDf, meanDf, subjName, meanDfName):
    roiOrder = ['LPFC', 'OFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC', 'OCC']

    def getRegion(roi):
        roiDict = get_cortical_rois()
        for region, roiList in roiDict.iteritems():
            if roi in roiList:
                return region

    infoDf['roi'] = infoDf.subroi.str[3:]
    infoDf['region'] = infoDf.roi.apply(getRegion)
    meanDf['roi'] = meanDf.subroi.str[3:]
    meanDf = meanDf.groupby(['roi','side']).mean().reset_index()
    meanDf['region'] = meanDf.roi.apply(getRegion)

    gb = infoDf.groupby('region')
    infoDf = pd.concat([gb.get_group('LPFC'),
                        gb.get_group('OFC'),
                        gb.get_group('MPFC'),
                        gb.get_group('LTC'),
                        gb.get_group('MTC'),
                        gb.get_group('SMC'),
                        gb.get_group('PC'),
                        gb.get_group('OCC')])

    gbmean = meanDf.groupby('region')
    meanDf = pd.concat([gbmean.get_group('LPFC'),
                        gbmean.get_group('OFC'),
                        gbmean.get_group('MPFC'),
                        gbmean.get_group('LTC'),
                        gbmean.get_group('MTC'),
                        gbmean.get_group('SMC'),
                        gbmean.get_group('PC'),
                        gbmean.get_group('OCC')])

    gb = infoDf.groupby('side')
    label = infoDf.subroi.str[3:].unique()

    fig = plt.figure(figsize=(22,12))
    fig.suptitle("Cortical thickness in all regions", fontsize=20)

    lh_g = plt.subplot2grid((2,2),(0, 0), colspan=2)
    rh_g = plt.subplot2grid((2,2),(1, 0), colspan=2)
    lh_g.plot(gb.get_group('lh')['thickavg'],'r',label=subjName)

    lh_g.plot(meanDf.groupby('side').get_group('lh')['thickavg'],'r--',label=meanDfName)

    # error bar

    eb1 = lh_g.errorbar(range(len(meanDf.roi.unique())),
                        meanDf.groupby('side').get_group('lh')['thickavg'],
                        meanDf.groupby('side').get_group('lh')['thickstd'],
                        linestyle='None',
                        marker='^')
    eb1[-1][0].set_linestyle('--')

    lh_g.set_xlabel('Left', fontsize=16)

    lh_g.set_ylim(1.0, 4.7)

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
                        gb.get_group('lh'),
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




    rh_g.plot(gb.get_group('rh')['thickavg'],'b',label=subjName)
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
    rh_g.set_ylim(1, 4.7)
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
                        gb.get_group('rh'),
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
    fig.savefig('/ccnc/mri_team/'+os.path.basename(subjName))

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



def getInfoFromLabel(freesurfer_dir,roiDict):
    infoDict={}

    pbar = ProgressBar().start()
    totalNum = 2 * len(roiDict.keys())
    num = 1
    for side in ['lh','rh']:
        for cortex, rois in roiDict.iteritems():
            if len(rois) > 1:
                command = 'mris_anatomical_stats \
                -l {loc}/{side}_{cortex} {name} {side} 2>/dev/null'.format(
                    loc=os.path.join(freesurfer_dir,'tmp'),
                    side=side,
                    cortex=cortex,
                    name=os.path.basename(freesurfer_dir)
                )
            else:
                command = 'mris_anatomical_stats \
                -l {loc}/{side}.{cortex}.label {name} {side} 2>/dev/null'.format(
                    loc=os.path.join(freesurfer_dir,'tmp'),
                    side=side,
                    cortex=cortex,
                    name=os.path.basename(freesurfer_dir)
                )

            output=os.popen(re.sub('\s+',' ',command)).read()
            pbar.update((num/totalNum) * 100)
            num+=1

            thickness = re.search('thickness\s+=\s+(\S+)\s+mm\s+\S+\s+(\S+)', output).group(1,2)
            numvert = re.search('number of vertices\s+=\s+(\S+)', output).group(1)
            surfarea = re.search('total surface area\s+=\s+(\S+)', output).group(1)
            grayvol = re.search('total gray matter volume\s+=\s+(\S+)', output).group(1)
            meancurv = re.search('average integrated rectified mean curvature\s+=\s+(\S+)', output).group(1)
            gauscurv = re.search('average integrated rectified Gaussian curvature\s+=\s+(\S+)', output).group(1)
            foldind = re.search('folding index\s+=\s+(\S+)', output).group(1)
            curvind = re.search('intrinsic curvature index\s+=\s+(\S+)', output).group(1)

            thickness = tuple([float(x) for x in thickness])
            for i in numvert, surfarea, grayvol, meancurv, gauscurv, foldind, curvind:
                i = float(i)

            infoDict[side+'_'+cortex] = [numvert, surfarea, grayvol,
                                              thickness[0], thickness[1],
                                              meancurv, gauscurv, foldind, curvind]
            #print cortex, rois, thicknessDict[side+'_'+cortex]
    pbar.finish()
    return infoDict


def mergeLabel(freesurfer_dir, roiDict):
    for side in ['lh','rh']:
        for cortex, rois in roiDict.iteritems():
            command = 'mri_mergelabels {inLabel} -o {outLabel} 2>/dev/null'.format(
                inLabel = ' '.join(['-i '+os.path.join(freesurfer_dir,'tmp',side+'.'+x+'.label') for x in rois]),
                outLabel = os.path.join(freesurfer_dir,'tmp',side+'_'+cortex))
            os.popen(command).read()


def makeLabel(freesurfer_dir):
    '''
    Run mri_annotation2label for lh and rh hemisphere.
    Creates labels in $freesurfer_dir/tmp
    '''
    for side in ['lh','rh']:
        command = 'mri_annotation2label \
            --subject {basename} \
            --hemi {side} --outdir {outDir} 2>/dev/null'.format(basename=os.path.basename(freesurfer_dir),
                                                    side=side,
                                                    outDir=os.path.join(freesurfer_dir,'tmp'))
        os.popen(re.sub('\s+',' ',command)).read()

        command = 'mri_annotation2label \
            --subject {basename} \
            --hemi {side} --outdir {outDir} --ctab {outDir}/{side}_ctab.txt 2>/dev/null'.format(basename=os.path.basename(freesurfer_dir),
                                                    side=side,
                                                    outDir=os.path.join(freesurfer_dir,'tmp'))
        os.popen(re.sub('\s+',' ',command)).read()

def collectStats_v2(freesurfer_dir):
    '''
    CollectStats version 2
    Summarise cortical thickness in more than one subjects.
    'background_subject_locs' should given in python list format.
    For each background,
    1. Creates labels using makeLabel
    2. Merges labels according to the roiDict using mergeLabel
    3. Estimates cortical thickness in each merged labels (dict)
    4. Converts dict to pandas Dataframe using dictWithTuple2df
    5. save df to freesurfer_dir/tmp/thick_kev_detailed.csv
    Than mean df is returned.
    '''
    #collectStats('ha','ho',ha')

    # cortical regions as a dictionary
    roiDict = get_cortical_rois_detailed()

    # freesurfer sub-directory description

    if not os.path.isfile(os.path.join(freesurfer_dir,'tmp','thick_kev_detailed_new.csv')):
        makeLabel(freesurfer_dir)
        mergeLabel(freesurfer_dir, roiDict)
        infoDict = getInfoFromLabel(freesurfer_dir, roiDict)
        infoDf = dictWithTuple2df(infoDict)
        infoDf.to_csv(os.path.join(freesurfer_dir,'tmp','thick_kev_detailed_new.csv'))

    else:
        infoDf = pd.read_csv(os.path.join(freesurfer_dir,'tmp','thick_kev_detailed_new.csv'),
                             index_col=0)

    # return mean
    return infoDf

def collectStats(background_subject_locs):
    '''
    Summarise cortical thickness in more than one subjects.
    'background_subject_locs' should given in python list format.
    For each background,
    1. Creates labels using makeLabel
    2. Merges labels according to the roiDcit using mergeLabel
    3. Estimates cortical thickness in each merged labels (dict)
    4. Converts dict to pandas Dataframe using dictWithTuple2df
    5. save df to freesurfer_dir/tmp/thick_kev.csv
    Than mean df is returned.
    '''
    #collectStats('ha','ho',ha')

    # 8 cortex dictionary
    roiDict = get_cortical_rois()

    subjectDict = {}
    for background in background_subject_locs:

        # freesurfer sub-directory description
        FS_description= ['bem','mri','scripts','src','stats','surf','tmp']
        freesurfer_dir = ccncpy.subDirSearch(FS_description, background)

        if len(freesurfer_dir) > 1:
            sys.exit(re.sub('\s+',' ',
            'There are more than 1 freesurfer directory \
                    under {0}'.format(subject_loc)))
        else:
            freesurfer_dir = ''.join(freesurfer_dir)

        if not os.path.isfile(os.path.join(freesurfer_dir,'tmp','thick_kev.csv')):
            makeLabel(freesurfer_dir)
            mergeLabel(freesurfer_dir, roiDict)
            infoDict = getInfoFromLabel(freesurfer_dir, roiDict)
            infoDf = dictWithTuple2df(infoDict)
            infoDf.to_csv(os.path.join(freesurfer_dir,'tmp','thick_kev.csv'))
            
        else:
            infoDf = pd.read_csv(os.path.join(freesurfer_dir,'tmp','thick_kev.csv'))

        subjectDict[background] = infoDf


    # sum up dataframes in a subjectDict dictionary
    finalDf = pd.concat([x for x in subjectDict.values()])

    # return mean
    return finalDf.groupby('subroi').mean().reset_index()



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
            '''.format(codeName=os.path.basename(__file__),
                       in_put = 'subjectLoc',
                       output = 'outLoc')))

    parser.add_argument(
        '-i', '--inputDir',
        help='Subject location',
        default=os.getcwd())

    parser.add_argument(
        '-c', '--createMeanFrom',
        help='Subject locations to create mean',
        nargs='+',
        default=False)

    parser.add_argument(
        '-x', '--gender',
        help='M or F',
        )

    parser.add_argument(
        '-s', '--saveMeanDf',
        help='output location of the meanDf created from -c',
        default=False)

    parser.add_argument(
        '-g', '--graph',
        help='Draw graph',
        default=True,
        action='store_true')

    parser.add_argument(
        '-m', '--meanDf',
        help='meanDf',
        default=True,
        action='store_true')
        #default = '/ccnc_bin/meanThickness/mean_thickness.csv')

    parser.add_argument(
        '-v', '--verbose',
        help='Use detailed ROIs',
        default=True,
        action='store_true')
        #default = '/ccnc_bin/meanThickness/mean_thickness.csv')

    parser.add_argument(
        '-p', '--brain',
        help='make brain picture',
        default=True,
        action='store_true')

    args = parser.parse_args()

    # Main subject
    main_freesurferDir = get_freesurferDir(os.path.abspath(args.inputDir))
    print main_freesurferDir
    infoDf = collectStats_v2(main_freesurferDir)
    print infoDf

    # Freesurfer environment Settings
    # os.environ["FREESURFER_HOME"] = '/Applications/freesurfer'
    # os.environ["SUBJECTS_DIR"] = '{0}'.format(os.path.dirname(main_freesurferDir))

    # Mean Df
    if args.createMeanFrom:
        freesurferList = subjDirs_to_fsDirs(args.createMeanFrom)
        concatDf = concatFsDf(freesurferList)
        meanDf = concatDf_to_meanDf(concatDf)
        if args.saveMeanDf:
            meanDf.to_csv(args.saveMeanDf)
    elif args.gender.upper() == 'M':
        meanDf = pd.read_csv('/ccnc_bin/meanThickness/male_df.csv', index_col=0)
    elif args.gender.upper() == 'F':
        meanDf = pd.read_csv('/ccnc_bin/meanThickness/female_df.csv', index_col=0)
    else:
        meanDf = pd.read_csv('/ccnc_bin/meanThickness/detailed_mean_2015_12_28.csv', index_col=0)

    print meanDf

    # Graph
    draw_thickness_detailed(infoDf,
                            meanDf,
                            os.path.basename(os.path.abspath(args.inputDir)),
                            'CCNC_mean')





