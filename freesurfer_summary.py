#!/ccnc/anaconda2/bin/python
from __future__ import division
import os
import re
import seaborn as sns
from os.path import join, basename, dirname, isfile
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import cm 
import argparse
import itertools
import textwrap
from progressbar import ProgressBar
import time
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')
#import tksurferCapture
pd.options.mode.chained_assignment = None  # default='warn'


__author__ = 'kcho'
plt.style.use('ggplot')

def demo_match(age, age_range, sex, all_data_Loc):
    '''
    by yb
    '''
    matching = pd.read_csv(all_data_Loc)

    upper = int(age) + age_range
    lower = int(age) - age_range

    matching_age = matching[(matching['age'] >= lower) & (matching['age'] <= upper)]
    matching_sex = matching_age[matching_age['sex'] == sex]

    # see column names
    if 'thickness' in matching.columns:
        matched = matching_sex[['side', 'roi', 'region', 'thickness','thicknessstd', 'volume', 'subject']]
        matched_mean = matched.groupby(['roi','side','region']).mean().reset_index()
        matched_std = matched.groupby(['roi','side','region']).std().reset_index()
        matched_std = matched_std.rename(columns={'thickness': 'thickstd', 
                                                  'thicknessstd': 'per_region_std',
                                                  'volume': 'volumestd'})
        matched_mean_std = pd.merge(matched_mean, matched_std, how='inner')

    elif 'thickness' not in matching.columns:
        matched = matching_sex[['roi', 'volume', 'region', 'subject']]
        matched_mean = matched.groupby(['roi', 'region']).mean().reset_index()
        matched_std = matched.groupby(['roi','region']).std().reset_index()
        matched_std = matched_std.rename(columns={'volume': 'volumestd'})
        matched_mean_std = pd.merge(matched_mean, matched_std, how='inner')
    
    return matched_mean_std

def getGroupMeanInfo(meanDfLoc):
    '''
    Read cortical thickness of all control subjects
    '''

    df = pd.read_csv(meanDfLoc, index_col=0)
    if 'roi' not in df.columns or 'region' not in df.columns:
        df['roi'] = df.subroi.str[3:]
        df['region'] = df.roi.apply(getRegion)
        df.to_csv(meanDfLoc)

    # meanDf to have averaged values for all groups
    meanDf = df.groupby(['roi','side','region']).mean().reset_index()
    meanDf.columns = ['roi','side','region', 'thickavg','thickstd']
    df_name = 'CCNC_mean'
    return meanDf, df_name

def makeMean(inputDirs):
    '''
    inputDirs is the freesurfer directories to be merged into 
    background information
    '''
    cortical_dfs = pd.DataFrame()
    subcortical_dfs = pd.DataFrame()
    pbar = ProgressBar().start()

    # mri spreadsheet
    mri_excel_loc = '/home/kangik/Dropbox/MRI/MRI.xls'
    mri_excel = pd.ExcelFile(mri_excel_loc)
    nor_df = mri_excel.parse('NOR')
    nor_df = nor_df[(nor_df.timeline=='baseline')][['folderName','age','sex']]

    for fsDirNum, fsDir in enumerate(inputDirs):
        pbar.update((fsDirNum/len(inputDirs)) * 100)
        # aparcstats2table
        cortical_df = aparcstats2table(os.path.abspath(fsDir), 'aparc')
        cortical_df['subject'] = basename(fsDir)
        cortical_dfs = pd.concat([cortical_dfs, cortical_df])

        # asegstats2table
        # df.columns = ['roi', 'volume', 'region']
        # regions are subocortex
        subcortical_df =  asegstats2table(os.path.abspath(fsDir))
        subcortical_df['subject'] = basename(fsDir)
        subcortical_dfs = pd.concat([subcortical_dfs, subcortical_df])
    pbar.finish()

    cortical_dfs = pd.merge(cortical_dfs, nor_df,
                            left_on='subject',
                            right_on='folderName',
                            how='left')
    subcortical_dfs = pd.merge(subcortical_dfs, nor_df,
                            left_on='subject',
                            right_on='folderName',
                            how='left')

    cortical_dfs.to_csv('all_cortical_dfs_{date}.csv'.format(
        date = time.strftime("%Y_%m_%d")))
    subcortical_dfs.to_csv('all_subcortical_dfs_{date}.csv'.format(
        date = time.strftime("%Y_%m_%d")))
    
    mean_cortical_dfs = cortical_dfs.groupby(['side', 'roi', 'region']).mean()
    mean_subcortical_dfs = subcortical_dfs.groupby(['roi', 'region']).mean()

    # thickstd : standard deviation of thickness across subjects
    # thicknessstd : mean of standard deviations of thicknness in each regions
    mean_cortical_dfs['volumestd'] = cortical_dfs.groupby(['side', 'roi', 'region']).std()['volume']
    mean_cortical_dfs['thickstd'] = cortical_dfs.groupby(['side', 'roi', 'region']).std()['thickness']
    mean_subcortical_dfs['volumestd'] = subcortical_dfs.groupby(['roi', 'region']).std()['volume']

    mean_cortical_dfs.to_csv('mean_cortical_dfs_{date}.csv'.format(
        date = time.strftime("%Y_%m_%d")))
    mean_subcortical_dfs.to_csv('mean_subcortical_dfs_{date}.csv'.format(
        date = time.strftime("%Y_%m_%d")))
    #print mean_subcortical_dfs

def freesurferSummary(inputDirs, nameList=False, ageList=False, genderList=False, ageRange=3, colorList=False, nobackground=False):
    '''
    Summarizes freesurfer outputs using matplotlib
    - Cortical thickness and volume
    - Subcortical volumes & ICV
    '''

    # Collect infoDfs and names
    cortical_dfs = []
    subcortical_dfs = []
    fsNames = []

    for fsDirNum, (fsDir, age, gender) in enumerate(zip(inputDirs,
                                                            ageList,
                                                            genderList)):
        print(fsDir)
        if nameList:
            subjNames = '{name} {gender} {age}'.format(
                name = nameList[fsDirNum],
                gender = gender,
                age = age)
        else:
            subjNames = raw_input('{0} Subject initial :'.format(fsDir))
        fsNames.append(subjNames)

        # Label approach
        #fsInfoDf.append(collectStats(os.path.abspath(fsDir)))

        # aparcstats2table
        cortical_df = aparcstats2table(os.path.abspath(fsDir), 'aparc')
        cortical_df['subject'] = subjNames
        cortical_dfs.append(cortical_df)

        # asegstats2table
        # df.columns = ['roi', 'volume', 'region']
        # regions are subocortex
        subcortical_df =  asegstats2table(os.path.abspath(fsDir))
        subcortical_df['subject'] =  subjNames
        subcortical_dfs.append(subcortical_df)

    if nobackground==False: # Mean graph option turned on
        # Read CCNC healthy control information
        for age, gender in zip(ageList, genderList):
            # csv with all subjects' data
            mean_cortical_df_loc = '/Volume/CCNC_BI_3T/freesurfer/NOR/all_cortical_dfs_2017_07_04.csv'
            mean_subcortical_df_loc = '/Volume/CCNC_BI_3T/freesurfer/NOR/all_subcortical_dfs_2017_07_04.csv'

            # Yoobin function added here
            mean_cortical_df = demo_match(age, ageRange, gender, 
                                          mean_cortical_df_loc)
            mean_subcortical_df = demo_match(age, ageRange, gender, 
                                          mean_subcortical_df_loc)

            # Add CCNC HCs informat
            cortical_dfs.append(mean_cortical_df)
            subcortical_dfs.append(mean_subcortical_df)

            # Mean graph name with age and gender information
            mgName = 'CCNC_mean {age} age range: {ageRange} {gender}'.format(age=age, 
                                                                    ageRange=ageRange, 
                                                                    gender=gender)
            fsNames.append(mgName)

    # Make line plots of cortical thickness for each hemisphere
    draw_thickness_list(cortical_dfs, fsNames, colorList)
    draw_volume_list(cortical_dfs, fsNames, colorList)
    draw_subcortical_volume_list(subcortical_dfs, fsNames, colorList)

    # tksurferCapture.main(fsDir, join(fsDir,
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

def draw_subcortical_volume_list(infoDfList, nameList, colorList):
    '''
    Draw graph that summarizes subcortical volumes from 
    freesurfer outputs.

    The last items in the given input lists are that of meandf.

    1. Separate out data that has side information.
    2. Divide the rest of volumes in
       to four graphs according to the volumes
    '''
    number_of_graphs = 4

    # mean infoDf
    meanDf = infoDfList[-1]
    # remove regions
    regions_to_remove = ['BrainSegVol', 'BrainSegVolNotVent', 
                         'BrainSegVolNotVentSurf', 'MaskVol',
                         'SupraTentorialVol', 'SupraTentorialVolNotVent',
                         'SupraTentorialVolNotVentVox',
                         'SubCortGrayVol', 'Brain-Stem', 
                         'Left-Cerebellum-Cortex',
                         'Right-Cerebellum-Cortex',
                         'Left-Cerebellum-White-Matter',
                         'Right-Cerebellum-White-Matter',
                         'BrainSegVol-to-eTIV',
                         'MaskVol-to-eTIV',
                         '5th-Ventricle',
                         'Optic-Chiasm',
                         'SurfaceHoles', 
                         'lhSurfaceHoles', 'rhSurfaceHoles',
                         'WM-hypointensities',
                         'non-WM-hypointensities',
                         'Right-non-WM-hypointensities',
                         'Left-non-WM-hypointensities',
                         'Right-WM-hypointensities',
                         'Left-WM-hypointensities',
                         'Right-vessel', 'Left-vessel']

    meanDf = meanDf[(~meanDf.roi.isin(regions_to_remove))]
    rois_with_side = meanDf.roi[meanDf.roi.str.contains('(Left|Right)')]
    meanDf_withside = meanDf[meanDf.roi.isin(rois_with_side)]
    meanDf_withoutside = meanDf[~meanDf.roi.isin(rois_with_side)]

    # kmeans
    # volume range divisions
    graphNum = 4
    # data to divide
    subcortical_volumes_mean = np.array(meanDf_withoutside.volume).reshape(-1,1)
    # settin up and training the kmeans model
    kmeans = KMeans(n_clusters=graphNum, 
                    random_state=0).fit(subcortical_volumes_mean)
    # chage here later
    # gives error without
    # pd.options.mode.chained_assignment = None  # default='warn'
    meanDf_withoutside['gtype'] = kmeans.labels_

    # Graphs
    fig = plt.figure(figsize=(22,12), facecolor='white')
    fig.suptitle("Sub-cortical Volume Summary", fontsize=20)

    # axes
    # Upper axes : Data with side information
    gs = gridspec.GridSpec(3,4)
    left_ax = plt.subplot(gs[0,:])
    right_ax = plt.subplot(gs[1,:])
    axes_side = [left_ax, right_ax]
    # Lower axes : Data without side information
    # K-means clustered
    ax1 = plt.subplot(gs[2,0])
    ax2 = plt.subplot(gs[2,1])
    ax3 = plt.subplot(gs[2,2])
    ax4 = plt.subplot(gs[2,3])
    axes = [ax1, ax2, ax3, ax4]

    # color definition
    color=iter(cm.rainbow(np.linspace(0,1,len(infoDfList))))

    sig_diff_region = []
    # For every subject information dataframes
    for infoDfNum, (infoDf, subjName) in enumerate(zip(infoDfList, nameList)):
        # Setting graph colours
        if colorList:
            # if meanDf is on, colorList is one shorter.
            try:
                c = colorList[infoDfNum]
            except:
                c = 'b'
        else:
            c = next(color)

        # For subcortical main structures, 
        # divided left and right
        df = infoDf[(~infoDf.roi.isin(regions_to_remove))]
        df_side = df[(df.roi.isin(rois_with_side))]
        df_side.loc[:, 'side'] = df_side.roi.str.extract('(Left|Right)', expand=False)
        df_side.loc[:, 'merged_roi'] = df.roi.str.extract('\w{4,5}-(\S*)', expand=False)
        df_without_side = df[~(df.roi.isin(rois_with_side))]

        for sideNum, side in enumerate(['Left', 'Right']):
            side_df = df_side.groupby('side').get_group(side)
            side_df = side_df.sort_values('merged_roi')
            ax = axes_side[sideNum]
            if 'CCNC_mean' in subjName:
                #ax.plot(infodf.thickavg, '--', c=c, label=subjName)
                eb = ax.errorbar(range(len(side_df.merged_roi.unique())),
                                  side_df.volume,
                                  side_df.volumestd*2,
                                  marker='^',
                                  label=subjName, 
                                  color='b',
                                  ecolor='b',
                                  alpha=0.7)
                plotline, caplines, barlinecols = eb
                barlinecols[0].set_linestyle('--')
                plotline.set_linestyle('--')
            else:
                ax.plot(range(len(side_df.merged_roi)), 
                        side_df.volume, 
                        c=c, 
                        marker='o',
                        label=subjName)

                ### annotation
                mergedDf = pd.merge(meanDf_withside,
                                    side_df,
                                    on='roi',
                                    how='right')
                mergedDf['mean_sub_indv'] = mergedDf.volume_x - mergedDf.volume_y

                # Reorder Dfs
                mergedDf = mergedDf.sort_values('roi')
                #mergedDf = reorder_df(mergedDf, 'region', roiOrder)

                #for sigNum, row in enumerate(mergedDf[abs(mergedDf['mean_sub_indv']) > 0.5].iterrows()):
                # greater than two standard deviation
                arrowSize = np.mean(side_df.volume.std())
                for sigNum, row in enumerate(mergedDf[abs(mergedDf['mean_sub_indv']) > mergedDf['volumestd']*2].iterrows()):
                    if (sigNum+1) % 2 == 0:
                        sign = 1
                        diff = arrowSize/2
                    else:
                        sign = -1
                        diff = 0

                    if row[1].mean_sub_indv < 0:
                        col = 'green'
                        sign = 1
                    else:
                        col = 'red'

                    textLoc_y = row[1].volume_y + (sign * arrowSize) + diff

                    if row[1].roi not in sig_diff_region:
                        ax.annotate(row[1].roi,
                                    xy = (row[0], row[1].volume_y), 
                                    xycoords='data',
                                    xytext = (row[0], textLoc_y),
                                    textcoords='data',
                                    arrowprops = dict(facecolor=col, 
                                                      shrinkB=5,
                                                      alpha=0.5), 
                                    horizontalalignment='center', 
                                    fontsize=15)
                        sig_diff_region.append(row[1].roi)
            ax.grid(False)
            ax.set_xlim(-.5, len(side_df.merged_roi)+0.5)
            ax.set_xticks(range(len(side_df.merged_roi)))
            ax.set_xticklabels(side_df.merged_roi)
            ax.patch.set_facecolor('white')
            ax.autoscale(enable=True, axis='x')

        # Bottom three graphs
        for axNum, ax in enumerate(axes):
            # List of rois in subgroups
            roi_order = meanDf_withoutside[(meanDf_withoutside.gtype==axNum)].roi

            # Subgroup roi df
            df = infoDf[(infoDf.roi.isin(roi_order))]
            df = df.sort_values('roi')
            #df = df.set_index('roi').reindex(roi_order).reset_index() # sort
            #df['gtype'] = axNum

            if 'CCNC_mean' in subjName:
                eb = ax.errorbar(range(len(df.roi.unique())),
                                  df.volume,
                                  df.volumestd*2,
                                  marker='^',
                                  label=subjName, 
                                  color='b',
                                  ecolor='b',
                                  alpha=0.7)
                plotline, caplines, barlinecols = eb
                barlinecols[0].set_linestyle('--')
                plotline.set_linestyle('--')
            else:
                ax.plot(range(len(roi_order)), 
                        df.volume, 
                        c=c, 
                        marker='o',
                        label=subjName)

                ### annotation
                meanDf_set = meanDf.loc[(meanDf.roi.isin(roi_order))]
                mergedDf = pd.merge(meanDf_set,
                                    df,
                                    on=['roi','region'],
                                    how='inner')
                mergedDf['mean_sub_indv'] = mergedDf.volume_x - mergedDf.volume_y

                # Reorder Dfs
                mergedDf = mergedDf.sort_values('roi')
                #mergedDf = reorder_df(mergedDf, 'region', roiOrder)

                #for sigNum, row in enumerate(mergedDf[abs(mergedDf['mean_sub_indv']) > 0.5].iterrows()):
                # greater than two standard deviation
                arrowSize = np.mean(df.volume.std())
                for sigNum, row in enumerate(mergedDf[abs(mergedDf['mean_sub_indv']) > mergedDf['volumestd']*2].iterrows()):
                    if (sigNum+1) % 2 == 0:
                        sign = 1
                        diff = arrowSize/2
                    else:
                        sign = -1
                        diff = 0

                    if row[1].mean_sub_indv < 0:
                        col = 'green'
                        sign = 1
                    else:
                        col = 'red'

                    textLoc_y = row[1].volume_y + (sign * arrowSize) + diff

                    if row[1].roi not in sig_diff_region:
                        ax.annotate(row[1].roi,
                                    xy = (row[0], row[1].volume_y), 
                                    xycoords='data',
                                    xytext = (row[0], textLoc_y),
                                    textcoords='data',
                                    arrowprops = dict(facecolor=col, 
                                                      shrinkB=5,
                                                      alpha=0.5), 
                                    horizontalalignment='center', 
                                    fontsize=15)
                        sig_diff_region.append(row[1].roi)

                # axis settings
                    #ax = axes[snum]
                ax.patch.set_facecolor('white')
                # Graph settings
                #ax.set_ylim(-5000, 35000)
                #ax.set_xlabel('Subcortical ROI', fontsize=16)
                ax.set_xticks(range(len(roi_order)))
                #ax.set_xticklabels(['' for x in roiList])
                ax.set_xlim(-.5, len(df.roi)+0.5) 
                ax.grid(False)
                ax.autoscale(enable=True, axis='x')

                ax.set_xticklabels(df.roi)
                labels = ax.get_xticklabels()
                if len(roi_order) > 3:
                    plt.setp(labels, rotation=30)
                plt.tight_layout(pad=7, w_pad=3, h_pad=0.2)


    ax2.set_xticklabels(['Total Intra-cranial Volume'], rotation=0)
    ax2.set_ylim(1000000, 2200000)
    left_ax.legend()
    legend = left_ax.legend(frameon = 1)
    frame = legend.get_frame()
    frame.set_facecolor('white')

    #fig.show()
    fname = 'subcortical_volume_summary.png'
    fig.savefig(fname)
    print('feh %s' %join(os.getcwd(), fname))

def draw_volume_list(infoDfList, nameList, colorList):
    # graph order from left
    roiOrder = ['LPFC', 'OFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC', 'OCC']

    fig, axes = plt.subplots(nrows=2, figsize=(22,12), facecolor='white')
    fig.suptitle("Cortical Volume Summary", fontsize=20)

    # color definition
    color=iter(cm.rainbow(np.linspace(0,1,len(infoDfList))))

    # mean infoDf
    meanDf = infoDfList[-1]
    roiList = meanDf.roi.unique()

    sig_diff_region = []
    for infoDfNum, (infoDf, subjName) in enumerate(zip(infoDfList, nameList)):
        infoDf_gb = infoDf.groupby('side')
        
        if colorList:
            # if meanDf is on, colorList is one shorter.
            try:
                c = colorList[infoDfNum]
            except:
                c = 'b'
        else:
            c = next(color)

        for snum, side in enumerate(['lh', 'rh']):
            infodf = infoDf_gb.get_group(side)

            # Reorder Dfs
            infodf = infodf.sort_values(['roi','side'])
            infodf = reorder_df(infodf, 'region', roiOrder)

            ax = axes[snum]
            if 'CCNC_mean' in subjName:
                #ax.plot(infodf.thickavg, '--', c=c, label=subjName)
                eb = ax.errorbar(range(len(roiList)),
                                  infodf.volume,
                                  infodf.volumestd*2,
                                  marker='^',
                                  label=subjName, 
                                  color='b',
                                  ecolor='b',
                                  alpha=0.7)
                plotline, caplines, barlinecols = eb
                barlinecols[0].set_linestyle('--')
                plotline.set_linestyle('--')
            else:
                ax.plot(infodf.volume, c=c, marker='o', label=subjName)

                ### annotation
                mergedDf = pd.merge(meanDf,
                                    infodf,
                                    on=['roi','side','region'],
                                    how='inner')
                mergedDf['mean_sub_indv'] = mergedDf.volume_x - mergedDf.volume_y

                # Reorder Dfs
                mergedDf = mergedDf.sort_values(['roi','side'])
                mergedDf = reorder_df(mergedDf, 'region', roiOrder)

                #for sigNum, row in enumerate(mergedDf[abs(mergedDf['mean_sub_indv']) > 0.5].iterrows()):
                # greater than two standard deviation
                for sigNum, row in enumerate(mergedDf[abs(mergedDf['mean_sub_indv']) > mergedDf['volumestd']*2].iterrows()):
                    if (sigNum+1) % 2 == 0:
                        sign = 1
                        diff = 2500
                    else:
                        sign = -1
                        diff = 0

                    if row[1].mean_sub_indv < 0:
                        col = 'green'
                        sign = 1
                    else:
                        col = 'red'

                    textLoc_y = row[1].volume_y + (sign * 5000) + diff

                    if row[1].roi not in sig_diff_region:
                        ax.annotate(row[1].roi,
                                    xy = (row[0], row[1].volume_y), 
                                    xycoords='data',
                                    xytext = (row[0], textLoc_y),
                                    textcoords='data',
                                    arrowprops = dict(facecolor=col, 
                                                      shrinkB=5,
                                                      alpha=0.5), 
                                    horizontalalignment='center', 
                                    fontsize=15)
                        sig_diff_region.append(row[1].roi)

    # axis settings
    for snum, side in enumerate(['lh', 'rh']):
        ax = axes[snum]
        ax.patch.set_facecolor('white')
        # Graph settings
        ax.set_ylim(-5000, 35000)
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

        regionToHighlight = roiOrder[1::2]
        xCoords = [np.where(roiOrder_array==x)[0] for x in regionToHighlight]

        for x in xCoords:
            ax.axvspan(x[0], x[-1], alpha=0.5)

        startNum = 0
        for region in roiOrder:
            x_starting_point = startNum
            startNum = startNum + len(roiDict[region])

            ax.text((x_starting_point-.5 + startNum-.5)/2, 1.2,
                    region,
                    horizontalalignment='center',
                    alpha=.4,
                    fontsize=15)

        ax.set_xticklabels(infodf.roi)
        labels = ax.get_xticklabels()
        plt.setp(labels, rotation=30)
        plt.tight_layout(pad=7, w_pad=3, h_pad=0.2)

    #fig.show()
    fname = 'volume_summary.png'
    fig.savefig(fname)
    print('feh %s' %join(os.getcwd(), fname))

def draw_thickness_list(infoDfList, nameList, colorList):
    # graph order from left
    roiOrder = ['LPFC', 'OFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC', 'OCC']

    fig, axes = plt.subplots(nrows=2, figsize=(22,12), facecolor='white')
    fig.suptitle("Cortical thickness in all regions", fontsize=20)

    # color definition
    color=iter(cm.rainbow(np.linspace(0,1,len(infoDfList))))

    # mean infoDf
    meanDf = infoDfList[-1]
    roiList = meanDf.roi.unique()

    sig_diff_region = []
    for infoDfNum, (infoDf, subjName) in enumerate(zip(infoDfList, nameList)):
        infoDf_gb = infoDf.groupby('side')
        
        if colorList:
            # if meanDf is on, colorList is one shorter.
            try:
                c = colorList[infoDfNum]
            except:
                c = 'b'
        else:
            c = next(color)

        for snum, side in enumerate(['lh', 'rh']):
            infodf = infoDf_gb.get_group(side)

            # Reorder Dfs
            infodf = infodf.sort_values(['roi','side'])
            infodf = reorder_df(infodf, 'region', roiOrder)

            ax = axes[snum]
            if 'CCNC_mean' in subjName:
                #ax.plot(infodf.thickavg, '--', c=c, label=subjName)
                eb = ax.errorbar(range(len(roiList)),
                                  infodf.thickness,
                                  infodf.thickstd*2,
                                  marker='^',
                                  label=subjName, 
                                  color='b',
                                  ecolor='b',
                                  alpha=0.7)
                plotline, caplines, barlinecols = eb
                barlinecols[0].set_linestyle('--')
                plotline.set_linestyle('--')
            else:
                ax.plot(infodf.thickness, c=c, marker='o', label=subjName)

                ### annotation
                mergedDf = pd.merge(meanDf,
                                    infodf,
                                    on=['roi','side','region'],
                                    how='inner')
                mergedDf['mean_sub_indv'] = mergedDf.thickness_x - mergedDf.thickness_y

                # Reorder Dfs
                mergedDf = mergedDf.sort_values(['roi','side'])
                mergedDf = reorder_df(mergedDf, 'region', roiOrder)

                # Greater than two standard deviation
                for sigNum, row in enumerate(mergedDf[abs(mergedDf['mean_sub_indv']) > mergedDf['thickstd'] *2].iterrows()):
                    if (sigNum+1) % 2 == 0:
                        sign = 1
                        diff = 0.5
                    else:
                        sign = -1
                        diff = 0

                    if row[1].mean_sub_indv < 0:
                        col = 'green'
                        sign = 1
                    else:
                        col = 'red'

                    textLoc_y = row[1].thickness_y + (sign*1) + diff

                    if row[1].roi not in sig_diff_region:
                        ax.annotate(row[1].roi,
                                    xy = (row[0], row[1].thickness_y), 
                                    xycoords='data',
                                    xytext = (row[0], textLoc_y),
                                    textcoords='data',
                                    arrowprops = dict(facecolor=col, 
                                                      shrinkB=5,
                                                      alpha=0.5), 
                                    horizontalalignment='center', 
                                    fontsize=15)
                        sig_diff_region.append(row[1].roi)

    # axis settings
    for snum, side in enumerate(['lh', 'rh']):
        ax = axes[snum]
        ax.patch.set_facecolor('white')
        # Graph settings
        ax.set_ylim(1, 5)
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

        regionToHighlight = roiOrder[1::2]
        xCoords = [np.where(roiOrder_array==x)[0] for x in regionToHighlight]

        for x in xCoords:
            ax.axvspan(x[0], x[-1], alpha=0.5)

        startNum = 0
        for region in roiOrder:
            x_starting_point = startNum
            startNum = startNum + len(roiDict[region])

            ax.text((x_starting_point-.5 + startNum-.5)/2, 1.2,
                    region,
                    horizontalalignment='center',
                    alpha=.4,
                    fontsize=15)

        ax.set_xticklabels(infodf.roi)
        labels = ax.get_xticklabels()
        plt.setp(labels, rotation=30)
        plt.tight_layout(pad=7, w_pad=3, h_pad=0.2)

    #fig.show()
    fname = 'thickness_summary.png'
    fig.savefig(fname)
    print('feh %s' %join(os.getcwd(), fname))

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

            #print output

            try:
                thickness = re.search('thickness\s+=\s+(\S+)\s+mm\s+\S+\s+(\S+)', output).group(1,2)
                numvert = re.search('number of vertices\s+=\s+(\S+)', output).group(1)
                surfarea = re.search('total surface area\s+=\s+(\S+)', output).group(1)
                grayvol = re.search('total gray matter volume\s+=\s+(\S+)', output).group(1)
                meancurv = re.search('average integrated rectified mean curvature\s+=\s+(\S+)', output).group(1)
                gauscurv = re.search('average integrated rectified Gaussian curvature\s+=\s+(\S+)', output).group(1)
                foldind = re.search('folding index\s+=\s+(\S+)', output).group(1)
                curvind = re.search('intrinsic curvature index\s+=\s+(\S+)', output).group(1)
            except:
                sys.exit('{0} : re.search with mris_antomical_stats not working. Check row data or FS environment'.format(fsDir))
                

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


        #print re.sub('\s+',' ',command)
        os.popen(re.sub('\s+',' ',command)).read()

def asegstats2table(fsDir):
    os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
    os.environ["SUBJECTS_DIR"] = dirname(fsDir)

    output_text = join(fsDir, 'tmp', 'subcortex_table.txt')
    command = 'asegstats2table \
            --subjects {dirName} \
            -t {output_text}'.format(dirName = basename(fsDir),
                                     output_text = output_text)

    if not isfile(output_text):
        os.popen(command).read()
    df = pd.read_table(output_text).T.reset_index()
    df = df.drop(0)
    df.columns = ['roi', 'volume']
    df['region'] = 'subcortex'
    df['volume'] = df['volume'].astype('float')

    return df

def aparcstats2table(fsDir, parc):
    os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
    os.environ["SUBJECTS_DIR"] = dirname(fsDir)

    allDf = pd.DataFrame()
    measures = 'volume', 'thickness', 'thicknessstd'

    for side in 'lh', 'rh':
        for measure in measures:
            output_text = join(fsDir, 'tmp', '{0}_{1}_table.txt'.format(side, measure))
            command = 'aparcstats2table \
                    --hemi {side} \
                    --subjects {dirName} \
                    --parc {parc} \
                    --meas {measure} \
                    -t {output_text}'.format(side = side,
                                             dirName = basename(fsDir),
                                             parc = parc,
                                             measure = measure,
                                             output_text = output_text)

            if not isfile(output_text):
                os.popen(command).read()

            df = pd.read_table(output_text).T
            df = df.drop(df.index[0])
            df.columns = ['value']
            df['roi'] = df.index.str.extract('[rl]h_(\S+)_{0}'.format(measure), expand=False)
            df['side'] = side
            df['measure'] = measure
            df['region'] = df.roi.apply(getRegion)

            df = df.reset_index()
            allDf = pd.concat([allDf, df[['roi', 'region', 'side', 'measure', 'value']]])

    allDf = pd.pivot_table(allDf, 
                           values='value', 
                           index=['side','roi', 'region'], 
                           columns='measure', 
                           aggfunc=np.sum)
    allDf = allDf.reset_index()
    for measure in measures:
        allDf[measure] = allDf[measure].astype('float')

    return allDf

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

    Returns infoDf
        ,subroi,numvert,surfarea,grayvol,thickavg,thickstd,meancurv,gauscurv,foldind,curvind,side
        0,lh_bankssts,1803.0,1207.0,2743.0,2.301,0.439,0.107,0.019,13.0,1.5,lh
    '''

    # cortical regions as a dictionary
    roiDict = get_cortical_rois_detailed()

    if not isfile(join(fsDir,'tmp','thick_kev_detailed_new.csv')):
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
            eg) {codeName} -i YB_FS KCHO_FS
            eg) {codeName} -i YB_FS KCHO_FS -n Yoobin Kangik
            eg) {codeName} -i YB_FS KCHO_FS -n Yoobin Kangik -g f m
            eg) {codeName} -i YB_FS KCHO_FS -n Yoobin Kangik -g f m -c red blue
            eg) {codeName} -i YB_FS KCHO_FS -n Yoobin Kangik -g f m -c red blue -nb
            '''.format(codeName=basename(__file__))))

    parser.add_argument(
        '-i', '--inputDirs',
        help='One or more locations of freesurfer directory',
        nargs='+',
        default=False)

    parser.add_argument(
        '-n', '--nameList',
        help='List of legend names for each freesurfer summary',
        nargs='+',
        default=False)

    parser.add_argument(
        '-a', '--ageList',
        help='List of age for each freesurfer summary',
        nargs='+',
        default=False)

    parser.add_argument(
        '-ar', '--ageRange',
        help='List of age for each freesurfer summary',
        default=3)

    parser.add_argument(
        '-g', '--genderList',
        help='List of genders for each freesurfer summary',
        nargs='+',
        default=False)

    parser.add_argument(
        '-c', '--colorList',
        help='List of colors for each freesurfer summary',
        nargs='+',
        default=False)

    parser.add_argument(
        '-nb', '--nobackground',
        help='Removes mean cortical thickness graph of the healthy controls in the background',
        action='store_true',
        default=False)

    parser.add_argument(
        '-m', '--makeMean',
        help='Make mean df from the given fsDirs',
        action='store_true',
        default=False)

    args = parser.parse_args()

    # if lengths of the lists do not match
    # edit this later
    for i in args.genderList, args.nameList, args.colorList:
        if i:
            if not len(args.inputDirs) == len(i):
                print('The number freesurfer directories : {0}'.format(len(args.inputDirs)))
                print('The number items given : {0}'.format(len(i)))
                sys.exit('Please make sure there are the same number of items in each list')

    # Make mean
    if args.makeMean:
        makeMean(args.inputDirs)
    else:
        freesurferSummary(args.inputDirs,
                         args.nameList,
                         args.ageList,
                         args.genderList,
                         args.ageRange,
                         args.colorList,
                         args.nobackground)
