#!/ccnc/anaconda2/bin/python
from __future__ import division
import os
import re
from os.path import join, basename, dirname, isfile
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
import argparse
import itertools
import textwrap
from progressbar import ProgressBar
import time
#import tksurferCapture

__author__ = 'kcho'
plt.style.use('ggplot')

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

def freesurferSummary(args):
    '''
    Output freesurfer summary
    '''
    # Read CCNC healthy control information
    meanDfLoc = '/ccnc_bin/meanThickness/detailed_mean_2017_06_16.csv'
    meanDf, meanDf_name = getGroupMeanInfo(meanDfLoc)

    # Collect infoDfs and names
    fsInfoDf = []
    fsNames = []
    for argsFsDir in args.inputDirs:
        argsSubjNames = raw_input('{0} Subject initial :'.format(argsFsDir))
        fsNames.append(argsSubjNames)
        fsInfoDf.append(collectStats(os.path.abspath(argsFsDir)))

    # Add CCNC HCs informat
    fsInfoDf.append(meanDf)
    fsNames.append(meanDf_name)

    # Make line plots of cortical thickness for each hemisphere
    draw_thickness_list(fsInfoDf, fsNames)

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

def draw_thickness_list(infoDfList, nameList):
    # graph order from left
    roiOrder = ['LPFC', 'OFC', 'MPFC', 'LTC', 'MTC', 'SMC', 'PC', 'OCC']

    fig, axes = plt.subplots(nrows=2, figsize=(22,12), facecolor='white')
    fig.suptitle("Cortical thickness in all regions", fontsize=20)

    # color definition
    color=iter(cm.rainbow(np.linspace(0,1,len(infoDfList))))

    # mean infoDf
    meanDf = infoDfList[-1]
    roiList = meanDf.roi.unique()

    for infoDf, subjName in zip(infoDfList, nameList):
        infoDf_gb = infoDf.groupby('side')
        c = next(color)
        for snum, side in enumerate(['lh', 'rh']):
            infodf = infoDf_gb.get_group(side)

            # Reorder Dfs
            infodf = infodf.sort_values(['roi','side'])
            infodf = reorder_df(infodf, 'region', roiOrder)

            ax = axes[snum]
            if subjName=='CCNC_mean':
                #ax.plot(infodf.thickavg, '--', c=c, label=subjName)
                eb = ax.errorbar(range(len(roiList)),
                                  infodf.thickavg,
                                  infodf.thickstd,
                                  marker='^',
                                  label=subjName, 
                                  color='b',
                                  ecolor='b',
                                  alpha=0.7)
                plotline, caplines, barlinecols = eb
                barlinecols[0].set_linestyle('--')
                plotline.set_linestyle('--')
            else:
                ax.plot(infodf.thickavg, c=c, label=subjName)

                ### annotation
                mergedDf = pd.merge(meanDf,
                                    infodf,
                                    on=['roi','side','region'],
                                    how='inner')
                mergedDf['mean_sub_indv'] = mergedDf.thickavg_x - mergedDf.thickavg_y

                # Reorder Dfs
                mergedDf = mergedDf.sort_values(['roi','side'])
                mergedDf = reorder_df(mergedDf, 'region', roiOrder)

                for sigNum, row in enumerate(mergedDf[abs(mergedDf['mean_sub_indv']) > 0.5].iterrows()):

                    if (sigNum+1) % 2 == 0:
                        sign = 1
                    else:
                        sign = -1
                    textLoc_y = row[1].thickavg_y + (sign*0.5)

                    if row[1].mean_sub_indv < 0:
                        col = 'green'
                    else:
                        col = 'red'

                    ax.annotate(row[1].roi,
                                xy = (row[0], row[1].thickavg_y), 
                                xycoords='data',
                                xytext = (row[0], textLoc_y),
                                textcoords='data',
                                arrowprops = dict(facecolor=col, 
                                                  shrinkB=5,
                                                  alpha=0.5), 
                                horizontalalignment=side, 
                                fontsize=15)

    # axis settings
    for snum, side in enumerate(['lh', 'rh']):
        ax = axes[snum]
        ax.patch.set_facecolor('white')
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
    print(fname)

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

            #print output

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
            eg) {codeName} --input {in_put} --output {output}
            '''.format(codeName=basename(__file__),
                       in_put = 'subjectLoc',
                       output = 'outLoc')))

    parser.add_argument(
        '-i', '--inputDirs',
        help='One or more locations of freesurfer directory',
        nargs='+',
        default=False)

    args = parser.parse_args()

    freesurferSummary(args)
