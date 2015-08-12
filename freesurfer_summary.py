#!/ccnc_bin/venv/bin/python
__author__ = 'kcho'

import re
import os
import pandas as pd
import matplotlib
matplotlib.use('GTK')
import matplotlib.pyplot as plt
import argparse
import textwrap



def main(subject_loc = '/Users/kcho/T1', locations=['/Users/kcho/T1','/Users/kcho/T1'], roi_list = ['ctx_lh_G_cuneus']):
    ##########################################################
    # Freesurfer setting
    ##########################################################
    os.environ["FREESURFER_HOME"] = '/Applications/freesurfer'
    # where is the freesurfer directory
    os.environ["SUBJECTS_DIR"] = '{0}'.format(subject_loc)

    freesurfer_dir = get_freesurfer_loc(subject_loc)

    ##########################################################
    # Get roi dict : 8 cortex each hemisphere
    ##########################################################
    roiDict = get_cortical_rois()

    ##########################################################
    # annotation2label --> merge labels --> freesurfer/tmp
    ##########################################################
    if not os.path.isfile(os.path.join(freesurfer_dir,'tmp','thick_kev.csv')):
        makeLabel(freesurfer_dir)
        mergeLabel(freesurfer_dir,roiDict)

        ##########################################################
        # thicknessDict[side_cortex] = (thickness, std) in mm
        ##########################################################
        thicknessDict = getThickness(freesurfer_dir,roiDict)
        thicknessDf = dictWithTuple2df(thicknessDict)
        #print thicknessDf
        thicknessDf.to_csv(os.path.join(freesurfer_dir,'tmp','thick_kev.csv'))
    else:
        thicknessDf = pd.read_csv(os.path.join(freesurfer_dir,'tmp','thick_kev.csv'))
    draw_thickness(thicknessDf,os.path.basename(subject_loc))


    volumeDf = openStatsTable(freesurfer_dir)
    volumeDf['name'] = volumeDf.side + '_' + volumeDf.ROI

    volumeDf = getSummary(volumeDf,roiDict)
    #print volumeDf


    # graph
    #draw_graph(volumeDf)

    # # collect stats
    # if len(locations) > 1:
    #     meanDf = collectStats(locations)



    #freesurfer_table = openTable('/Applications/freesurfer/FreeSurferColorLUT.txt')
    #print freesurfer_table

    #for roi in roi_list:
    #    print freesurfer_table[roi]


def draw_thickness(thicknessDf,subjName):
    thicknessDf['roi'] = thicknessDf.subroi.str[3:]
    gb = thicknessDf.groupby('roi')
    thicknessDf = pd.concat([gb.get_group('LPFC'),
                        gb.get_group('OFC'),
                        gb.get_group('MPFC'),
                        gb.get_group('LTC'),
                        gb.get_group('MTC'),
                        gb.get_group('SMC'),
                        gb.get_group('PC'),
                        gb.get_group('OCC')])

    gb = thicknessDf.groupby('side')
    label = thicknessDf.subroi.str[3:].unique()

    fig = plt.figure(figsize=(10,5))
    fig.suptitle("Cortical thickness in eight regions", fontsize=20)
    #plt.ylabel('Cortical thickness', fontsize=16)
    #plt.xticks(range(len(label)), label)

    lh_g = plt.subplot2grid((2,2),(0, 0), rowspan=2)
    rh_g = plt.subplot2grid((2,2),(0, 1), rowspan=2)
    #ax1 = fig.add_subplot(211)
    #ax2 = fig.add_subplot(212)
    lh_g.plot(gb.get_group('lh')['thickness'],'r',label=subjName)
    lh_g.set_xticklabels(label)
    lh_g.set_xlabel('Left', fontsize=16)
    lh_g.set_ylabel('Cortical thickness in mm', fontsize=16)
    lh_g.set_ylim(1.8, 3.2)
    lh_g.legend()

    rh_g.plot(gb.get_group('rh')['thickness'],'b',label=subjName)
    rh_g.set_xticklabels(label)
    rh_g.set_xlabel('Right', fontsize=16)
    rh_g.set_ylim(1.8, 3.2)
    rh_g.legend()

    plt.show()


def dictWithTuple2df(thicknessDict):
    df = pd.DataFrame.from_dict(thicknessDict)
    df = df.stack().reset_index()
    df = df[df.level_0==0].merge(df[df.level_0==1], on='level_1', how='inner')
    df.columns = ['__','subroi','thickness','_','std']
    df['side'] = df['subroi'].str[:2]

    return df[['subroi','thickness','side','std']]




def getThickness(freesurfer_dir,roiDict):
    thicknessDict={}
    for side in ['lh','rh']:
        for cortex, rois in roiDict.iteritems():
            command = 'mris_anatomical_stats \
            -l {loc}/{side}_{cortex} FREESURFER {side} 2>/dev/null'.format(
                loc=os.path.join(freesurfer_dir,'tmp'),
                side=side,
                cortex=cortex
            )
            output=os.popen(re.sub('\s+',' ',command)).read()
            thickness = re.search('thickness\s+=\s+(\S+)\s+mm\s+\S+\s+(\S+)', output).group(1,2)
            thickness = tuple([float(x) for x in thickness])
            thicknessDict[side+'_'+cortex] = thickness
            print thickness
            os.remove('{loc}/{side}_{cortex}'.format(
                loc=os.path.join(freesurfer_dir,'tmp'),
                side=side,
                cortex=cortex
            ))
    return thicknessDict


def mergeLabel(freesurfer_dir, roiDict):
    for side in ['lh','rh']:
        for cortex, rois in roiDict.iteritems():
            command = 'mri_mergelabels {inLabel} -o {outLabel} 2>/dev/null'.format(
                inLabel = ' '.join(['-i '+os.path.join(freesurfer_dir,'tmp',side+'.'+x+'.label') for x in rois]),
                outLabel = os.path.join(freesurfer_dir,'tmp',side+'_'+cortex))
            os.popen(command).read()

            for roi in [os.path.join(freesurfer_dir,'tmp',side+'.'+x+'.label') for x in rois]:
                os.remove(roi)


def makeLabel(freesurfer_dir):
    for side in ['lh','rh']:
        command = 'mri_annotation2label \
            --subject {basename} \
            --hemi {side} --outdir {outDir} 2>/dev/null'.format(basename=os.path.basename(freesurfer_dir),
                                                    side=side,
                                                    outDir=os.path.join(freesurfer_dir,'tmp'))
        os.popen(re.sub('\s+',' ',command)).read()


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

    print rh_volume_sums
    plt.plot(lh_volume_sums.values(),'r')
    plt.plot(rh_volume_sums.values(),'b')
    plt.xticks(range(len(cortexList)), cortexList)
    plt.show()
    # plt.plot(gb.get_group('lh')['Volume'],'r-')
    # plt.plot(gb.get_group('rh')['Volume'],'b-')
    # plt.xticks(range(len(gb.get_group('lh')['ROI'])),gb.get_group('lh')['ROI'])
    # plt.show()

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


def get_cortical_rois():
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
        '-i', '--inputDir',
        help='Subject location',
        default=os.getcwd())

    #parser.add_argument(
        #'-o', '--output',
        #help='Output',
        #default=os.getcwd())

    parser.add_argument(
        '-l', '--locations',
        help='subject inputs in python list format',
        default=[x for x in os.listdir(os.getcwd()) if os.path.isdir(x)])


    parser.add_argument(
        '-r', '--rois',
        help='roi inputs in python list format',
        default="ctx_lh_G_cuneus")

    args = parser.parse_args()

    main(args.inputDir, args.locations, args.rois)
