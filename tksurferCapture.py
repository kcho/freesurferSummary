#!/ccnc/anaconda2/bin/python
import shutil
import os
from os.path import join, isfile, basename, dirname
import re
from nilearn import datasets
from nilearn import plotting
import nipype.interfaces.freesurfer as fs          # fsl
import pandas as pd
import argparse
import textwrap
import sys
pd.set_option('max_rows', 5000)

def labelValueSwap(labelLoc, newNum):
    #print labelLoc
    with open(labelLoc, 'r') as f:
        lines = f.read()

    newLines = lines.replace('0.0000000000', str(newNum))
    
    newName = os.path.dirname(labelLoc)+'/' +\
            os.path.basename(labelLoc)[:-6] + \
            '_new.label'
    #print newName
    with open(newName, 'w') as newF:
        newF.write(newLines)


def tclWrite(location, thicknessAsc):

    toWrite = '''set curv {thicknessAsc}
set tiff_directory {directory}
read_binary_curv
curv_to_val
set overlayflag 1
#set fthresh .001
#set fmid .2
sclv_set_current_threshold_from_percentile .92 .93 .99
set colscalebarflag 1
set curvflag 0
UpdateAndRedraw
redraw
make_lateral_view
save_tiff lateral.tif

rotate_brain_x 90
redraw
save_tiff inferior.tif

make_lateral_view
rotate_brain_y 180 
redraw
save_tiff medial.tif

# will cause FS to exit
exit 0'''.format(thicknessAsc = thicknessAsc, 
        directory = os.path.dirname(thicknessAsc),)

    with open(location, 'w') as f:
        f.write(toWrite)

def ctedit(ctab, labelName, newNum):
    with open(ctab, 'r') as f:
        lines = f.readlines()

    newLines = []
    for line in lines:
        if labelName in line:
            #try:
                #print re.search(line,'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)').group(0)
                #print re.search('(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line).group(0)
                #newLine = re.sub('(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', 
                        #['haho','hoho','..','kk'],
                        #line)
            newLine = re.sub('(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', 
                    '60  60  60  {value}'.format(value=int(newNum)), 
                    line)
            #newLine = re.sub('(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', 
                    #'250 300 250 {value}'.format(value=newNum), 
                    #line)
            newLines.append(newLine)
            #except:
                #pass
        else:
            newLines.append(line)
    with open(ctab, 'w') as f:
        for line in newLines:
            f.write(line)


def makeBrainPic(freesurfer_dir):
    for side in ['lh','rh']:
        tcl_script = os.path.join(freesurfer_dir,
                'tmp',side+'_tksurfer.tcl')
        #annotFile = os.path.join(freesurfer_dir,
                #'mri/aparc+aseg.mgz')
        print freesurfer_dir


        #print os.path.dirname(freesurfer_dir)
        # shots = fs.SurfaceSnapshots(
                # subject_id=os.path.basename(freesurfer_dir),
                # subjects_dir=os.path.dirname(freesurfer_dir),
                # hemi=side,
                # tcl_script=tcl_script,
                # #annot_file = annotFile, 
                # #overlay_range = (-4,4),
                # #screenshot_stem = freesurfer_dir+'_'+side,
                # six_images = True,
                # surface="inflated")
        command = 'tksurfer {subjid} {side} inflated \
                -tcl {tcl}'.format(
                    subjid=os.path.basename(freesurfer_dir),
                    side=side,
                    tcl=tcl_script)

        os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
        os.environ["SUBJECTS_DIR"] = os.path.dirname(freesurfer_dir)

        # print(shots.cmdline)
        #print "export SUBJECTS_DIR={0}".format(os.path.dirname(freesurfer_dir))
        # os.popen(shots.cmdline).read()
        os.popen(command).read()
        #print '*'*80
        #try:
            ##res = shots.run()
            #pass

        #except:
            #pass

        createdList = [side+'_'+x+'.tif' for x in ["inferior", "lateral", "medial"]]
        for img in createdList:
            #print img
            if 'baseline' in freesurfer_dir:
                folderName = os.path.dirname(freesurfer_dir).split('/baseline')[0]
            else:
                folderName = os.path.dirname(freesurfer_dir)

            shutil.move(img.split('_')[1], 
                    os.path.join('/ccnc/mri_team/',
                        os.path.basename(folderName) + '_' + os.path.basename(img)))




def mergeLabel(freesurfer_dir, labelNames):
    for side in ['lh','rh']:
        #command = 'mri_mergelabels {inLabel} -o {outLabel} 2>/dev/null'.format(
        command = 'mri_mergelabels {inLabel} -o {outLabel}'.format(
            inLabel = ' '.join(['-i '+os.path.join(freesurfer_dir,'tmp',side+'.'+x+'_new.label') for x in labelNames]),
            outLabel = os.path.join(freesurfer_dir,'tmp',side+'_all.label'))

        #print command
        os.popen(command).read()


def cleanMean(meanCSV, indCSV):
    # read mean CSV
    meanDf = pd.read_csv(meanCSV, index_col=0)
    meanDf = meanDf.groupby(['roi','side']).mean().reset_index()

    df = pd.read_csv(indCSV, index_col=0)
    df['roi'] = df.subroi.str[3:]
    df['side'] = df.subroi.str[:2]

    mergedDf = pd.merge(meanDf,
                        df,
                        on=['roi','side'],
                        how='inner')

    # individual subject thickness - mean thickness
    mergedDf['mean_sub_indv'] = (mergedDf.thickavg_y - mergedDf.thickavg_x)
    mergedDf['mean_sub_indv_cov'] =  mergedDf['mean_sub_indv']
    return mergedDf

#def rescale(npArray, new_min = 70, new_max = 100):
    #values = npArray.tolist()
    #print values
    #output = []
    #old_min, old_max = min(values), max(values)

    #for v in values:
        #new_v = (new_max - new_min) / (old_max - old_min) * (v - old_min) + new_min
        #output.append(new_v)

    #return output

def makeAsc(freesurferLoc):
    '''
    Converts thickness and white matter surface map to a ascii file
    '''
    for side in ['lh','rh']:
        thickSurf = os.path.join(freesurferLoc,
                'surf',side+'.thickness')
        whiteSurf = os.path.join(freesurferLoc,
                'surf',side+'.white')
        output = os.path.join(freesurferLoc,
                'tmp', side+'.thickness.asc')

        toAsc = fs.MRIsConvert(
                scalarcurv_file = thickSurf,
                in_file = whiteSurf,
                out_datatype = 'ico',
                #out_file = output,
                )
        command = 'mris_convert -c {thickF} {whiteF} {outF}'.format(
                thickF = thickSurf,
                whiteF = whiteSurf,
                outF = output)

        out = os.popen(command).read()

def main(fsDir,csv):
    '''
    Make brain slice captures using tksurfer
    '''

    mean_csv_loc = '/ccnc_bin/meanThickness/detailed_mean_2017_06_16.csv'

    # returns mean df + individual df
    mergedDf = cleanMean(mean_csv_loc, csv)

    #standard_ctab = '/ccnc_bin/meanThickness/standard_ctab.txt'
    standard_ctab = '/usr/local/freesurfer/FreeSurferColorLUT.txt'
    subject_ctab = join(fsDir,'ctab.txt')
    shutil.copyfile(standard_ctab, subject_ctab)

    # converts surface files into ascii files
    makeAsc(fsDir)

    # for each side
    for side in ['lh','rh']:
        mergedDfSide = mergedDf.groupby('side').get_group(side)

        # load thickness ascii file
        asciiFile = join('{0}/tmp/{1}.thickness.asc'.format(fsDir, side))
        asciiFileNew = re.sub('{0}.thickness'.format(side),
                              'new.{0}.thickness'.format(side),
                             asciiFile)

        asciiDf = pd.read_csv(asciiFile, sep='\s+')
        asciiDf.columns = ['vertexNum', 'x','y','z','value']
        asciiDf.set_index('vertexNum', inplace=True)
        asciiDf['value'] = 0

        # for each labels
        for row in mergedDfSide.iterrows():
            # get locations
            roiName = row[1].side + '.' + row[1].roi + '.label'
            roiLoc = join(fsDir, 'tmp',roiName)

            # get newNumber
            newNum = row[1].mean_sub_indv

            # get label df
            df = pd.read_csv(roiLoc, skiprows=[0,1], sep='\s+')
            df.columns = ['vertexNum', 'x','y','z','value']
            df['value'] = newNum
            #print df[df['value'] != 0]

            df.set_index('vertexNum', inplace=True)

            #update asciiDf
            asciiDf.update(df)

            df.roi = row[1].roi

            #labelValueSwap(roiLoc, row[1].mean_sub_indv)
        #asciiDf.index = asciiDf.index.apply(lambda x: x.zfill(15))
        asciiDf.to_csv(asciiFileNew,sep=' ', header=False)
        tclWrite(join('{0}/tmp/{1}_tksurfer.tcl'.format(fsDir,side)),
                asciiFileNew)

    #mergeLabel(fsDir, mergedDf.roi.unique())
    makeBrainPic(fsDir)



#if __name__=='__main__':
    #main(sys.argv[1], sys.argv[2])

    ##main('/Volumes/promise/CCNC_MRI_3T/NOR/NOR103_SHS/baseline/FREESURFER',
            ##'/Volumes/promise/CCNC_MRI_3T/NOR/NOR103_SHS/baseline/FREESURFER/tmp/thick_kev_detailed.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=textwrap.dedent('''\
            {codeName} : 
            ========================================
            eg) {codeName} --fsDir {inputLoc} 
                {codeName} --fsDir {inputLoc} --csv ~/thick_kev_detailed.csv
            '''.format(codeName=basename(__file__),
                       inputLoc = 'subjectLoc')))

    parser.add_argument(
        '-f', '--fsDir',
        help='Freesurfer directory location',
        default=os.getcwd())

    parser.add_argument(
        '-c', '--csv',
        help='Output csv file from freesurfer_summary.py',
        default=join(os.getcwd(), 'tmp/thick_kev_detailed_new.csv'))

    args = parser.parse_args()
    os.environ["FREESURFER_HOME"] = '/usr/local/freesurfer'
    os.environ["SUBJECTS_DIR"] = dirname(dirname(args.fsDir))

    if isfile(args.csv):
        main(args.fsDir, args.csv)
    else:
        print(args.csv)
        sys.exit('CSV is missing - Please learn freesurfer_summary.py')

