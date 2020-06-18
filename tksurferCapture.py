#!/usr/bin/env python3
import shutil
import os
from os.path import join, isfile, basename, dirname
from pathlib import Path
import re
from nilearn import datasets
from nilearn import plotting
import nipype.interfaces.freesurfer as fs          # fsl
import pandas as pd
import argparse
import textwrap
import sys
from plumbum import local
pd.set_option('max_rows', 5000)

# def labelValueSwap(labelLoc:str, newNum:int) -> None:
    # '''Replace zeros with newNum'''
    # #print labelLoc

    # # read lines from labelLoc
    # with open(labelLoc, 'r') as f:
        # lines = f.read()

    # newLines = lines.replace('0.0000000000', str(newNum))
    
    # newName = os.path.dirname(labelLoc)+'/' +\
            # os.path.basename(labelLoc)[:-6] + \
            # '_new.label'

    # #print newName
    # with open(newName, 'w') as newF:
        # newF.write(newLines)


def tclWrite(thicknessAsc:str, location:str) -> None:
    '''Write tcl script that saves screenshots from the thicknessAsc file'''

    directory = Path(location).parent
    toWrite = f'''set curv {thicknessAsc}
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
exit 0'''

    with open(location, 'w') as f:
        f.write(toWrite)

def ctedit(ctab:str, labelName:str, newNum:int) -> None:
    '''Change color for the given label in the color tab file'''

    # read color table file
    with open(ctab, 'r') as f:
        lines = f.readlines()

    newLines = []
    # for each line in the color table
    for line in lines:
        if labelName in line:
            newLine = re.sub(r'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', 
                             f'60  60  60  {int(newNum)}',
                             line)
            newLines.append(newLine)
        else:
            newLines.append(line)

    # wirte the new color file
    with open(ctab, 'w') as f:
        for line in newLines:
            f.write(line)


def makeBrainPic(freesurfer_dir:str):
    '''Save brain images into files'''
    freesurfer_dir = Paht(freesurfer_dir)

    for side in ['lh','rh']:
        tcl_script = freesurfer_dir / 'tmp' / side + '_tksurfer.tcl'
        #annotFile = os.path.join(freesurfer_dir,
                #'mri/aparc+aseg.mgz')
        # print freesurfer_dir

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

        command = f'tksurfer {freesurfer_dir.name} {side} inflated \
                -tcl {tcl_script}'
        print(command)
        os.environ["FREESURFER_HOME"] = local.path(os.getenv('FREESURFER_HOME'))
        os.environ["SUBJECTS_DIR"] = freesurfer_dir.parent

        # print(shots.cmdline)
        #print "export SUBJECTS_DIR={0}".format(os.path.dirname(freesurfer_dir))
        # os.popen(shots.cmdline).read()
        os.popen(command).read()

        # createdList = [side+'_'+x+'.tif' for x in ["inferior", "lateral", "medial"]]
        # for img in createdList:
            # #print img
            # if 'baseline' in freesurfer_dir:
                # folderName = os.path.dirname(freesurfer_dir).split('/baseline')[0]
            # else:
                # folderName = os.path.dirname(freesurfer_dir)

            # shutil.move(img.split('_')[1], 
                    # os.path.join('/ccnc/mri_team/',
                        # os.path.basename(folderName) + '_' + os.path.basename(img)))


def mergeLabel(freesurfer_dir:str, labelNames:list) -> None:
    '''Merge labels into single file using Freesurfer mri_mergelabels'''
    freesurfer_dir = Path(freesurfer_dir)
    out_dir = freesurfer_dir / 'tmp'

    in_label = ' '.join(
            [f'-i {out_dir}/{side}.{x}_new.label' for x in labelNames]
            )
    for side in 'lh', 'rh':
        out_label = f'{out_dir}/{side}_all.label'
        command = f'mri_mergelabels {in_label} -o {out_label} 2>/dev/null'

        os.popen(command).read()


def cleanMean(meanCSV:str, indCSV:str) -> pd.DataFrame:
    '''Merge mean csv with individual csv file'''
    # read mean CSV
    meanDf = pd.read_csv(meanCSV, index_col=0)
    meanDf = meanDf.groupby(['roi', 'side']).mean().reset_index()

    df = pd.read_csv(indCSV, index_col=0)
    df['roi'] = df.subroi.str[3:]
    df['side'] = df.subroi.str[:2]

    mergedDf = pd.merge(meanDf, df,
                        on=['roi', 'side'],
                        how='inner')

    # individual subject thickness - mean thickness
    mergedDf['mean_sub_indv'] = mergedDf.thickavg_y - mergedDf.thickavg_x

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

def makeAsc(freesurferLoc:str, output_dir:str) -> None:
    '''Converts thickness and white matter surface map to a ascii file'''
    freesurferLoc = Path(freesurferLoc)
    output_dir = Path(output_dir)

    for side in 'lh', 'rh':
        thickSurf = freesurferLoc  / 'surf' / f'{side}.thickness'
        whiteSurf = freesurferLoc  / 'surf' / f'{side}.white'
        # output = freesurferLoc  / 'tmp' / f'{side}.thickness.asc'
        output = output_dir / f'{side}.thickness.asc'

        command = f'mris_convert -c {thickSurf} {whiteSurf} {output}'
        out = os.popen(command).read()


def main(fsDir:str, csv:str, mean_csv_loc:str):
    '''Make brain slice screenshots using tksurfer'''
    fsDir = Path(fsDir)

    # mean_csv_loc = '/ccnc_bin/meanThickness/detailed_mean_2017_06_16.csv'

    # get mean df + individual df
    mergedDf = cleanMean(mean_csv_loc, csv)

    #standard_ctab = '/ccnc_bin/meanThickness/standard_ctab.txt'
    fs_home = local.path(os.getenv('FREESURFER_HOME'))
    standard_ctab = f'{fs_home}/FreeSurferColorLUT.txt'

    subject_ctab = fsDir / 'ctab.txt'
    shutil.copyfile(standard_ctab, subject_ctab)

    # converts surface files into ascii files
    makeAsc(fsDir, fsDir / 'tmp')

    # for each side
    for side in 'lh', 'rh':
        mergedDfSide = mergedDf.groupby('side').get_group(side)

        # load thickness ascii file
        ascii_file = fsDir / 'tmp' / f'{side}.thickness.asc'
        ascii_file_new = fsDir / 'tmp' / f'new.{side}.thickness.asc'

        ascii_df = pd.read_csv(ascii_file, sep='\s+')
        ascii_df.columns = ['vertexNum', 'x','y','z','value']
        ascii_df.set_index('vertexNum', inplace=True)
        ascii_df['value'] = 0

        # for each labels
        for index, row in mergedDfSide.iterrows():
            # get locations
            roi_name = f'{row.side}.{row.roi}.label'
            roi_loc = fsDir / 'tmp' / roi_name

            # get label df
            df = pd.read_csv(roi_loc, skiprows=[0,1], sep='\s+')
            df.columns = ['vertexNum', 'x', 'y', 'z', 'value']

            # here replacing value
            df['value'] = row.mean_sub_indv

            df.set_index('vertexNum', inplace=True)

            #update asciiDf
            ascii_df.update(df)

            df.roi = row.roi

            #labelValueSwap(roiLoc, row[1].mean_sub_indv)
        #asciiDf.index = asciiDf.index.apply(lambda x: x.zfill(15))
        ascii_df.to_csv(ascii_file_new,sep=' ', header=False)
        tcl_loc = fsDir / 'tmp' / f'{side}_tksurfer.tcl'
        tclWrite(asciiFileNew, tcl_loc)

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
        '-m', '--mean_csv',
        help='Location of mean csv file',
        default=os.getcwd())

    parser.add_argument(
        '-c', '--csv',
        help='Output csv file from freesurfer_summary.py',
        default=join(os.getcwd(), 'tmp/thick_kev_detailed_new.csv'))

    args = parser.parse_args()
    os.environ["FREESURFER_HOME"] = local.path(os.getenv('FREESURFER_HOME'))
    os.environ["SUBJECTS_DIR"] = Path(args.fsDir).parent

    if Path(args.mean_csv).is_file():
        pass
    # else:
        # fsaverage_dir = Path(local.path(os.getenv('FREESURFER_HOME'))) \
                # / 'subjects/fsaverage'
        # makeAsc(fsaverage_dir, Path(args.mean_csv).parent)
        # ascii_file = fsDir / 'tmp' / f'{side}.thickness.asc'

        # for 

    if Path(args.csv).is_file():
        main(args.fsDir, args.csv, args.mean_csv)
    else:
        print(args.csv)
        sys.exit('CSV is missing - Please learn freesurfer_summary.py')

