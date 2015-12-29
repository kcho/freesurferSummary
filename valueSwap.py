import shutil
import os
import re
from nilearn import datasets
from nilearn import plotting
import nipype.interfaces.freesurfer as fs          # fsl
import pandas as pd

def labelValueSwap(labelLoc, newNum):
    print labelLoc
    with open(labelLoc, 'r') as f:
        lines = f.read()

    newLines = lines.replace('0.0000000000', str(newNum))
    
    newName = os.path.dirname(labelLoc)+'/' +\
            os.path.basename(labelLoc)[:-6] + \
            '_new.label'
    print newName
    with open(newName, 'w') as newF:
        newF.write(newLines)


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


def makeBrainPic(freesurfer_dir, label, ctab):
    bgBrain = os.path.join(freesurfer_dir,'mri/brain.mgz')
    bgBrainNii = os.path.dirname(bgBrain) + 'brain.nii.gz'

    #converter = fs.MRIConvert(backgroundin_file = bgBrain,
            #out_type='niigz',
            #out_file=bgBrainNii)
    #converter.run()

    for side in ['lh','rh']:
        labelName = os.path.join(freesurfer_dir,
                'tmp',
                side+'_all.label')
        #annotFile = os.path.join(freesurfer_dir,
                #'mri/aparc+aseg.mgz')
        annotFile = os.path.join(freesurfer_dir,
                'label/{side}.aparc.annot'.format(side=side))


        shots = fs.SurfaceSnapshots(
                subject_id=os.path.basename(freesurfer_dir),
                subjects_dir=os.path.dirname(freesurfer_dir),
                hemi=side,
                label_file = labelName,
                #label_file = annotFile,
                #annot_file = annotFile, 
                colortable = ctab,
                #overlay_range = (-4,4),
                #screenshot_stem = freesurfer_dir+'_'+side,
                six_images = True,
                surface="pial")
        res = shots.run()
        createdList = res.outputs.snapshots
        for img in createdList:
            shutil.move(img, '/ccnc/'+os.path.basename(img))




def mergeLabel(freesurfer_dir, labelNames):
    for side in ['lh','rh']:
        #command = 'mri_mergelabels {inLabel} -o {outLabel} 2>/dev/null'.format(
        command = 'mri_mergelabels {inLabel} -o {outLabel}'.format(
            inLabel = ' '.join(['-i '+os.path.join(freesurfer_dir,'tmp',side+'.'+x+'_new.label') for x in labelNames]),
            outLabel = os.path.join(freesurfer_dir,'tmp',side+'_all.label'))

        print command
        os.popen(command).read()


def cleanMean(meanCSV, indCSV):
    meanDf = pd.read_csv(meanCSV, index_col=0)
    meanDf['roi'] = meanDf.subroi.str[3:]
    meanDf['side'] = meanDf.subroi.str[:2]
    meanDf = meanDf.groupby(['roi','side']).mean().reset_index()



    df = pd.read_csv(indCSV, index_col=0)
    df['roi'] = df.subroi.str[3:]
    df['side'] = df.subroi.str[:2]




    mergedDf = pd.merge(meanDf,
                        df,
                        on=['roi','side'],
                        how='inner')

    mergedDf['mean_sub_indv'] = (mergedDf.thickness_x - mergedDf.thickness_y)
    mergedDf['mean_sub_indv_cov'] =  rescale(mergedDf['mean_sub_indv'])
    return mergedDf

def rescale(npArray, new_min = 70, new_max = 100):
    values = npArray.tolist()
    print values
    output = []
    old_min, old_max = min(values), max(values)

    for v in values:
        new_v = (new_max - new_min) / (old_max - old_min) * (v - old_min) + new_min
        output.append(new_v)

    return output

def main(freesurferLoc,indcsv):
    mergedDf = cleanMean('/ccnc_bin/meanThickness/detailed_mean_2015_12_28.csv', indcsv)
    #standard_ctab = '/ccnc_bin/meanThickness/standard_ctab.txt'
    standard_ctab = '/Applications/freesurfer/FreeSurferColorLUT.txt'
    subject_ctab = os.path.join(freesurferLoc,'ctab.txt')
    shutil.copyfile(standard_ctab, subject_ctab)


    for row in mergedDf.iterrows():
        roiName = row[1].side + '.' + row[1].roi + '.label'
        roiLoc = os.path.join(freesurferLoc,
                'tmp',roiName)

        newNum = row[1].mean_sub_indv
        ctedit(subject_ctab,
                row[1].side+'-'+row[1].roi,
                row[1].mean_sub_indv_cov)
        #labelValueSwap(roiLoc, row[1].mean_sub_indv)

    #mergeLabel(freesurferLoc, mergedDf.roi.unique())
    makeBrainPic(freesurferLoc, 'all', subject_ctab)


if __name__=='__main__':
    main('/Volumes/CCNC_3T_2/kcho/ccnc/GHR_project/NOR04_JJW/FREESURFER',
            '/Volumes/CCNC_3T_2/kcho/ccnc/GHR_project/NOR04_JJW/FREESURFER/tmp/thick_kev_detailed.csv')

