#__author__ = 'kcho'
##import matplotlib.pyplot as plt
##plt.style.use('ggplot')
### from freesurfer_summary import openStatsTable

##import read_freesurfer_stats

###df = read_freesurfer_stats.openStatsTable('/Users/kcho/T1/FREESURFER')
###print df

##ICV = read_freesurfer_stats.getICV('/Users/kcho/T1/FREESURFER')
##print ICV

##df = read_freesurfer_stats.openStatsTable_big('/Users/kcho/T1/FREESURFER')
##print df
##fig = read_freesurfer_stats.graph_ind(df, 'surfarea')
##plt.show()





#from freesurfer_summary import *
#import os

## freesurfer_dir = '/Users/kcho/T1/FREESURFER'
## os.environ["FREESURFER_HOME"] = '/Applications/freesurfer'
## os.environ["SUBJECTS_DIR"] = '{0}'.format(os.path.dirname(freesurfer_dir))
##
## print freesurfer_summary.collectStats_v2('/Users/kcho/T1/FREESURFER')
## roiDict = freesurfer_summary.get_cortical_rois()
## print roiDict
## roiDict2 = freesurfer_summary.get_cortical_rois_detailed()
## print roiDict2
##
## print 'range(len(roiDict.keys())*2)', range(len(roiDict.keys())*2)
##
## freesurfer_summary.makeLabel(freesurfer_dir)
## freesurfer_summary.mergeLabel(freesurfer_dir, roiDict)
##
##
## freesurfer_summary.mergeLabel(freesurfer_dir,roiDict)
## infoDict = freesurfer_summary.getInfoFromLabel(freesurfer_dir, roiDict)
## infoDf = freesurfer_summary.dictWithTuple2df(infoDict)
##
## print infoDf


##main_freesurferDir = '/Users/kcho/T1/FREESURFER'
#main_freesurferDir = '/Volumes/promise/nas_BackUp/CCNC_MRI_3T/DNO/DNO46_KJU'

#main_freesurferDir = get_freesurferDir(main_freesurferDir)

#print main_freesurferDir
#os.environ["FREESURFER_HOME"] = '/Applications/freesurfer'
#os.environ["SUBJECTS_DIR"] = '{0}'.format(os.path.dirname(main_freesurferDir))


#meanDf = pd.read_csv('/ccnc_bin/meanThickness/detailed_mean_2015_12_28.csv', index_col=0)
#roiDict = get_cortical_rois()

#infoDf = collectStats_v2(main_freesurferDir)#background_subject_locs)
#print infoDf
#subjName = 'ha'


## print infoDf

#freesuferList = subjDirs_to_fsDirs(['/Volumes/promise/nas_BackUp/CCNC_MRI_3T/DNO/DNO46_KJU',
                                    #'/Volumes/promise/nas_BackUp/CCNC_MRI_3T/DNO/DNO46_KJU'])
#print freesuferList
#concatDf = concatFsDf(freesuferList)
#meanDf = concatDf_to_meanDf(concatDf)

#draw_thickness_detailed(infoDf,
                        #meanDf,
                        #'prac_subj',
                        #'CCNC_mean')


##pd.concat((df1, df2), axis=1).mean(axis=1)
## draw_thickness_detailed(infoDf,
##                         meanDf,
##                         subjName,
##                         'HCs')
## valueSwap.main(main_freesurferDir,
##                os.path.join(main_freesurferDir,
##                             'tmp/thick_kev_detailed_new.csv'))

