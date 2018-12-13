# freesurferSummary

---

## How to use in D1

```sh
/ccnc_bin/anaconda2/bin/python /ccnc_bin/freesurfer_summary.py -f ${FREESURFER_DIR}

#eg) 
/ccnc_bin/anaconda2/bin/python /ccnc_bin/freesurfer_summary.py -f /home/kangik/KANGIK/FREESURFER

# Then prompt will ask for the enter the initial of the subject
```

- This will create a graph of thickness summary in `${FREESURFER_DIR}`

### To do
- freesurferDir search : os.walk add



## Figure

```py
f = FreesurferFigureVolume(cortex_df, ['NOR', 'FEP']).linegraph()
f = FreesurferFigureVolume(cortex_df, ['NOR', 'FEP']).linegraph_side_mean()
f = FreesurferFigureVolume(cortex_df, ['NOR', 'FEP']).catplot()
f = FreesurferFigureVolume(cortex_df, ['NOR', 'FEP']).catplot_side_mean()

f = FreesurferFigureThickness(cortex_df, ['NOR', 'FEP']).linegraph()
f = FreesurferFigureThickness(cortex_df, ['NOR', 'FEP']).linegraph_side_mean()
f = FreesurferFigureThickness(cortex_df, ['NOR', 'FEP']).catplot()
f = FreesurferFigureThickness(cortex_df, ['NOR', 'FEP']).catplot_side_mean()
```
