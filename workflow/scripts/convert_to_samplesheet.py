#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
import numpy as np

# pylint: disable=undefined-variable

# logging
sys.stderr = open(snakemake.log[0], "w")

def convert_to_samplesheet(inp, out):
    
    # read extended sample sheet
    ex_sheet = pd.read_excel(
        io = inp
    )
    # validate if mandatory variables are not missing
    if not ex_sheet.loc[:, ['project_id', 'sample_id', 'index_seq']].isnull().any().any():
       
        # hardcode header of sample sheet
        header_dict={
            '0': ['[Header]', 'IEMFileVersion', 'Investigator Name', 'Experiment Name',
                  'Date', 'Workflow', 'Application', 'Assay', 'Description', 'Chemistry',
                  '', '[Reads]', '', '', '', '[Settings]', '', '[Data]', 'Sample_ID'],
            '1': np.repeat(['', 'Sample_Name'], [18,1]).tolist(),
            '2': np.repeat(['', 'Sample_Plate'],[18,1]).tolist(),
            '3': np.repeat(['', 'Sample_Well'], [18,1]).tolist(),
            '4': np.repeat(['', 'I7_Index_ID'], [18,1]).tolist(),
            '5': np.repeat(['', 'index'], [18,1]).tolist(),
            '6': np.repeat(['', 'Sample_Project'], [18,1]).tolist(),
            '7': np.repeat(['', 'Description'], [18,1]).tolist(),
            '8': np.repeat('', [19]).tolist(),
            '9': np.repeat('', [19]).tolist()
        }
        header_df = pd.DataFrame(header_dict)
        
        # add contents of extended sheet
        ex_sheet_dict = {
            '0': ex_sheet['sample_id'].tolist(),
            '5': ex_sheet['index_seq'].tolist(),
            '6': ex_sheet['project_id'].tolist()
        }
        sheet_df = pd.DataFrame(ex_sheet_dict)
        
        # combine and write sample sheet
        sheet = pd.concat([header_df, sheet_df], ignore_index=True).fillna('')
        sheet.to_csv(
            out,
            header = None,
            index = None
        )
    else:
        sys.exit('error: check if there are no empty entries in the extended_sample_sheet.xlsx')


convert_to_samplesheet(
    inp = snakemake.input[0],
    out = snakemake.output[0]
)

