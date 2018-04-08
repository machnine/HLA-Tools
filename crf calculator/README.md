Calculated reaction frequency (cRF) is an indicator of the level of sensitisation for a patient used by [NHSBT-ODT](https://www.odt.nhs.uk) in the UK. It is calculated by finding the percentage of blood group compatible, HLA-incompatible donors in a pool of 10,000 donors. In all of the versions released so far, only HLA-A, B, C, DRB1/3/4/5 and DQB1 are taken into account; HLA-DQA1, DPA1 and DPB1 have not yet been included. The latest version of the calculator in the form of a Microsoft Excel binary file can be downloaded from the [Tools page](https://www.odt.nhs.uk/transplantation/tools-policies-and-guidance/calculators/) at ODT: 

[Calculated Reaction Frequency tool (XLSB)](https://nhsbtdbe.blob.core.windows.net/umbraco-assets-corp/5948/hla-mm-and-crf.xlsb)

Although the tool itself is useful, it is very slow and is not compatible with iPads and other tablets. This module extracts all of the 10,000 donor HLA types from the xlsb file and store them in a database for faster calculation and integration into apps without the requirement of Excel support.

## Dependancies
    pandas
    pyxlsb

## Download the donor information into a database (run once):
```python
    python xlsb_to_db_converter.py
```

## Calculate crf:
```python
import crf

#recipient blood group
blood_group = 'A'

#recipient unacceptable antigen specificities
unacceptable_antigens = ['A2', 'A3', 'DR7', 'DQ1']

crf_val = crf_cal(bg=blood_group, ua=unacceptable_antigens)
# returns
# 0.9202195018995357
```