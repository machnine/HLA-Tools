# proline
Read an HLA protein sequence flat text files from IMGT HLA's Github (https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments) and convert them to Pandas.DataFrame objects.

## download the lastest release:
```
$python download_latest.py
```

## import the Protein_Alignment object
```python
from proline import Protein_Alignment
```

## create a protein alignment object:
```python
hla_dpb1 = Protein_Alignment('DPB1_prot.txt')

#dataframe object
hla_dpb1.aligned
```

## only align alleles in a list:
```python
list_of_hla_a = ['A*01:01:01:01', 'A*02:01:01:01','A*33:01:01:01']
some_hla_a = Protein_Alignment('A_prot.txt', alleles = list_of_hla_a)

#dataframe object
some_hla_a.aligned
```

## displaying meta data:
```python
hla_dpb1 = Protein_Alignment('DPB1_prot.txt')
#dictionary
print(hla_dpb1.meta)
```

## return unique protein sequences:
```python
hla_dpb1 = Protein_Alignment('DPB1_prot.txt')
#dataframe object
hla_dpb1.unique_seq()
```

### or unique sequences within a given range:
```python
#dataframe object
hla_dpb1.unique_seq(aa_range=[4, 84])
```

