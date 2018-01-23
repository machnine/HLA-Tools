import re
import os
from datetime import datetime
import pandas as pd

'''
This module processes and parses the IMGT HLA protein alignments
files. The latest release can be downloaded from:   
https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments

Lines 0 - 4 are the 'header' lines containing all of the meta info
Line 5 breaks up the 'header' and the main 'body'
Line 6 and other similar lines indicate the position of a particular
amino acid. Positions < 1 indicate the leader peptide sequences.
There is always an empty following the 'Prot    [position]' line
Line 8 is the actual beginning of a section of sequence alignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[0]HLA-DPA1 Protein Sequence Alignments\n
[1]IPD-IMGT/HLA Release: 3.30.0\n
[2]Sequences Aligned: 2017 October 27\n
[3]Steven GE Marsh, Anthony Nolan Research Institute.\n
[4]Please see http://hla.alleles.org/terms.html for terms of use.\n
[5]\n
[6] Prot              -40                                         1\n
[7]                   |                                           |\n
[8] DPA1*01:03:01:01           M RPEMIR AVLS FLLSLRGAGA... \n (max length=130)

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


To do:

The amino acid position numbering is not 'correct' if there is allele with
insertion(s) in the alignment file. C_Prot.txt for example. This does not seem
to affect usage, but the output position number can be confusing.

'''
META = {'ROWS' : 10,
        'LOCUS' : 0,
        'VERSION' : 1,
        'DATE' : 2,
        'SEQ_START' : 6,
        'FIRST_ALLELE' : 8}

WIDTHS = {'ALLELE':18, 'SEQ':122}

class Protein_Alignment:
    '''
    HLA protein alignment in Python Pandas Dataframe format.
    
    Parameters
    ----------
    alignment_file: alignment flat file (*.txt) file downloaded from IMGT
    alleles: keyword argument, List like
             Use this to limit the number of alleles to process to reduce time
             in computing. Processing 'DPB1_Prot.txt' (v3.30) takes about 800ms, 
             when limiting the number of alleles, the time cab be reduced to as
             low as 300ms. Processing 'ClassI_Prot.txt' (v3.30) takes 3200ms, 
             when limiting the number of alleles, the processing time can be 
             reduced to as low as 800ms. 
             Tested on PC with i5 4460 CPU and 8Gb RAM. 
             
    Examples
    --------
    >>> all_dpb = Protein_alignment('DPB1 Prot.txt')
    >>> 
    ...
    >>> DP_List = ['DPB1*01:01:01:01', 'DPB1*03:01:01:01', ... 'DPB1*135:01']
    >>> some_dpb = Protein_alignment('DPB1 Prot.txt', alleles=DP_List)
    >>>
    ...  
    '''
  
    
    def __init__(self, alignment_file, **kwargs):
        #this may require more input validations
        self.__allele_list = kwargs.get('alleles', None)
        if os.path.exists(alignment_file):
            self.__alignment_file = alignment_file
        else:
            raise Exception(FileNotFoundError("The file '{}' "
                            "cannot be found!".format(alignment_file)))

        self.__meta_data = self.__get_meta_data()
    
              
    ###Properties###    
    @property
    def meta(self):
        '''return meta data'''
        return self.__meta_data


    @property
    def aligned(self):
        '''
        Return the alignment in a pandas dataframe
        Speed could potentially be an issue here as processing the alignment
        file for all of the HLA class I proteins (Class I Prot.txt) can 
        take up to 160s on a PC with an i5 processor
        '''
        return self.__get_alignment(self.__alignment_file)

    
    ###Functions###  

    def __get_meta_data(self):
        '''
        Return the meta data of the protein alignment file
        Line numbers for the meta data are hard coded in
        Although this is not a good practice, because it's format is not 
        specified in any of the IMGT official documents, they have not 
        changed for at least 10 years based on all previous versions of 
        alignment files.
        '''

        patterns = {'LOCUS' : '(.*) Protein Sequence Alignments',
                    'VERSION' : 'IPD-IMGT/HLA Release: (\d+\.\d+.\d+)\n',
                    'DATE' : 'Sequences Aligned: (.*)\n'}
        
        temp = {}
        
        data = self.__read_raw_text_rows(self.__alignment_file)
        
        #protein start
        temp['PROT_START'] = self.__protein_start(data)
        #change data index
        
        for k, v in patterns.items():
            try:
                temp[k] = (re.match(v, data[META[k]])[1])
                if k == 'DATE':
                    temp[k] = self.__convert_date(temp[k])
                elif k == 'VERSION':
                    temp[k] = self.__convert_version(temp[k])
                
            except TypeError:
                print('This is not a protein alignment '
                                'file or it is corrupted!')

        return temp
    

    def __convert_version(self, versionstring):
        '''
            convert version string to tuple
        '''
        version = [int(x) for x in versionstring.split('.')]

        return tuple(version)
        
    

    def __convert_date(self, datestring):
        '''
            convert datestring to datetime.datetime()
        '''
        try:
            return datetime.strptime(datestring, '%Y %B %d')
        except Exception as err:
            print('Date time convertion error:', err)

            

    def __protein_start(self, meta_chunk: str):
        '''
        grab the start of mature protein (1 indexed)
        '''
        #' Prot             -30                     1\n' 
        loc = meta_chunk[6][WIDTHS['ALLELE']:]
        #'                   |                      |\n',
        #' A*01:01:01:01     MAVM APRTLLLS GALALTWA GSH
        seq = meta_chunk[8][WIDTHS['ALLELE']:]

        start_loc = list(loc).index('\n')
        start_seq = seq[:start_loc - 1]
        
        return len([x for x in start_seq if x !=' '])
    
    
    
    def __read_raw_text_rows(self, alignment_file):
        with open(alignment_file, 'r') as fin:
             data = [fin.readline() for x in range(META['ROWS'])]
        return data

    
    def __is_row_to_keep(self, input_string: str):
        '''
        Try to match each row to find an HLA allele like string:
            e.g.           A*01:
                        DRB1*01:
                        TAP1*01:
        if found return True to keep the row
        '''
        if re.search('\w+\*\d{2,3}:', input_string):
            return True
        else:
            return False
    

    def __relabel_columns(self, df: pd.DataFrame):
        '''
        relabel columns so position 1 of seq is the beginning amino acid of
        the mature protein, minus positions being the lead peptide
        '''
        start = self.__meta_data['PROT_START']
        length = len(df.columns)
        columns = list(range(-start, length - start + 1))
        columns.remove(0)
        return columns



    def __get_alignment(self, alignment_file):
        d = pd.read_fwf(alignment_file, 
                widths=[WIDTHS['ALLELE'], WIDTHS['SEQ']],  #total = 130 max
                header=None,                    
                names=['allele', 'seq'],
                converters={'seq': lambda x: x.replace(' ', '') 
                }
            )
        
        #fill all na with specific strings so later data manipulation 
        #can use string operations
        d.fillna('', inplace=True)
        
        #filter rows to keep
        kept_alleles = d.allele.apply(self.__is_row_to_keep)
        d = d[kept_alleles] 
        
        
        #limit the number of alleles to the given list to save time        
        if self.__allele_list:
            #the reference seq
            alleles_to_process = [d['allele'].iloc[0]]
            
            #extend to the user specified list of seq
            alleles_to_process.extend(self.__allele_list) 
            
            
            #filter based on above list
            d = d[d['allele'].isin(alleles_to_process)]
        
        
        #amalgamate (group) all rows with the same allele name
        # !!! sort = False !!! important, else the reference seq
        # will be incorrect because the default sort order is
        # alpha numeric, which means in HLA-DQB1, DQ2 is the reference
        # instead of DQB1*05:01:01:01 in the alignment file
        d = d.groupby('allele', sort=False).sum()
        
            
        
        #convert seq string to list
        d['seq'] = d['seq'].apply(lambda x: list(x))
        
        d = pd.DataFrame(data=d['seq'].values.tolist(), # d['seq'] to list 
                         index=d.index                  #use the same index
                        )
        
        #go through each column, if a cell is '-' replace it with the ref seq
        #this is at least 10 time faster than DataFrame.replace()
        for col in d.columns:
            d[col][d[col] == '-'] = d[col].iloc[0]
        
        d.columns = self.__relabel_columns(d)

        return d
