import re
import os
from collections import Iterable
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

- The amino acid position numbering is not 'correct' if there is allele with
  insertion(s) in the alignment file. C_Prot.txt for example. This does not seem
  to affect usage, but the output position number can be confusing.

- Ability to choose any reference seq for alignment instead of just the default.

'''

#version number
VERSION = 'proline: HLA protein aligment utilities version: ' + '0.0.2'

#meta data to parse
META = {'ROWS' : 10,
        'LOCUS' : 0,
        'VERSION' : 1,
        'DATE' : 2,
        'SEQ_START' : 6,
        'FIRST_ALLELE' : 8}

#read the alignment file in as fixed width text
WIDTHS = {'ALLELE':18, 'SEQ':122}

#suffixes for Null/Question?/Secreted proteins
EXPRESSION_EXCLUSION = 'NQS'


class Protein_Alignment:
    '''
    HLA protein alignment in Python Pandas Dataframe format.
    
    Parameters
    ----------
    alignment_file      : alignment flat file (*.txt) file downloaded from IMGT
    alleles             : limit to the list of alleles to align
    ignore_non_expressed: flag to exclude alleles with sufixes defined in EXPRESSION_EXCLUSION
    
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
        self.__ignore_non_expressed = kwargs.get('ignore_non_expressed', False)
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
        Return the alignment as a pandas dataframe
        Speed could potentially be an issue here as processing the alignment
        file for all of the HLA class I proteins (Class I Prot.txt) can 
        take up to 160s on a PC with an i5 processor
        '''
        return self.__get_alignment(self.__alignment_file)
    
    
    ###Public Functions###
    def unique_seq(self, aa_range=None):
        '''
        uniquefy sequeneces in the given aa_range 
        return them in a dataframe without allele indices
        '''
        if not aa_range:
            #no aa range specified
            pass
        elif isinstance(aa_range, Iterable):
            #only start and finish given
            if len(aa_range) == 2:
                aa_range = list(range(aa_range[0], aa_range[1] + 1))
            #wrong param
            elif len(aa_range) < 2:
                raise ValueError('aa_range size < 2, please specify the start and finish of AA positions')
            #individual positions given in a list
            else:
                pass
        else:
            raise ValueError('aa_range must be an Iterable object')
            
        aligned = self.__get_alignment(self.__alignment_file)
        
        return self.__unique_protein_seq(aligned, aa_range)

    
    ###Private Functions###  

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
        '''
        read the raw text alignment file row by row into a list
        '''
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
        if ignore_non_expressed flag = True also exclude all of the
        unknown 'Q', secreted 'S' and null 'N' alleles defined in EXPRESSION_EXCLUSION.
        '''
        
        looks_like_hla = re.compile('(\w+\*\d{{2,3}}.*:\d{{2,3}})([{}]?)'.format(EXPRESSION_EXCLUSION))
        
        matches = looks_like_hla.search(input_string)
        if matches:
            if self.__ignore_non_expressed:
                hla_name, suffix = matches.groups()
                if len(suffix) >0:
                    return False
                return True
            else:
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

    
    def __unique_protein_seq(self, df, aa_range):
        '''
        uniquefy a given protein sequence alignment dataframe
        return a new dataframe without allele labelling
        but with the same column (aa position) labels
        '''
        if aa_range:
            columns = aa_range
            df = df[aa_range]
        else:
            columns = df.columns
        
        #remove null alleles and 'Q', and 'S' alleles
        df = df.loc[[x for x in df.index if not (x[-1] in ('S', 'Q', 'N'))]]

        #merge individual seq into a string
        df['seq'] = df.apply(lambda x: ''.join(x), axis=1)

        #uniquefy the datafame data 
        df = pd.DataFrame([pd.Series(list(x)) for x in set(df['seq'])])

        #relabelling columns
        df.columns = columns

        return df