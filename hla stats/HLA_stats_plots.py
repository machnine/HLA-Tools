from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import requests

class DataGrowthPlot():
    '''
    Plot the data growth stats from EBI/IMGT/HLA
    '''
    def __init__(self, url):
        stats = pd.read_html(url, parse_dates=True)
        if stats[1].iloc[0, 0] != 'Database Growth':
            raise Exception('The page layout of {url} has changed - unable to parse data.')
        else:
            d = stats[1].iloc[2:].copy()

            #reset columns
            d.columns = list(stats[1].iloc[1])

            #reset index
            d.set_index('Release', drop=True, inplace=True)

            #convert datatypes
            d.loc[:, ['Alleles','Component Entries']] \
            = d.loc[:, ['Alleles','Component Entries']].astype(int)

            self.__stats = d



    def plot(self, to_file='HLA_data_growth.png', savefig=False):

        #shorthand
        d = self.__stats

        #plot
        d.plot(figsize=(16, 8), 
               alpha=0.6, 
               linewidth=6, 
               rot=90)
        
        ax = plt.gca()
        
        #display one tick label every 3 releases (9 calendar months in between)
        #displaying every release makes it too 'busy'
        
        _intervals = 3
        _len = len(d)
        _month = d['Month']
        
        xticks = [x for x in range(_len)]
        xticklabels = _month
        
        #set ticks = interval plus the latest release
        ax.set_xticks(xticks[::_intervals] + xticks[-1:])
        ax.set_xticklabels(xticklabels[::_intervals].tolist() + xticklabels[-1:].tolist())
        
        ax.set_title(label='IMGT HLA Data Growth', fontsize=16)       
        plt.xlabel('Release Month', {'fontsize': 12})
        plt.legend(prop={'size':12})
        
        x_offset = -1.5
        y_offset = 500
        for tick, label in zip(xticks, xticklabels):
            if (tick % (_intervals ** 2) == 0) or (tick == xticks[-1]):
                ya = d.loc[d[_month == label].index, 'Alleles']
                plt.text(x=tick+x_offset, y=ya+y_offset, s=ya[0], alpha=0.6, fontsize=12)

                yb = d.loc[d[_month == label].index, 'Component Entries']
                plt.text(x=tick+x_offset, y=yb+y_offset, s=yb[0], alpha=0.6, fontsize=12)
        if savefig:
            plt.savefig(to_file)


    def release_dates(self):
        return self.__stats['Month']



class LocusStackingPlot():
    '''
    Plot the 'Allelelist_history.txt' file of the latest release
    '''
    def __init__(self, url_or_text_file, version_to_date_map=None):
        if 'http' in url_or_text_file:
            print('Reading directly from url may take a while...')
        
        header_line, separator = self.__get_history_file_header_type(url_or_text_file)
        
        #read in the data
        d = pd.read_csv(url_or_text_file,
                        sep=separator,                        
                        skiprows=header_line,
                        index_col='HLA_ID',
                        dtype=str)
        

        #fill the gaps with empty strings
        d.fillna('', inplace=True)
        #keep only the locus names
        d = d.applymap(lambda x: x.split('*')[0].replace('w', ''))

        #if a conversion mapping is available
        if  version_to_date_map is None:
            #don't do any conversion
            f = lambda x:x 
        else:
            #lambda to convert version numbers
            # e.g. 3310 -> 3.31, 1050 -> 1.5 etc.
            f = lambda x: version_to_date_map['.'.join( \
                                               [x[0], str(int(x[1:3]))]) \
                                             ]
            

        #count each locus
        c = d.apply(lambda x: Counter(x))

        #convert c to dataframe and label according to version_to_date_map
        new_d = pd.DataFrame()

        for value in c.index[::-1]:
            new_index = f(value)
            new_d = new_d.append(pd.DataFrame(c.loc[value], index=[new_index]))

        self.__stats = new_d


    def __get_history_file_header_type(self, url_or_file):
        '''
        function to find where the file really starts
        after release 3.32 the allele hisotry file changes
        from tsv to csv with meta information on top
        '''
        header = self.__read_header(url_or_file)
        for n, line in enumerate(header):
            if 'HLA_ID' in line:
                return n, line[len('HLA_ID'):len('HLA_ID')+1]
    
    
    def __read_header(self, url_or_file):
        if 'http://' in url_or_file or 'https://' in url_or_file:
            return self.__read_header_from_url(url_or_file)
        else:
            return self.__read_header_from_file(url_or_file)
                
    
    
    
    def __read_header_from_file(self, file_name):
        '''
        helper function read the first 100 lines from 
        allele history text file
        ''' 
        with open(file_name, 'r') as f:
            hundred_lines = f.readlines()[:100]
        return hundred_lines        
    
    
    def __read_header_from_url(self, url):
        '''
        helper function read the first 100 lines from
        allele history text file URL
        '''
        with requests.get(url) as r:
            hundred_lines = r.text.splitlines()[:100]
        return hundred_lines

    
    def plot(self, to_file='Allele_growth_by_locus.png', savefig=False):
        #only interested in the 'classic' HLA
        d = self.__stats[['A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1']]
        #ignore minor version changes take the max values
        d = d.reset_index().groupby('index').agg('max')
        plt.figure(figsize=(16, 8))
        plt.stackplot(d.index, d.T, alpha=.8)
        plt.xlim([0, len(d.index)-1])
        plt.xlabel('Release Version', {'fontsize': 12})
        ax = plt.gca()
        ax.set_title(label='HLA Allele Number Growth By Locus', fontsize=16)
        ax.set_xticks(range(len(d.index)))
        ax.set_xticklabels(d.index, rotation=90)

        y = 0
        x = len(d)
        for k, v in d.iloc[-1].items():
            y += v
            plt.text(s=f'{k}: {v}', 
                     x=x-1, 
                     y=y-v*2/3, 
                     alpha=.8, 
                     ha='right',
                     fontdict={'size':14, 'alpha':.5},
                     bbox=dict(boxstyle="round, pad=0.1",
                       ec='none',
                       fc=(1, 1, 1, 0.5),
                       )) 
        if savefig:
            plt.savefig(to_file)