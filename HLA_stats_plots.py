from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt

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



    def plot(self, to_file='HLA_data_growth.png'):

        #shorthand
        d = self.__stats

        #plot
        d.plot(x=['Month'], 
            y=['Alleles', 'Component Entries'], 
            figsize=(16, 8), 
            alpha=0.6, 
            linewidth=6, 
            rot=90)
        
        ax = plt.gca()
        
        
        #display one tick label every 3 releases (9 calendar months in between)
        #displaying every release makes it too 'busy'
        
        _release_intervals = 3
        
        xticks = [x for x in range(0, len(d['Month']), _release_intervals)]
        xticklabels = d['Month'][::_release_intervals]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_title(label='IMGT HLA Data Growth', fontsize=16)
        
        xticks = xticks[::3]
        xticklabels = xticklabels[::3]
        
        plt.xlabel('Release Month', {'fontsize': 12})
        
        plt.legend(prop={'size':14})
        
        for tick, label in zip(xticks, xticklabels):
            ya = d.loc[d[d['Month'] == label].index, 'Alleles']
            plt.text(x=tick, y=ya - 1500, s=ya[0], alpha=0.6, fontsize=12)

            yb = d.loc[d[d['Month'] == label].index, 'Component Entries']
            plt.text(x=tick, y=yb + 2000, s=yb[0], alpha=0.6, fontsize=12)
        
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

        #read in the data
        d = pd.read_csv(url_or_text_file,
                sep='\t',
                header=0,
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


    
    def plot(self, to_file='Allele_growth_by_locus.png'):
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
            plt.text(s=k, 
                    x=x - 5, 
                    y=y-v*2/3, 
                    alpha=.8, 
                    fontdict={'size':14, 'alpha':.5}) 
            
        plt.savefig(to_file)