from pyxlsb import open_workbook
import pandas as pd
import io
import sqlite3
import requests


URL = 'https://nhsbtdbe.blob.core.windows.net/umbraco-assets-corp/5948/hla-mm-and-crf.xlsb'

def read_crf_xlsb(xlsb_file):
    '''
    read the crf calculator file from odt and convert it to a pandas dataframe
    '''
    data = []
    with open_workbook(xlsb_file) as wb:
        #get the sheet with donor types (hopefully only 1)
        donor = [s for s in wb.sheets if 'Donor' in s][0]
        #open the donor type sheet
        with wb.get_sheet(donor) as sheet:
            #find the header
            header = [x.v for x in next(sheet.rows())]
            #find the begining and the end of HLA type with Blood Group in front
            #assuming the first 'BG' and 'DQ4' are start and finish
            column_range = [header.index(x) for x in header if 'BG' in x or 'DQ4' in x][:2]
            if len(column_range) != 2 or (column_range[1] - column_range[0] <= 0):
                raise Exception('Input xlsb seems to be in an unrecognisable format!')
            
            for row in sheet.rows():
                data.append([x.v for x in row[column_range[0]:column_range[1] + 1]])
    df = pd.DataFrame(data, dtype=int)
    df.rename(columns=df.iloc[0].apply(lambda x:x.upper()), inplace=True)
    df.drop(df.index[0], inplace=True)
    df.fillna(0, inplace=True)
    return df

def crf_xlsb_to_db(url=URL, db='ten_k_donors.db'):
    '''
    output the converted dataframe to the donor database
    '''
    #create download connection
    dl = requests.get(url)
    #download the file, convert to byte stream then dataframe
    df = read_crf_xlsb(io.BytesIO(dl.content))
    #output to database
    con = sqlite3.connect(db)
    df.to_sql(con=con, name='donors')


if __name__ == '__main__':
    crf_xlsb_to_db()
