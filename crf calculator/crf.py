import pandas as pd
import sqlite3

def compatible_blood_groups(bg: str):
    '''
    return a list of compatible blood groups
    '''
    bg = bg.upper()
    blood_groups = ['A', 'AB', 'B', 'O']
    if bg not in blood_groups:
        raise ValueError(f'{bg} is not a valid blood group!')

    compatible = {
        'A': ['A', 'O'],
        'B': ['B', 'O'],
        'AB': ['A', 'AB', 'B', 'O'],
        'O': ['O']
    }

    return compatible[bg]


def get_donors(db='ten_k_donors.db'):
    '''
    read the donor type info from the database
    '''
    with sqlite3.connect(db) as con:
        donors = pd.read_sql_query('select * from donors', con=con)
    return donors


def crf_cal(bg: str, ua: list):
    '''
    calculate crf based on the given blood group:bg
    and list of unacceptable antigens:ua
    '''
    #read the database
    donors = get_donors()
    #get donor with the compatible blood groups
    comp_bg = compatible_blood_groups(bg)
    comp_donors = donors[donors.BG.isin(comp_bg)]
    #donor with unacceptable antigens
    query = ' or '.join([f'{x}==1' for x in ua])
    comp = len(comp_donors.query(query))
    #total compatible donors
    total = len(comp_donors)
    return comp / total




