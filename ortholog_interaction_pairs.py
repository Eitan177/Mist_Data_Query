import numpy as np
import pandas as pd
import streamlit as st
from os.path import exists


import zipfile


texttomatch1=st.text_input('first gene to find',value='')

texttomatch2=st.text_input('second gene to find',value='')

filtonnum=st.radio('should we filter on minimum ortholog interaction number',options=['yes','no'],index=1)
numorth=st.selectbox('select a gene that is paired with another gene having a minimum of this many interactions in orthologs',options=[2,3,4,5,6,7,8,9,10],index=0)
showmim=st.checkbox('show mim annotations',value=False)

if texttomatch1 != '' or texttomatch2 != '' or filtonnum =='yes':
    file_exists = exists('mist_interlog_match_genes_sorted.csv.zip')
    if not file_exists:
        with zipfile.ZipFile('mist_interlog_match_genes_sorted.csv.zip', 'r') as zip_ref:   
            zip_ref.extractall('.')    
    ff=pd.read_csv('mist_interlog_match_genes_sorted.csv')
    if filtonnum == 'yes':
        ff=ff[ff['num']>=numorth]
    if texttomatch1 != '' and texttomatch2 != '':
        ff=ff[np.logical_or(ff['GeneA'].str.contains(texttomatch1,case=False,regex=False),ff['GeneB'].str.contains(texttomatch1,case=False,regex=False))]
        ff=ff[np.logical_or(ff['GeneA'].str.contains(texttomatch2,case=False,regex=False),ff['GeneB'].str.contains(texttomatch2,case=False,regex=False))]
        
    elif texttomatch1 != '':
        ff=ff[np.logical_or(ff['GeneA'].str.contains(texttomatch1,case=False,regex=False),ff['GeneB'].str.contains(texttomatch1,case=False,regex=False))]
    st.write('Protein-Protein Interaction file MIST result')            
    st.write('number of interactions found: '+str(ff[['GeneA','GeneB']].drop_duplicates().shape[0]))
    st.write(ff)
    

    if showmim:
            
        kk=pd.read_table('mimannot.txt',sep='\t')
        if texttomatch1 != '' and texttomatch2 != '':
            foundrow=np.where(((kk['OmimPresent'].str.contains(texttomatch1,case=False,regex=True)) | (kk['OmimAbsent'].str.contains(texttomatch1,case=False,regex=True))) & \
                    ((kk['OmimPresent'].str.contains(texttomatch2,case=False,regex=True)) | (kk['OmimAbsent'].str.contains(texttomatch2,case=False,regex=True)) ))[0]
        elif texttomatch1 != '':
            foundrow=np.where((kk['OmimPresent'].str.contains(texttomatch1,case=False,regex=True)) | (kk['OmimAbsent'].str.contains(texttomatch1,case=False,regex=True)) )[0]
        kk=kk.iloc[foundrow]
        st.write('OMIM')
        st.write('found',len(kk),'rows')

        st.write(kk)
