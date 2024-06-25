import os
import zipfile

import pandas as pd
import streamlit as st

st.set_page_config(page_title='DIANN Output PCQ Converter')

st.title('DIANN Output PCQ Converter')
st.caption('This app converts DIANN output to PCQ format')

csv_files = st.file_uploader('Upload a CSV file', type=['csv'], accept_multiple_files=True)

if not csv_files:
    st.warning('Please upload at least one CSV file')
    st.stop()

dfs = []
for csv_file in csv_files:
    df = pd.read_csv(csv_file, index_col=0)
    dfs.append(df)

df = pd.concat(dfs)
file_names = df['File.Name'].unique()

group1_files = st.multiselect('Select group 1 files', file_names)
remaining_files = [file for file in file_names if file not in group1_files]
group2_files = st.multiselect('Select group 2 files', remaining_files)

# add sequential number to the dataframe
df['PSM.ID'] = range(1, len(df) + 1)
df['Channel.Score'] = 1 - df['Channel.Q.Value']

# set psm id to be indexed
cols_to_keep = ['PSM.ID', 'Stripped.Sequence', 'Light/Heavy.Ratio', 'Quantity.Quality', 'Protein.Group']
group1_df = df[df['File.Name'].isin(group1_files)][cols_to_keep].set_index('PSM.ID')
group2_df = df[df['File.Name'].isin(group2_files)][cols_to_keep].set_index('PSM.ID')

t1, t2 = st.tabs(['Group 1', 'Group 2'])
t1.subheader('Group 1')
t1.dataframe(group1_df, use_container_width=True)
t2.subheader('Group 2')
t2.dataframe(group2_df, use_container_width=True)

try:
    # Save the two dfs to files and then zip them together and have this as a download file (st.download_button)
    group1_df.to_csv('group1.tsv', sep='\t')
    group2_df.to_csv('group2.tsv', sep='\t')

    with zipfile.ZipFile('groups.zip', 'w') as zipf:
        zipf.write('group1.tsv')
        zipf.write('group2.tsv')

    with open('groups.zip', 'rb') as f:
        st.download_button('Download ZIP', f, file_name='groups.zip', type='primary', use_container_width=True)
except Exception as e:
    st.error(f'An error occurred: {e}')

finally:
    # Cleanup the CSV files after zipping
    os.remove('group1.tsv')
    os.remove('group2.tsv')

    # Cleanup the ZIP file after download
    os.remove('groups.zip')
