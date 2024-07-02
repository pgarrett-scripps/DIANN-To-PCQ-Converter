from collections import Counter
import os
import zipfile

import numpy as np
import pandas as pd
import streamlit as st
import scipy.stats as stats
import plotly.graph_objs as go
import statsmodels.stats.multitest as smm
import ntpath

st.set_page_config(page_title='DIANN Output PCQ Converter', initial_sidebar_state='expanded')

with st.sidebar:

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

    group_type = st.selectbox(label='Grouping', options=['Manual Selection', 'Auto File Split'])
    file_to_group = {}
    missing_selection = False
    if group_type == 'Manual Selection':
        num_groups = st.number_input(label='Number of Groups', value=2)

        groups = {}
        remaining_files = set(file_names)
        for i in range(1, num_groups+1):
            file_selection = st.multiselect(f'Select group {i} files', options=list(remaining_files), key=f'{i}_file_selection')
            remaining_files = set.difference(remaining_files, file_selection)
            file_to_group.update({file: i for file in file_selection})
            if len(file_selection) == 0:
                missing_selection = True
    else:
        num_groups = len(file_names)
        for i, file in enumerate(file_names, 1):
            file_to_group[file] = i

    if missing_selection is True:
        st.warning("Not all groups are filled out. Ensure that All groups have at least 1 file specified")
        st.stop()

# add sequential number to the dataframe
df['PSM.ID'] = range(1, len(df) + 1)
df['Channel.Score'] = 1 - df['Channel.Q.Value']
df['group'] = df['File.Name'].apply(lambda x: file_to_group.get(x))

# drop rows with missing group
df = df.dropna(subset=['group'])

# set psm id to be indexed
cols_to_keep = ['PSM.ID', 'Stripped.Sequence', 'Light/Heavy.Ratio', 'Quantity.Quality', 'Protein.Group', 'File.Name']

group_dfs = {}
for i in range(1, num_groups+1):
    group_dfs[i] = df[df['group'] == i][cols_to_keep].set_index('PSM.ID')

summary_group_dfs = {}
for i, group_df in group_dfs.items():
    summary_group_df = group_df.groupby(['Stripped.Sequence', 'Protein.Group'])['Light/Heavy.Ratio'].agg(
        mean_ratio='mean',
        median_ratio='median',
        std_ratio='std',
        sem_ratio=lambda x: x.std() / (len(x) ** 0.5),
        count='count'
    ).reset_index()

    summary_group_df['PSM.ID'] = list(range(len(summary_group_df)))
    summary_group_dfs[i] = summary_group_df.set_index('PSM.ID')

t1, t2, t3 = st.tabs(['Download', 'Groups', 'Volcano Plot'])
with t1:

    c1, c2 = st.columns([3,1])
    zip_basename = c1.text_input('Enter the name of the ZIP file', 'groups')
    zip_extension = c2.text_input(f'File Extension', f'.zip', disabled=True, key=f'zip_extension')
    zip_name = zip_basename + zip_extension
    if group_type == 'Manual Selection':
        group_names = {}
        summary_group_names = {}
        for i in range(1, num_groups+1):
            c1, c2 = st.columns([3,1])
            file_basename = c1.text_input(f'Enter the name of the group {i} file', f'group{i}')
            file_extenstion = c2.text_input(f'File Extension', f'.tsv', disabled=True, key=f'{i}_file_extension')
            group_names[i] = file_basename + file_extenstion
            summary_group_names[i] = group_names[i].replace('.tsv', '_summary.tsv')

    else:
        group_names = {}
        summary_group_names = {}
        for i in range(1, num_groups+1):
            fname = os.path.basename(file_names[i-1])
            fname = ''.join(fname.split('.')[0].split('/')[-1].split('\\')[-1])
            c1, c2 = st.columns([3,1])
            file_basename = c1.text_input(f'Enter the name of the group {i} file', f'{fname}')
            file_extenstion = c2.text_input(f'File Extension', f'.tsv', disabled=True, key=f'{i}_file_extension')
            group_names[i] = file_basename + file_extenstion
            summary_group_names[i] = group_names[i].replace('.tsv', '_summary.tsv')
    
    try:
        for i in range(1, num_groups + 1):
            group_dfs[i].to_csv(group_names[i], sep='\t')
            summary_group_dfs[i].to_csv(summary_group_names[i], sep='\t')

            unique_files = group_dfs[i]['File.Name'].unique()
            for file_path in unique_files:
                file_df = group_dfs[i][group_dfs[i]['File.Name'] == file_path]
                file_path = file_path.strip()
                file_path = ntpath.normpath(file_path)
                file_name = ntpath.basename(file_path)
                file_name_without_ext, file_ext = ntpath.splitext(file_name)
                file_df_name = f'group_{i}_{file_name_without_ext}.tsv'
                file_df.to_csv(file_df_name, sep='\t')
                group_names[f'{i}_{file_name_without_ext}'] = file_df_name

        with zipfile.ZipFile(zip_name, 'w') as zipf:
            for i in range(1, num_groups + 1):
                group_folder = f'group_{i}/'
                zipf.write(group_names[i], arcname=f'{group_folder}{os.path.basename(group_names[i])}')
                zipf.write(summary_group_names[i], arcname=f'{group_folder}{os.path.basename(summary_group_names[i])}')

                for file_path in group_dfs[i]['File.Name'].unique():
                    file_path = file_path.strip()
                    file_path = ntpath.normpath(file_path)
                    file_name = ntpath.basename(file_path)
                    file_name_without_ext, file_ext = ntpath.splitext(file_name)
                    file_df_name = group_names[f'{i}_{file_name_without_ext}']
                    zipf.write(file_df_name, arcname=f'{group_folder}{os.path.basename(file_df_name)}')

        with open(zip_name, 'rb') as f:
            st.download_button('Download ZIP', f, file_name=zip_name, type='primary', use_container_width=True)

    except Exception as e:
        st.error(f'An error occurred: {e}')
        raise e

    finally:
        for i in range(1, num_groups + 1):
            os.remove(group_names[i])
            os.remove(summary_group_names[i])
            for file_path in group_dfs[i]['File.Name'].unique():
                file_path = file_path.strip()
                file_path = ntpath.normpath(file_path)
                file_name = ntpath.basename(file_path)
                file_name_without_ext, file_ext = ntpath.splitext(file_name)
                key = f'{i}_{file_name_without_ext}'
                if key in group_names:
                    os.remove(group_names[key])

        os.remove(zip_name)


with t2:
    for i in range(1, num_groups+1):
        st.header(f'Group {i}', divider=True)
        for fname, j in file_to_group.items():
            if i == j:
                st.caption(fname)
        #st.dataframe(group_dfs[i], use_container_width=True)
        st.subheader(f'Summary Data', divider=True)
        st.caption('Click on a row to view all peptides.')
        selection = st.dataframe(summary_group_dfs[i], use_container_width=True, on_select='rerun', selection_mode='single-row')
        selected_row = selection['selection']['rows']
        if selected_row:
            st.subheader('All Data', divider=True)
            selected_row = selected_row[0]
            selected_sequence = summary_group_dfs[i].iloc[selected_row]['Stripped.Sequence']
            all_peptide_data = group_dfs[i][group_dfs[i]['Stripped.Sequence'] == selected_sequence]
            st.dataframe(all_peptide_data, use_container_width=True)
            

with t3:
    st.subheader('Volcano Plot')
    c1, c2 = st.columns(2)
    index_col = c1.selectbox('Select index column', ['Modified.Sequence', 'Stripped.Sequence'], index=1)
    values_col = c2.selectbox('Select values column', ['Light/Heavy.Log2Ratio', 'Light/Heavy.Ratio'], index=0)
    agg_func = st.selectbox('Select aggregation function', ['mean', 'median'], index=0)

    group_options = list(range(1, num_groups+1))
    group1_index = int(st.selectbox(label='Select First Group', options=group_options))
    group2_index = int(st.selectbox(label='Select Second Group', options=group_options))

    if group1_index == group2_index:
        st.warning('Ensure groups are not the same')
        st.stop()

    # Handle missing values
    df = df.dropna(subset=[index_col, values_col, 'group'])

    pivot_table = df.pivot_table(index=index_col, columns='group', values=values_col, aggfunc=agg_func)
    
    with st.expander('Show Full pivot Table'):
        st.dataframe(pivot_table)

    indexes_to_keep = list(set([group1_index, group2_index]))
    pivot_table = pivot_table[indexes_to_keep]

    # drop rows with missing values
    pivot_table = pivot_table.dropna()

    # Calculate p-values
    p_values = []
    for seq in pivot_table.index:
        group1 = df[(df[index_col] == seq) & (df['group'] == group1_index)][values_col].dropna()
        group2 = df[(df[index_col] == seq) & (df['group'] == group2_index)][values_col].dropna()
        t_stat, p_val = stats.ttest_ind(group1, group2)
        p_values.append(p_val)

    pivot_table['p_value'] = p_values

    valid_pivot_table = pivot_table.dropna()

    # Calculate q-values for valid rows
    valid_pivot_table['q_value'] = smm.multipletests(valid_pivot_table['p_value'], method='fdr_bh')[1]

    # Create a dataframe for rows with None p-values, setting their p_value and q_value to 1
    non_pvalue_rows = pivot_table[pivot_table['p_value'].isna()].copy()
    non_pvalue_rows['p_value'] = 1
    non_pvalue_rows['q_value'] = 1

    # Concatenate the valid and non-valid rows
    pivot_table = pd.concat([valid_pivot_table, non_pvalue_rows], ignore_index=False)

    # add back in Protein.Group (lookup Protein.Group based on the Sequence)
    sequence_to_group = df.set_index(index_col)['Protein.Group'].to_dict()
    pivot_table['Protein.Group'] = pivot_table.index.map(sequence_to_group)

    

    

    # Plotting the volcano plot
    fig = go.Figure()

    # Assuming group 1 log2 ratios as x-axis and -log10 of q-values as y-axis
    x = pivot_table.iloc[:, 0]
    y = -np.log10(pivot_table['q_value'])

    fig.add_trace(go.Scatter(x=x,
                             y=y,
                             mode='markers',
                             hovertext=pivot_table.index,
                             ))

    # Adding labels to the plot for better understanding
    fig.update_layout(
        title='Volcano Plot',
        xaxis_title='Log2 Fold Change',
        yaxis_title='-Log10 Q-value'
    )

    selection = st.plotly_chart(fig, on_select='rerun')
    selected_points = selection['selection']['point_indices']

    if len(selected_points) == 0:
        st.dataframe(pivot_table, use_container_width=True)

    else:
        st.dataframe(pivot_table.iloc[selected_points], use_container_width=True)
