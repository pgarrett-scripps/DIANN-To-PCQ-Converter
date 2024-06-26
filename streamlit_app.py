from collections import Counter
import os
import zipfile

import numpy as np
import pandas as pd
import streamlit as st
import scipy.stats as stats
import plotly.graph_objs as go
import statsmodels.stats.multitest as smm
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
cols_to_keep = ['PSM.ID', 'Stripped.Sequence', 'Light/Heavy.Ratio', 'Quantity.Quality', 'Protein.Group']

group_dfs = {}
for i in range(1, num_groups+1):
    group_dfs[i] = df[df['group'] == i][cols_to_keep].set_index('PSM.ID')

t1, t2, t3 = st.tabs(['Download', 'Groups', 'Volcano Plot'])
with t1:

    zip_name = st.text_input('Enter the name of the ZIP file', 'groups.zip')
    if group_type == 'Manual Selection':
        group_names = {}
        for i in range(1, num_groups+1):
            group_names[i] = st.text_input(f'Enter the name of the group {i} file', f'group{i}.tsv')
    else:
        group_names = {}
        for i in range(1, num_groups+1):
            fname = os.path.basename(file_names[i-1])
            fname = ''.join(fname.split('.')[0].split('/')[-1].split('\\')[-1])
            group_names[i] = st.text_input(f'Enter the name of the group {i} file', f'{fname}.tsv')
    
    try:
        for i in range(1, num_groups+1):
            group_dfs[i].to_csv(group_names[i], sep='\t')

        with zipfile.ZipFile(zip_name, 'w') as zipf:
            for i in range(1, num_groups+1):
                zipf.write(group_names[i])

        with open(zip_name, 'rb') as f:
            st.download_button('Download ZIP', f, file_name=zip_name, type='primary', use_container_width=True)
    except Exception as e:
        st.error(f'An error occurred: {e}')

    finally:
        for i in range(1, num_groups+1):
            os.remove(group_names[i])

        # Cleanup the ZIP file after download
        os.remove(zip_name)

with t2:
    for i in range(1, num_groups+1):
        st.subheader(f'Group {i}')
        for fname, j in file_to_group.items():
            if i == j:
                st.caption(fname)
        st.dataframe(group_dfs[i], use_container_width=True)
with t3:
    st.subheader('Volcano Plot')
    c1, c2 = st.columns(2)
    index_col = c1.selectbox('Select index column', ['Modified.Sequence', 'Stripped.Sequence'], index=1)
    values_col = c2.selectbox('Select values column', ['Light/Heavy.Log2Ratio', 'Light/Heavy.Ratio'], index=0)
    agg_func = st.selectbox('Select aggregation function', ['mean', 'median'], index=0)

    group_options = list(range(1, num_groups+1))
    group1_index = int(st.selectbox(label='Select First Group', options=group_options))
    group2_index = int(st.selectbox(label='Select Second Group', options=group_options))

    if not st.button('run'):
        st.stop()

    # Handle missing values
    df = df.dropna(subset=[index_col, values_col, 'group'])

    pivot_table = df.pivot_table(index=index_col, columns='group', values=values_col, aggfunc=agg_func)

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

    # if value is None replace with 1
    pivot_table = pivot_table.fillna(1)

    # Calculate q-values (FDR correction)
    pivot_table['q_value'] = smm.multipletests(pivot_table['p_value'], method='fdr_bh')[1]

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


