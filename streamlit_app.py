import os
import zipfile

import numpy as np
import pandas as pd
import streamlit as st
import scipy.stats as stats
import plotly.graph_objs as go
import statsmodels.stats.multitest as smm
st.set_page_config(page_title='DIANN Output PCQ Converter')

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

    group1_files = st.multiselect('Select group 1 files', file_names)
    remaining_files = [file for file in file_names if file not in group1_files]
    group2_files = st.multiselect('Select group 2 files', remaining_files)

    if not group1_files or not group2_files:
        st.warning('Please select at least one file for each group')
        st.stop()

file_to_group = {file: 1 for file in group1_files}
file_to_group.update({file: 2 for file in group2_files})

# add sequential number to the dataframe
df['PSM.ID'] = range(1, len(df) + 1)
df['Channel.Score'] = 1 - df['Channel.Q.Value']
df['group'] = df['File.Name'].apply(lambda x: file_to_group.get(x))

# drop rows with missing group
df = df.dropna(subset=['group'])

# set psm id to be indexed
cols_to_keep = ['PSM.ID', 'Stripped.Sequence', 'Light/Heavy.Ratio', 'Quantity.Quality', 'Protein.Group']
group1_df = df[df['group'] == 1][cols_to_keep].set_index('PSM.ID')
group2_df = df[df['group'] == 2][cols_to_keep].set_index('PSM.ID')

t1, t2, t3, t4 = st.tabs(['Download', 'Group 1', 'Group 2', 'Volcano Plot'])
with t1:

    zip_name = st.text_input('Enter the name of the ZIP file', 'groups.zip')
    c1, c2 = st.columns(2)
    group1_name = c1.text_input('Enter the name of the group 1 file', 'group1.tsv')
    group2_name = c2.text_input('Enter the name of the group 2 file', 'group2.tsv')

    try:
        # Save the two dfs to files and then zip them together and have this as a download file (st.download_button)
        group1_df.to_csv(group1_name, sep='\t')
        group2_df.to_csv(group2_name, sep='\t')

        with zipfile.ZipFile(zip_name, 'w') as zipf:
            zipf.write(group1_name)
            zipf.write(group2_name)

        with open(zip_name, 'rb') as f:
            st.download_button('Download ZIP', f, file_name=zip_name, type='primary', use_container_width=True)
    except Exception as e:
        st.error(f'An error occurred: {e}')

    finally:
        # Cleanup the CSV files after zipping
        os.remove(group1_name)
        os.remove(group2_name)
        # Cleanup the ZIP file after download
        os.remove(zip_name)

with t2:
    st.subheader('Group 1')
    st.dataframe(group1_df, use_container_width=True)
with t3:
    st.subheader('Group 2')
    st.dataframe(group2_df, use_container_width=True)
with t4:
    st.subheader('Volcano Plot')
    c1, c2 = st.columns(2)
    index_col = c1.selectbox('Select index column', ['Modified.Sequence', 'Stripped.Sequence'], index=1)
    values_col = c2.selectbox('Select values column', ['Light/Heavy.Log2Ratio', 'Light/Heavy.Ratio'], index=0)
    agg_func = st.selectbox('Select aggregation function', ['mean', 'median'], index=0)

    # Handle missing values
    df = df.dropna(subset=[index_col, values_col, 'group'])

    pivot_table = df.pivot_table(index=index_col, columns='group', values=values_col, aggfunc=agg_func)

    # drop rows with missing values
    pivot_table = pivot_table.dropna()

    # Calculate p-values
    p_values = []
    for seq in pivot_table.index:
        group1 = df[(df[index_col] == seq) & (df['group'] == 1)][values_col].dropna()
        group2 = df[(df[index_col] == seq) & (df['group'] == 2)][values_col].dropna()
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


