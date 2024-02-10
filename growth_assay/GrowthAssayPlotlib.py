
EXTENSION = 'png'

YLABEL = 'Growth rate (PEMR)'

speed_order = ['Static', '45-rpm', '90-rpm', '180-rpm']
plate_order = ['50ml flask', '6-well', '24-well', '96-well']

colors_light = {
    'NF54': '#829AA3',
    'Heparin 50':'#829AA3',
    
    '3D7':'#B9BDC5',
    '3D7 (viola)': '#B9BDC5',

    'KOP230P': '#91b1be',
    'KOPfs25': '#9fc9d9',

    'KOEBA140': '#DAEBC6',
    'KOEBA181': '#84B06D',
    'KOEBA175': '#C1E2D4',
    'KOEBA175+Anti-GYPA':'#A5D0D4',
    'KORH1': '#F9DA9A',
    'KORH4': '#F3A59B',
    'KORh2a': '#F7CAA1',
    'KORH2a': '#F7CAA1',

    

    'Barseq': '#BEBFBC',

    'cKOMSP1 DMSO':'#B9BDC5',
    'cKOMSP1 Rap':'#FFE3B8',

    'Neuraminidase':'#D6EE96',
    
    'Anti-CD55':'#97B6DA',
    'Anti-GYPA':'#B6E2DF',
    'R1':'#F1AAD1',
    'Anti-GYPC':'#93E1DB',
    'Anti-CR1 ab25':'#DEAC7F',
    'Anti-CR1 ABIN':'#DEAC7F',

    'Anti-Basigin':'#C6C5D3',
    
    'cKOGAP45 DMSO': '#B9BDC5',
    'cKOGAP45 Rap': '#E5486C',

    'cKOAMA1 DMSO':'#B9BDC5',
    'cKOAMA1 Rap':'#EBD1D9',
    
 
}

colors_dark = {
    'NF54': '#043546',
    'Heparin 50':'#818897',

    '3D7':'#818897',
    '3D7 (viola)': '#818897',

    'KOP230P': '#406a7d',
    'KOPfs25': '#78a3b9',

    'KOEBA140': '#9ECA69',
    'KOEBA181': '#486737',
    'KOEBA175': '#6AB997',
    'KOEBA175+Anti-GYPA':'#04A0AF',
    'KORH1': '#F1AF23',
    'KORH4': '#E8513E',
    'KORh2a': '#ED8626',
    'KORH2a': '#ED8626',


    'Barseq': '#BEBFBC',

    'cKOMSP1 DMSO':'#FFD086',
    'cKOMSP1 Rap':'#FFD086',

    'Neuraminidase':'#ABDC29',
    
    'Anti-CD55':'#97B6DA',
    'Anti-GYPA':'#91D4CF',
    'R1':'#E459A5',
    'Anti-GYPC':'#35BEB4',
    'Anti-CR1 ab25':'#C6545E',
    'Anti-CR1 ABIN':'#C6545E',

    'Anti-Basigin':'#7E7C9C',
   
    
  
    'cKOGAP45 DMSO': '#81132C',
    'cKOGAP45 Rap': '#81132C',

    
    'cKOAMA1 DMSO':'#CC889D',
    'cKOAMA1 Rap':'#CC889D',
}

line_order = [
    'NF54', 
    'KOP230P',
    'KOPfs25', 
    'KOEBA140',
    'KOEBA175', 
    'KOEBA181', 
    'KORH1', 
    'KORH2a', 
    'KORH4',
]

gene_order = [
    'Pfcyp87', # generally not included
    
    'AMA1', # housekeeping genes
    'RH5',

    'EBA140', # EBA family
    'EBA175', 
    'EBA181', 

    'RH1', # RH family
    'RH2a', 
    'RH2b',
    'RH4', 

    'EBA165', # pseudo genes
    'RH3', 
]

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import natsort

def make_pairs(lst):
    # Use list comprehension with slicing
    return [lst[i:i + 2] for i in range(0, len(lst), 2)]

def load_file(file_name, column_title_row):
    print('\nLoading', file_name)
    df = pd.read_excel(file_name, header=column_title_row)
    df.columns = [col.strip() for col in df.columns]
    for col in df.select_dtypes(include=['object']):
        df[col] = df[col].str.strip()

    # compute repeats
    # Find all the 'Day' columns
    day_columns = [col for col in df if col.startswith('Day')]
    print(day_columns)
    # Calculate ratios for each consecutive day pair
    for i, (a,b) in enumerate(make_pairs(day_columns)):
        new_repeat_name = f'Repeat {i+1}'
        new_para_name = f'Parasitemia {i+1}'

        print(f'{new_repeat_name} = {b} / {a}, {new_para_name} = {a}')
        df[new_repeat_name] = df[b] / df[a]
        df[new_para_name] = df[a]

    
    df['Line'] = df['Line'].str.strip()
    df['Clone'] = df['Line'].str.split(' ').apply(lambda x: x[1] if len(x) > 1 and x[1].startswith('C') else '')
    df['Strain'] = df['Line'].str.split(' ').apply(lambda x: x[0])

    return df

def melt_df(df, id_vars):

    repeat_columns = [col for col in df if col.startswith('Repeat ')]
    parasitemia_columns = [col for col in df if col.startswith('Parasitemia ')]

    # Melt the DataFrame for Repeat columns
    df_melted_repeat = pd.melt(
        df,
        id_vars=id_vars,
        value_vars=repeat_columns,
        var_name='Repeat',
        value_name='Invasion rate'
    )

    # Melt the DataFrame for Parasitemia columns
    df_melted_parasitemia = pd.melt(
        df,
        id_vars=id_vars,
        value_vars=parasitemia_columns,
        var_name='Repeat',
        value_name='Parasitemia'
    )

    df_melted_parasitemia['Repeat'] = 'Repeat ' + df_melted_parasitemia['Repeat'].str.split(' ').str[1]

    # # Merge the two melted DataFrames based on common columns
    return pd.merge(df_melted_repeat, df_melted_parasitemia, on=id_vars + ['Repeat'])

def group_repeats(df_melted_replicas_filtered, group_keys):
    df_melted = df_melted_replicas_filtered.groupby(group_keys, as_index=False).agg({'Invasion rate': 'mean', 'Parasitemia': 'mean'})
    return df_melted.dropna()

def get_output_folder(experiment, folder):
    output_folder = os.path.join(folder, f'Graphs {experiment}')
    if not os.path.exists(output_folder): os.makedirs(output_folder)

    return output_folder


MAX_PARASITEMIA = 1
MAX_INVASION_RATE = 25
def filter_dataframe(df, ignore_repeats, filter_list_filename):    
    dfq = df.copy()
    dfq = dfq.query('Repeat != @ignore_repeats')
    
    df_removed = dfq.query('Parasitemia >= @MAX_PARASITEMIA or `Invasion rate` >= @MAX_INVASION_RATE')
    df_removed.to_csv(filter_list_filename)
    for i, r in df_removed.iterrows():
        repeat = r["Repeat"].split(" ")[1]
        day = 2*int(repeat)
        print(f'({day}, "{r["Well position"]}"),')

    dfq = dfq.query('Parasitemia < @MAX_PARASITEMIA')
    dfq = dfq.query('`Invasion rate` < @MAX_INVASION_RATE')
    return dfq

def mark_dataframe(df, ignore_repeats):    
    dfq = df.copy()
    dfq['IgnoreRepeat'] = dfq['Repeat'].apply(lambda x: x in ignore_repeats)
    dfq['AboveMaxParasitemia'] = dfq['Parasitemia'] >= MAX_PARASITEMIA
    dfq['AboveMaxInvasionRate'] = dfq['Parasitemia'] >= MAX_INVASION_RATE
    return dfq


from scipy import stats

# returns true if normally distributed with 5% significance level
def adtest(x):
    test = stats.anderson(x)
    return test.statistic < test.significance_level[2]

def test(a, b):
    # a is reference stample, wt if applicable

    normally_ditributed = adtest(a)

    if normally_ditributed and adtest(b):
        test_label = 'ttest p-value'
        p_val = stats.ttest_ind(a, b).pvalue

    else:
        test_label = 'ranksums'
        p_val = stats.ranksums(a, b).pvalue

    return f'{test_label} = {p_val:.4f}, {"significant" if p_val < 0.05 else "not"}'


from itertools import combinations

def significance_testing(df, group_keys, wildtype, line_key='Line'):
    
    # Output: 
    # - for lines, compare to wildtype (for each speed and plate)
    # - for speeds, compare to each other
    # - for plates, compare to each other
    # - for hematocrit, compare to each other

    for k in group_keys:
        remaining_keys = [g for g in group_keys if g!= k]
        print()

        for conditions, dfg in df.groupby(remaining_keys):

            k_values = dfg[k].unique()
            if len(k_values) == 1: 
                print(f'Only one value of {k} for {", ".join(conditions)}')

            print(f'Testing {k} for {", ".join(conditions)}.')
            
            if k == line_key:
                k1 = wildtype
                v1 = dfg.query(f'{line_key} == @k1')['Invasion rate']
                if len(v1.index) == 0:
                    print(f'   - No match for wildtype {k1}')
                    continue
                for k2 in [v for v in k_values if v != k1]:
                    v2 = dfg.query(f'{line_key} == @k2')['Invasion rate']
                    print(f'   - {k1} and {k2:12s}: {test(v1, v2)}')
                
            else:
                for k1, k2 in combinations(k_values, 2):
                    v1 = dfg.query(f'{k} == @k1')['Invasion rate']
                    v2 = dfg.query(f'{k} == @k2')['Invasion rate'] 
                    print(f'   - {k1:8s} and {k2:8s}: {test(v1, v2)}')

        # # line, speed, plate

        # # for line, speed, test plate:
        # # for line, plate, test speed:
        # # for speed, plate, test line:



# Plots

## Plot parasimetia vs invasion rate 
def plot_parasimetia_vs_invasion_rate(df_melted_replicas, group_keys, output_folder):
    for (line, speed), dfg in df_melted_replicas.groupby(group_keys):
        plt.figure(figsize=(6,4), dpi=300)
        sns.scatterplot(
            data=dfg,
            x='Parasitemia',
            y='Invasion rate',
            hue='Repeat',
            palette='rocket',
        )
        plt.title(f'GA2 {line} {speed}')
        plt.savefig(os.path.join(output_folder, f'Parasetimia vs invasion rate {line} {speed}.{EXTENSION}'), bbox_inches='tight')
        plt.close()


def plot_line(df_melted, label, max_y, output_folder):
    for line, dfg in df_melted.groupby('Line'):
        plt.figure(figsize=(6,4), dpi=300)
        sns.boxplot(x='Speed', y='Invasion rate', data=dfg, order=speed_order, palette='rocket')
        sns.swarmplot(x='Speed', y='Invasion rate', data=dfg, order=speed_order, color='black', size=5)
        plt.title(line)
        plt.ylim((0,max_y))
        plt.savefig(os.path.join(output_folder, f'Line {line} {label}.{EXTENSION}'), bbox_inches='tight')
        plt.close()

def plot_speed(df_melted, label, max_y, output_folder):
    for speed, dfg in df_melted.groupby('Speed'):
        plt.figure(figsize=(6,4), dpi=300)
        lines = dfg['Line'].unique()
        line_order = natsort.natsorted(lines)
        sns.boxplot(x='Line', y='Invasion rate', data=dfg, order=line_order, palette='rocket')
        sns.swarmplot(x='Line', y='Invasion rate', data=dfg, order=line_order, color='black', size=5)
        plt.title(speed)
        # rotate x labels
        plt.xticks(rotation=90)
        plt.ylim((0,max_y))
        plt.savefig(os.path.join(output_folder, f'Speed {speed} {label}.{EXTENSION}'), bbox_inches='tight')
        plt.close()

def plot_plate(df_melted, label, max_y, output_folder):
    for plate, dfg in df_melted.groupby('plate'):
        plt.figure(figsize=(6,4), dpi=300)
        sns.boxplot(x='Plate', y='Invasion rate', data=dfg, order=plate_order, palette='rocket')
        sns.swarmplot(x='Plate', y='Invasion rate', data=dfg, order=plate_order, color='black', size=5)
        plt.title(plate)
        # rotate x labels
        plt.xticks(rotation=90)
        plt.ylim((0,max_y))
        plt.savefig(os.path.join(output_folder, f'Plate {plate} {label}.{EXTENSION}'), bbox_inches='tight')
        plt.close()


def plot_repeats_for_line(df_melted_replicas_filtered, max_y, output_folder):
    for line, dfg in df_melted_replicas_filtered.groupby('Line'):
        plt.figure(figsize=(6,4), dpi=300)
        sns.lineplot(x='Repeat', y='Invasion rate', data=dfg, hue='Speed', palette='rocket', marker='o')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title(line)
        plt.ylim((0,max_y))
        plt.savefig(os.path.join(output_folder, f'Repeats for line {line}.{EXTENSION}'), bbox_inches='tight')
        plt.close()

# Plot how each well behaves per replica
def plot_repeats_for_wells(df_melted_replicas_filtered, max_y, output_folder):
    for (line, speed), dfg in df_melted_replicas_filtered.groupby(['Line', 'Speed']):
        plt.figure(figsize=(6,4), dpi=300)
        sns.lineplot(x='Repeat', y='Invasion rate', hue='Well position', data=dfg, palette='rocket', marker='o')
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title(f'{line} {speed}')
        plt.ylim((0,max_y))
        plt.savefig(os.path.join(output_folder, f'Repeats for wells {line} {speed}.{EXTENSION}'), bbox_inches='tight')
        plt.close()
        



strain_order = ['NF54', 'KOP230P', 'KOPfs25', 'KOEBA140', 'KOEBA175', 'KOEBA181', 'KORH1', 'KORH2a', 'KORH4', 'Barseq']
line_order = ['NF54', 'KOP230P C3','KOP230P C5', 'KOPfs25 C1', 'KOPfs25 C3', 'KOEBA140 C3','KOEBA140 C4', 'KOEBA175 C6', 'KOEBA181 C1', 'KOEBA181 C2', 'KORH1 C1' , 'KORH2a C1', 'KORH2a C3', 'KORH4 C1', 'Barseq']


def plot_pretty_barplot(df, plot_clones, output_folder, bar_width=0.2):

    xlabel = "Line" if plot_clones else "Strain"
    x_sort_label = line_order if plot_clones else strain_order

    x_values = df[xlabel].unique()
    misisng_labels = [l for l in x_values if l not in x_sort_label]
    x_sort_label = [l for l in x_sort_label if l in x_values]

    if len(misisng_labels) > 0:
        print(f'Warning: missing order for {xlabel} for {misisng_labels}')

    for colors, (condition, dfq) in zip([colors_light, colors_dark], df.groupby('Speed')):

        print(condition)

        plt.figure()

        sns.set_theme(#context='notebook', 
                    #style='ticks', 
                    style='darkgrid', 
                    font='arial', 
                    font_scale=2) 
                    #   color_codes=True, 
                    #   rc=None)
        sns.set_style(rc = {'axes.facecolor': '#F5F5F9'})

        ax = sns.boxplot(
            data=dfq, 
            x=xlabel,
            order=x_sort_label, 
            y='Invasion rate', 
            hue='Strain', 
            palette=colors, 
            #saturation=.5, 
            dodge=False, 
            showfliers=False, 
            width=bar_width,
        )
        
        #sns.swarmplot(data=df, x='Strain', y='Mean', order=x_sort_label)
        plt.xticks(rotation=-90)
        #sns.stripplot(data=dfq, x='Strain', order=x_sort_label, y='Mean', marker='o', color='.0', alpha=0.3)
        
        sns.stripplot(
            data=dfq, 
            x=xlabel, 
            y='Invasion rate', 
            order=x_sort_label, 
            hue='Strain', 
            jitter=True, 
            linewidth=0.5, 
            palette=colors
        )
    
        plt.gcf().set_size_inches(12,9)
        plt.gcf().set_dpi(200)

        plt.ylim((0,14))
        plt.legend([],[], frameon=False)
        
        plt.xlabel('')
        plt.ylabel(YLABEL)

        plt.rcParams['svg.fonttype'] = 'none'
        plt.savefig(f'{output_folder}/FIG Invasion rate separate ({condition}).svg', bbox_inches='tight')


def ga1_plot_together(df, plot_clones, output_folder):

    x_label = "Line" if plot_clones else "Strain"
    x_sort_label = line_order if plot_clones else strain_order

    x_values = df[x_label].unique()
    misisng_labels = [l for l in x_values if l not in x_sort_label]
    x_sort_label = [l for l in x_sort_label if l in x_values]

    if len(misisng_labels) > 0:
        print(f'Warning: missing order for {x_label} for {misisng_labels}')


    # Separate box plots for shaking and static conditions
    plt.figure(figsize=(15, 12))

    OFFSET = 0.152
    #for boxes to sit next to each other offset should be half the width of the boxes

    # Shaking condition
    sns.set_theme(#context='notebook', 
                #style='ticks', 
                style='darkgrid', 
                font='arial', 
                font_scale=3.5) 
                #   color_codes=True, 
                #   rc=None)
    sns.set_style(rc = {'axes.facecolor': '#F5F5F9'})
    #sns.set_style(rc = {'axes.facecolor': '#FBF9F6'})

    edge_color='black'

    PROPS = {
        'boxprops':{'edgecolor':edge_color},
        'medianprops':{'color':'red'},
        'whiskerprops':{'color':edge_color},
        'capprops':{'color':edge_color}
    }

    ax_shaking = sns.boxplot(
        data=df[df['Speed'] == '90-rpm'],
        x=x_label,
        order=x_sort_label,
        y='Invasion rate',
        hue='Strain',
        palette=colors_dark,
        showfliers=False,
        dodge=False, 
        width=0.3,
        box_offset=OFFSET,
        **PROPS,
    )

    # Static condition
    sns.set_theme(context='notebook', style='darkgrid', font='Candara', font_scale=2, color_codes=True, rc=None)
    #sns.set_style(rc = {'axes.facecolor': '#FBF9F6'})
    sns.set_style(rc = {'axes.facecolor': '#F5F5F9'})
    edge_color='grey'
    PROPS = {
        'boxprops':{'edgecolor':edge_color},
        'medianprops':{'color':'red'},
        'whiskerprops':{'color':edge_color},
        'capprops':{'color':edge_color}
    }
    ax_static = sns.boxplot(
        data=df[df['Speed'] == 'Static'],
        x=x_label,
        order=x_sort_label,
        y='Invasion rate',
        hue='Strain',
        palette=colors_light,
        dodge=False, 
        showfliers=False,
        width=0.3,
        box_offset=-OFFSET,
        **PROPS,
    )

    plt.xticks(rotation=-90)
    plt.ylim((0, 14))
    plt.legend([], [], frameon=False)
    plt.xlabel(' ')
    plt.ylabel('Mean invasion rate')

    plt.tight_layout()
    plt.savefig(f'{output_folder}/invasion_rate_separate.pdf')

    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(f'{output_folder}/growthassay_fig3c.svg', bbox_inches='tight')