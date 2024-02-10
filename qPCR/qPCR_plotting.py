
import json 
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as tick
import sigfig
from matplotlib.ticker import ScalarFormatter, FuncFormatter
from EmmaPlotColors import *
import scipy.optimize
    

import warnings
# Suppress all warnings
warnings.filterwarnings("ignore")

OUTPUT_FOLDER = 'Graphs'

def get_output_folder(experiment):
    output_folder = os.path.join(OUTPUT_FOLDER, experiment)
    if not os.path.exists(output_folder): 
        os.makedirs(output_folder)
    return output_folder


from sklearn.linear_model import HuberRegressor
from sklearn.metrics import r2_score

def linear_regression(xs, ys):

    xs = np.array(xs)
    ys = np.array(ys)

    # Initialize the HuberRegressor model
    model = HuberRegressor()

    # Fit the model to the data
    model.fit(xs.reshape(-1, 1), ys)

    # Get the slope and intercept of the linear regression
    slope = model.coef_[0]
    intercept = model.intercept_

    # Predict the values using the fitted model
    predicted_ys = model.predict(xs.reshape(-1, 1))

    # Calculate the R-squared value
    r_squared = r2_score(ys, predicted_ys)

    return slope, intercept, r_squared







# from sklearn.linear_model import HuberRegressor

# def linear_regression(xs, ys):

#     model = HuberRegressor
#     results = evaluate_model(xs, ys, model)

#     # todo

#     return slope, intercept, r_squared

# def linear_regression(xs, ys):
#     linear_function = lambda x, slope, intercept: slope * x + intercept
#     # Perform linear fit using Scipy's curve_fit
#     popt, _ = scipy.optimize.curve_fit(linear_function, xs, ys)
#     slope, intercept = popt
#     fit_line = linear_function(xs, slope, intercept)

#     # Calculate R-squared
#     residuals = ys - fit_line
#     ss_res = np.sum(residuals**2)
#     ss_tot = np.sum((ys - np.mean(ys))**2)
#     r_squared = 1 - (ss_res / ss_tot)

#     # Calculate linear fit coefficients using log base 10
#     # fit_coefficients = np.polyfit(np.log10(dfg['sample']), dfg['Cq'], 1)
#     # slope = fit_coefficients[0]
#     # intercept = fit_coefficients[1]
    
#     # Create the fit line
#     # fit_line = np.polyval(fit_coefficients, np.log10(dfg['sample']))
    
#     # # Calculate R-squared
#     # correlation_matrix = np.corrcoef(np.log10(dfg['sample']), dfg['Cq'])
#     # r_squared = correlation_matrix[0, 1] ** 2
    
#     return slope, intercept, r_squared
    

def find_outlier(points, threshold):
    num_points = len(points)
    points = np.array(points)

    max_outlier = np.abs(np.max(points - np.mean(points)))

    if num_points < 3: return None
    if max_outlier < threshold: return None

    # Calculate pairwise Euclidean distances
    distances = np.zeros((num_points, num_points))
    for i in range(num_points):
        for j in range(i + 1, num_points):
            distances[i][j] = abs(points[i] - points[j])
            distances[j][i] = distances[i][j]

    # Find the point with the largest sum of distances to other points
    sum_distances = np.sum(distances, axis=1)
    outlier_index = np.argmax(sum_distances)

    return points[outlier_index]

def add_outlier_column(df, threshold):
    for _, dfg in df.groupby(['sample', 'gene', 'type']):
        Cqs = dfg['Cq']
        max_value = find_outlier(Cqs, threshold)
        df.loc[dfg.index, 'is_outlier'] = False if max_value == None else (Cqs == max_value) 


def fit_to_standard(df, plot=True, output_folder=None):

    if plot:
        sns.set_theme(context='notebook', 
                    style='darkgrid', 
                    #font='Candara', 
                    font_scale=1, 
                    color_codes=True, 
                    rc=None)
        sns.set_style(rc = {'axes.facecolor': '#F5F5F9'})

    # Plot standards
    dfq = df.query('type == "standard"')

    standard_results = {}

    def efficiency_for_slope(slope):
        E = 10**(-1/slope)
        return (E - 1) * 100

    for gene, dfg in dfq.groupby('gene'):

        dfg = pd.DataFrame(dfg)
        dfg['sample'] = dfg['sample'].astype(float)

        # print(f'Removed {len(dfg.loc[dfg["is_outlier"]].index)} outliers from {gene}')
        dfg = dfg.query('is_outlier == False')
        # dfg = dfg.loc[~dfg['is_outlier']]
        
        dfg = dfg.dropna(subset=['Cq'])

        slope, intercept, r_squared = linear_regression(np.log10(dfg['sample']), dfg['Cq'])
        
        # compute efficiency
        efficiency = efficiency_for_slope(slope)

        # Create the fit line
        fit_line = np.polyval((slope, intercept), np.log10(dfg['sample'])) 

        label = f'Fit: y = {slope:.2f} * log10(x) + {intercept:.2f}\nR-squared = {r_squared:.2f}\nEfficiency = {efficiency:.0f}%'
        standard_results[gene] = {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared,
            'efficiency': efficiency,
        }

        if plot:
            # Scatter plot
            plt.figure(figsize=(6,4), dpi=300)
            sns.scatterplot(data=dfg, x='sample', y='Cq')
            # Add the fit line and R-squared to the plot
            plt.plot(dfg['sample'], fit_line, color='red', label=label)
            
            title = f'Linear Fit and Scatterplot for {gene}'

            plt.xlabel('Concentration')
            plt.xscale('log')  # Set x-axis to log scale
            plt.ylabel('Cq')
            plt.title(title)
            plt.legend()
            plt.grid(True)
            plt.ylim((15,30))
            plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
            plt.tight_layout()

            if output_folder:
                plt.savefig(os.path.join(output_folder, f'{title}.png'))
    
    return standard_results


def plot_nort_controls(df, output_folder=None):
    dfq = df.query('type == ["cDNA", "NO-RT control"]')

    for sample, dfg in dfq.groupby('sample'): 
        
        sns.set_theme(context='notebook', 
                    style='darkgrid', 
                    #font='Candara', 
                    font_scale=1, 
                    color_codes=True, 
                    rc=None)
        sns.set_style(rc = {'axes.facecolor': '#F5F5F9'})

        plt.figure(figsize=(6,4), dpi=300)
        sns.stripplot(
            data=dfg,
            x='gene',
            y='Cq',
            hue='type_outlier',
        )
        title = f'No-RT control for {sample}'
        plt.title(sample)
        plt.xticks(rotation=20)
        plt.ylim((15,40))
        
        legend = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()

        if output_folder:
            plt.savefig(
                os.path.join(output_folder, f'{title}.png'),
                bbox_extra_artists=[legend],
            )

def make_fold_change_plot(df, x_key, y_key, hue_key, title, ylabel, output_folder, figsize=(12,10), errorbar=None):

    sns.set_theme(context='notebook', 
                style='darkgrid', 
                #font='Candara', 
                font_scale=2, 
                color_codes=True, 
                rc=None)
    sns.set_style(rc = {'axes.facecolor': '#F5F5F9'})

    plt.figure(figsize=figsize, dpi=300)
    ax = sns.barplot(
        y=df[y_key]-1, 
        x=df[x_key], 
        hue=df[hue_key], 
        errorbar=errorbar,
        bottom=1, 
        palette=colors_dark
    )
    ax.axhline(y=1, c='k', linestyle='--', lw=0.8) # Make the line dashed and thinner

    plt.yscale('log', base=2)
    plt.ylim((0.0001,60))
    # Improved y-tick formatting
    #major_ticks = np.append(np.arange(0.01, 0.8, 0.01), np.arange(1, 2, 0.1))
    #minor_ticks = np.arange(0.3, 1, 0.1)
    #ax.set_yticks(major_ticks)
    #ax.set_yticks(minor_ticks, minor=True)
    #ax.yaxis.set_major_formatter(lambda x, pos: f'{x:.0f}' if x >= 1 else f'{x:.1f}')

    ax.yaxis.set_major_formatter(tick.FuncFormatter(lambda x, y: sigfig.round(x, sigfigs=1, type=str)))
    plt.xticks(rotation=45)

    # Remove the outline of the plot
    sns.despine()

    # Add title and labels for clarity
    ax.set_title(title)
    ax.set_xlabel("Gene")
    ax.set_ylabel(ylabel)

    legend = plt.legend(title='Sample', loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    
    if output_folder:
        plt.savefig(
            os.path.join(output_folder, f'{title}.png'),
            bbox_extra_artists=[legend],
        )



def load_dataframe(experiment_name):

    input_file = f'{experiment_name}.json'

    with open(input_file, 'r') as json_file:
        data = json.load(json_file)

    MM_actions = data['MM_actions']
    S_actions = data['S_actions']
    S_type = data['S_type']
    S_location = data['S_location']
    MM_location = data['MM_location']
    sample_number = data['sample_number']
    repeat = data['repeat']
    cq_file= data['cq_file']

    # This links the names of the wells to the master mixes that were placed in them
    master_mix_mapping = { destination_well: MM_location[source_well] for source_well, destination_wells in MM_actions.items() for destination_well in destination_wells }

    # This links the names of the wells to the samples that were placed in them
    sample_mapping = { destination_well: S_location[source_well] for source_well, destination_wells in S_actions.items() for destination_well in destination_wells }

    # This links the names of the wells to the samples type that were placed in them
    type_mapping = { destination_well: S_type[source_well[0]] for source_well, destination_wells in S_actions.items() for destination_well in destination_wells }
    
    if '165' in experiment_name:
        def s_tyepe_with_mixup(source_well):
            if source_well == 'A5' or source_well == 'A6': return "NO-RT control"
            if source_well == 'B5' or source_well == 'B6': return "cDNA"
            return S_type[source_well[0]]
        
        type_mapping = { destination_well: s_tyepe_with_mixup(source_well) for source_well, destination_wells in S_actions.items() for destination_well in destination_wells }


    # opens the file with the data in it and selects the columns with the relavent information in it

    if not os.path.exists(cq_file):
        raise Exception(f"The file '{cq_file}' does not exist.")
    else:
        df_full = pd.read_excel(cq_file)

    df = pd.DataFrame(df_full[['Well', 'Cq']])
    #This reformates the well names to remove the 0 so they match the formate of the lists above 
    df['Well'] = df['Well'].str[0] + df['Well'].str[1:].astype(str).str.lstrip('0')

    # This adds columes for the master mix, type and sample 
    df['gene'] = df['Well'].map(master_mix_mapping)
    df['type'] = df['Well'].map(type_mapping)
    df['sample'] = df['Well'].map(sample_mapping)
    df['sample_number'] = sample_number
    df['repeat'] = repeat
    df['experiment'] = experiment_name

    df.dropna(inplace=True)
    add_outlier_column(df, threshold=0.3)
    df['type_outlier'] = df['type'].str.cat(df['is_outlier'].astype(str), sep='_')

    return df


def add_concentration_and_normalizations(df, standard_results, normalization_gene, normalization_sample):
    df_cdna_replicates = pd.DataFrame(df.query('type == "cDNA"'))
    # compute ceoncentration from Cq
    for gene, dfg in df_cdna_replicates.groupby('gene'):
        standard_result = standard_results[gene]
        slope = standard_result['slope']
        intercept = standard_result['intercept']
        # print(gene, slope, intercept, len(dfg.index))

        df_cdna_replicates.loc[dfg.index, 'relative_concentration'] = 10 ** ((dfg['Cq'] - intercept) / slope)

    # compute mean and std for each combination of gene and sample
    df_cdna = df_cdna_replicates.groupby(['gene', 'sample'], as_index=False).agg(
        gene=('gene', 'first'),
        line=('sample', 'first'),
        mean_relative_concentration=('relative_concentration', 'mean'),
        std_relative_concentration=('relative_concentration', 'std'),
        sample=('repeat', 'first'),
        plate=('sample_number', 'first'),
    )
    # df_cdna = df_cdna_replicates.groupby(['gene', 'sample'], as_index=False).agg(
    #     gene=('gene', 'first'),
    #     line=('sample', 'first'),
    #     mean_relative_concentration=('relative_concentration', 'mean'),
    #     std_relative_concentration=('relative_concentration', 'std'),
    #     sample=('repeat', 'first'),
    #     plate=('sample_number', 'first'),
    # )

    gene_normalization_title = f'Normalized relative concentration ({normalization_gene})'
    sample_normalization_title = f'Expression fold change to {normalization_sample} (normalised to {normalization_gene})'

    # normalize by housekeeping genes
    for _, dfg in df_cdna.groupby('line'):
        gene_concentration = float(dfg.query('gene == @normalization_gene')['mean_relative_concentration'].iloc[0]) # should be just one before iloc[0]
        df_cdna.loc[dfg.index, gene_normalization_title] = dfg['mean_relative_concentration'] / gene_concentration

    # relative difference compared to NF54 (WT)
    if normalization_sample != None:
        for _, dfg in df_cdna.groupby('gene'):
            nf_concentration = float(dfg.query('line == @normalization_sample')[gene_normalization_title].iloc[0])
            df_cdna.loc[dfg.index, sample_normalization_title] = dfg[gene_normalization_title] / nf_concentration

    return df_cdna, sample_normalization_title, gene_normalization_title


def analse_qpcr_experiment(experiment_name):

    output_folder = get_output_folder(experiment=experiment_name)

    df = load_dataframe(experiment_name=experiment_name)
   
    standard_results = fit_to_standard(df, output_folder=output_folder)
    plot_nort_controls(df, output_folder=output_folder)

    normalization_genes = ['Actin1', 'Pfcyp87']
    normalization_samples = ['NF54', 'KOP230P']

    all_genes = df['gene'].unique()
    all_samples = df['sample'].unique()

    normalization_genes = [n for n in normalization_genes if n in all_genes]
    normalization_samples = [n for n in normalization_samples if n in all_samples]
    all_samples = [n for n in all_samples if n in all_samples]

    for normalization_gene, normalization_sample in [(g,s) for g in normalization_genes for s in normalization_samples]:
        
        print(f'Analysing for {normalization_gene} + {normalization_sample}')

        df_cdna, sample_normalization_title, gene_normalization_title = add_concentration_and_normalizations(
            df=df,
            standard_results=standard_results,
            normalization_gene=normalization_gene,
            normalization_sample=normalization_sample,
        )

        df_cdna_q = df_cdna.query('line != @normalization_sample and gene != @normalization_gene')
        
        make_fold_change_plot(
            df=df_cdna_q,
            x_key='gene',
            y_key=sample_normalization_title,
            title=f"Gene Expression Fold Change to {normalization_sample} ({normalization_gene})",
            ylabel=f"Fold Change (normalized to {normalization_gene})",
            hue_key='line',
            output_folder=output_folder,
        )

    return df

# group_ratio is ratio of space used for groups vs spacing of groups
# bar_ratio is ratio of space used for bars vs spacing of bars
def log_barplot(df, x_key, y_key, hue_key, hue_order, x_order, colors, ylabel, output_folder, title='', figsize=(20,10), group_ratio = 0.8, bar_ratio = 0.8, font_scale=2):

    categories = df[x_key].unique()
    subcategories = df[hue_key].unique()
    
    if x_order: categories = [l for l in x_order if l in categories]
    if hue_order: subcategories = [l for l in hue_order if l in subcategories]
    
    N = len(subcategories)
    group_width = group_ratio
    bar_width = group_width * bar_ratio / N
    bar_spacing = 0 if N == 1 else group_width * (1-bar_ratio) / (N-1)

    category_offset = [(i-(N-1)/2)*(bar_width+bar_spacing) for i in range(N)]

    values = np.array([[df[(df[x_key] == category) & (df[hue_key] == sub_category)][y_key].mean() for sub_category in subcategories] for category in categories])
    errors = np.array([[df[(df[x_key] == category) & (df[hue_key] == sub_category)][y_key].std() for sub_category in subcategories] for category in categories])

    x = np.arange(len(categories))  # the label locations

    # Define your color dictionary

    # Step 2: Create a bar plot

    plt.figure(figsize=figsize, dpi=200)
    sns.set_theme(context='notebook', 
              #style='darkgrid', 
              font='arial', 
              font_scale=font_scale, 
              color_codes=True, 
              rc=None)
    sns.set_style(rc = {'axes.facecolor': '#F5F5F9'})
    
    ax = plt.gca()

    for i, subcat in enumerate(subcategories):
        vs = values[:, i]
        es = errors[:, i]

        tops = np.where(vs > 1, vs-1, 1-vs)
        bottoms = np.where(vs > 1, 1, vs)
        error_bar_positions = np.where(vs > 1, vs, bottoms)

        yerr_down = es/2
        yerr_up = es

        xpos = x + category_offset[i]
        
        ax.bar(xpos, tops, bar_width, label=subcat, capsize=5, color=colors[subcat], bottom=bottoms)
        ax.errorbar(xpos, error_bar_positions, yerr=(yerr_down, yerr_up), fmt='none', color='black', capsize=2, label=None, capthick=0.6, elinewidth=0.6)

        # Translate error bars up or down based on 'tops' and 'bottoms'

    # Step 3: Customize the plot
    ax.set_xlabel(x_key.capitalize())
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    if len(subcategories) > 1:
        ax.legend(title=hue_key.capitalize(), loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_yscale('log', base=2)
    plt.xticks(rotation=90)
    ax.axhline(y=1, color='black', linewidth=0.8)
    
    ax.set_ylim(0.004, 64)
    
    y_formatter = ScalarFormatter()
    y_formatter.set_scientific(False)  # Disable scientific notation
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x,y: sigfig.round(x, 2, type=str)))

    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(os.path.join(output_folder, title + '.svg'), bbox_inches='tight')
    plt.savefig(os.path.join(output_folder, title + '.png'), bbox_inches='tight')
