import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import brewer2mpl as mpl
import statsmodels.api as sm
import math
from os.path import dirname, abspath
import pandas as pd
from sklearn.metrics import mean_squared_error, mean_absolute_error

from me_validation_configurations import cogs_to_skip
import me_validation_resources

resource_dir = dirname(abspath(me_validation_resources.__file__))

mpl_map = mpl.get_map('Set1', 'Qualitative', 9).mpl_colors
plt.style.use('ggplot')

COG_df = pd.read_csv('%s/cogs_ecoli_mg1655.csv' % resource_dir,
                     encoding="ISO-8859-1").set_index('locus')

skip = cogs_to_skip


def fit_line(x, y):
    """Return slope, intercept of best fit line."""
    X = sm.add_constant(x)
    model = sm.OLS(y, X, missing='drop')  # ignores entires where x or y is NaN
    fit = model.fit()
    return fit.params[1], fit.params[0]


def plot_grouped_by_cog(df, comparison_columns, savepath):
    i = 0
    fig, ax = plt.subplots(figsize=(5, 5))
    for key, values in df.iterrows():
        plt.plot([1e-11, 1], [1e-11, 1], linewidth=1, color='r')
        marker = 's'
        if i < 7:
            marker = 'o'
        elif i < 14:
            marker = '*'
        plt.loglog(values[comparison_columns[0]],
                   values[comparison_columns[1]], label=key,
                   markersize=12, marker=marker)
        plt.xlim([1e-7, 1])
        plt.ylim([1e-7, 1])
        i += 1

    plt.title('%s_%s by COG' % (comparison_columns[1], comparison_columns[0]))
    plt.ylabel(r'Simulated mole fraction', fontsize='xx-large')
    plt.xlabel(r'Measured mole fraction', fontsize='xx-large')

    fig.patch.set_facecolor('white')
    ax.set_axis_bgcolor('#f9f9f9')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()

    mse = mean_squared_error(df['Measured'], df['Simulated'])
    plt.annotate(r'$MSE$ = %.3f' % mse, (1e-6, 1e-2), fontsize='x-large')

    plt.savefig(savepath)


def format_datapoints_based_on_cog(dataframe):
    COG_format = {}
    for i, COG in enumerate(set(dataframe.COG.values)):

        if COG in skip:
            continue
        if type(COG) == float:
            COG = ''
        print(i, COG)
        COG_format[COG] = {}
        if i < 9:
            COG_format[COG]['Color'] = i
            COG_format[COG]['Marker'] = 's'
        elif i < 18:
            COG_format[COG]['Color'] = i - 9
            COG_format[COG]['Marker'] = 'o'
        elif i < 27:
            COG_format[COG]['Color'] = i - 18
            COG_format[COG]['Marker'] = '*'

    for COG in set(dataframe.COG.values):
        if COG in skip:
            continue
        if type(COG) == float:
            COG = ''
        dataframe.loc[dataframe.COG == COG, 'Color'] = COG_format[COG]['Color']
        dataframe.loc[dataframe.COG == COG, 'Marker'] = \
            COG_format[COG]['Marker']


def plot_pairwise_comparison(dataframe, comparison_columns, savepath,
                             only_metabolic=False, only_cytosolic=False):

    dataframe.COG = dataframe.COG.fillna('')
    if only_metabolic:
        dataframe = dataframe[dataframe.Metabolic == True]
    if only_cytosolic:
        dataframe = dataframe[dataframe.Membrane_complex_associated == False]

    format_datapoints_based_on_cog(dataframe)
    column1 = comparison_columns[0]
    column2 = comparison_columns[1]
    for column in comparison_columns:
        dataframe = dataframe[dataframe[column] > 1e-15]

    i = 0
    plt.figure()

    label_list = []
    plt.plot([1e-11, 1], [1e-11, 1], linewidth=1, color='r')
    min = 1
    for key, values in dataframe.iterrows():
        y = values[column2]
        x = values[column1]
        min = x if x < min else min
        min = y if y < min else min
        COG = values['COG']
        if COG in skip:
            continue
        if COG not in label_list:
            plt.loglog(x, y,
                       marker=values['Marker'], label=COG, alpha=0.66,
                       markersize=5, color=mpl_map[int(values['Color'])])
        else:
            plt.loglog(x, y,
                       marker=values['Marker'], alpha=0.66,
                       markersize=5, color=mpl_map[int(values['Color'])])
        i += 1
        label_list.append(COG)
    plt.xlim([min, 1])
    plt.ylim([min, 1])

    corr_df = dataframe[[comparison_columns[0], comparison_columns[1]]]
    r2 = corr_df.applymap(math.log10).corr().values[0][1] ** 2
    mse = mean_absolute_error(corr_df.applymap(math.log10)[column1],
                              corr_df.applymap(math.log10)[column2])
    plt.annotate(r'MSE = %.7f' % mse, (min*2, 2e-1), ha='left',
                 fontsize='xx-large')
    plt.annotate(r'$R^2$ = %.5f' % r2, (min*2, 2e-2),
                 ha='left',
                 fontsize='xx-large')
    plt.ylabel(r'Simulated mole fraction')
    plt.xlabel(r'Measured mole fraction')
    plt.title('%s_%s' % (column1, column2))
    lgd = plt.legend(bbox_to_anchor=(1.05, 1.0))

    plt.savefig(savepath, bbox_inches='tight', bbox_extra_artist=[lgd])

    return r2, mse
