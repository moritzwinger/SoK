import math
from typing import List
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from brokenaxes import brokenaxes
from matplotlib.ticker import FuncFormatter


def plot(labels: List[str], pandas_dataframes: List[pd.DataFrame], fig=None) -> plt.Figure:
    """

    :param labels:
    :param pandas_dataframes:
    :param fig:
    :return:
    """
    # Save current figure to restore later
    previous_figure = plt.gcf()

    # Set the current figure to fig
    # figsize = (int(len(labels) * 0.95), 6)

    inches_per_pt = 1.0 / 72.27 * 2  # Convert pt to inches
    golden_mean = ((np.math.sqrt(5) - 1.0) / 2.0) * .8  # Aesthetic ratio
    fig_width = 300 * inches_per_pt  # width in inches
    # fig_height = (fig_width * golden_mean)  # height in inches
    fig_height = 2.15
    figsize = [fig_width * 0.67, fig_height / 1.22]

    # figsize = (5.5, 4)
    config_dpi = 100
    if fig is None:
        fig = plt.figure(figsize=figsize, dpi=config_dpi)
    else:
        plt.rcParams["figure.figsize"] = figsize
        plt.rcParams["figure.dpi"] = config_dpi

    plt.figure(fig.number)

    # change size and DPI of resulting figure
    # change legend font size
    plt.rcParams["legend.fontsize"] = 8
    # NOTE: Enabling this requires latex to be installed on the Github actions runner
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = 'serif'

    # TODO: For the purpose of quickly producing a plot, I have added EVA MLP and LeNet-5 as a single bar
    #  Later, however, we need to adjust the plotting code to handle grouped bars (see, for example, cardio)

    # Nice names for labels: maps folder name -> short name
    # and adds linebreaks where required
    subs = {'SEAL-CKKS-Batched': 'SEAL-CKKS\n(MLP)',
            'SEALion': 'SEALION\n(MLP)',
            'nGraph-HE-MLP': 'nGraph-HE\n(MLP)',
            'nGraph-HE-Cryptonets': 'nGraph-HE\n(Cryptonets)',
            'nGraph-HE-LeNet5': 'nGraph-HE\n(LeNet-5)',
            'EVA-MLP': 'EVA\n(MLP, LeNet-5)',  # TODO: This is just a temporary hack
            }
    labels = [subs.get(item, item) for item in labels]

    positions = {
        'SEAL-CKKS\n(MLP)': 0,
        'SEALION\n(MLP)': 1,
        'nGraph-HE\n(MLP)': 2,
        'nGraph-HE\n(Cryptonets)': 3,
        'nGraph-HE\n(LeNet-5)': 4,
        'EVA\n(MLP, LeNet-5)': 5  # TODO: This is just a temporary hack
    }

    sorted_labels = ('SEAL-CKKS\n{\\fontsize{7pt}{3em}\\selectfont{}(MLP)}',
                     'SEALion\n{\\fontsize{7pt}{3em}\\selectfont{}(MLP)}',
                     'nGraph-HE\n{\\fontsize{7pt}{3em}\\selectfont{}(MLP)}',
                     'nGraph-HE\n{\\fontsize{7pt}{3em}\\selectfont{}(Cryptonets)}',
                     'nGraph-HE\n{\\fontsize{7pt}{3em}\\selectfont{}(LeNet-5)}',
                     'EVA\n{\\fontsize{7pt}{3em}\\selectfont{}(MLP, LeNet-5)}'  # TODO: This is just a temporary hack
                     )

    # Setup brokenaxes
    # TODO: Make break depend on value of nGraph-HE LeNet-5 runtime?
    # hspace controls how much space is in between the broken axes left=0.15
    # bax = brokenaxes(ylims=((0, 10), (120, 130)), hspace=.3, despine=False, left = 0.25, bottom = 0.25)
    bax = brokenaxes(ylims=((0, 10), (120, 140)), hspace=.4, despine=False, left=0.115, bottom=0.115)
    # bax = brokenaxes(ylims=((0, 10), (120, 130)), despine=False)

    # Setup Grids (0 is top, 1 is bottom part)
    bax.axs[0].grid(which='major', axis='y', linestyle=':')
    bax.axs[1].grid(which='major', axis='y', linestyle=':')

    def ms_to_sec(num):
        return num / 1_000

    colors = ['#15607a', '#ffbd70', '#e7e7e7', '#ff483a']

    # Plot Bars
    width = 0.002
    for i, label in enumerate(labels):
        df = pandas_dataframes[i]
        if len(df) == 0:
            continue
        if not labels[i] in positions:
            continue
        else:

            # TODO: Remove this dirty
            if 'EVA' in label:
                width = 2 * width + 0.1
                for k in df.columns:
                    df[k] = 0
            else:
                width = 0.30

            x_pos = positions[label]
            d1 = ms_to_sec(df['t_keygen'].mean())
            d1_err = 0 if math.isnan(df['t_keygen'].std()) else df['t_keygen'].std()
            p1 = bax.bar(x_pos, d1, width, color=colors[0])

            d2 = ms_to_sec(df['t_input_encryption'].mean())
            d2_err = 0 if math.isnan(df['t_input_encryption'].std()) else df['t_input_encryption'].std()
            p2 = bax.bar(x_pos, d2, width, bottom=d1, color=colors[1])

            d3 = ms_to_sec(df['t_computation'].mean())
            d3_err = 0 if math.isnan(df['t_computation'].std()) else df['t_computation'].std()
            p3 = bax.bar(x_pos, d3, width, bottom=d1 + d2, color=colors[2])

            d4 = ms_to_sec(df['t_decryption'].mean())
            d4_err = 0 if math.isnan(df['t_decryption'].std()) else df['t_decryption'].std()
            total_err = ms_to_sec(d1_err + d2_err + d3_err + d4_err)
            p4 = bax.bar(x_pos, d4, width, yerr=total_err, ecolor='black', capsize=3, bottom=d1 + d2 + d3,
                         color=colors[3])
            print(labels[i].replace('\n', ' '), ": \n", d1, '\t', d2, '\t', d3, '\t', d4, '\t( total: ',
                  d1 + d2 + d3 + d4,
                  ')')

    # Setup axes
    bax.axs[1].tick_params(axis='x', which='major', labelsize=9)
    bax.axs[0].tick_params(axis='y', which='major', labelsize=8)
    bax.axs[1].tick_params(axis='y', which='major', labelsize=8)
    bax.axs[0].get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
    bax.axs[1].get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
    bax.axs[0].set_xticks(np.arange(len(sorted_labels)))
    bax.axs[1].set_xticks(np.arange(len(sorted_labels)))
    bax.axs[1].set_xticklabels(sorted_labels)
    bax.set_ylabel('Time [s]', labelpad=26)

    # plt.title('Runtime for Neural Network Benchmark', fontsize=10)

    # Add Legend
    plt.legend((p4[0], p3[0], p2[0], p1[0]),
               ('Decryption', 'Computation', 'Encryption', 'Key Generation'),
               loc='upper left')

    # Restore current figure
    plt.figure(previous_figure.number)

    return fig


if __name__ == '__main__':
    print("Testing ploting with nn example")
    data = [pd.read_csv(
        's3://sok-repository-eval-benchmarks/20200830_125813__231137930/nGraph-HE-LeNet5/ngraph-he-lenet5-learned_nn.csv')]
    labels = ['nGraph-HE-LeNet5']
    fig = plot(labels, data)
    fig.savefig("nn_plot_test.pdf")
    fig.show()  # has some alignment issues, but pdf is fine