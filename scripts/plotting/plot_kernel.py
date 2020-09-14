import math
from typing import List
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter
import itertools
import operator
from operator import add


def human_format(num):
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])


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
    fig_width = 252 * inches_per_pt  # width in inches
    # fig_height = (fig_width * golden_mean)  # height in inches
    fig_height = 2.5
    figsize = [fig_width * 0.67, fig_height / 1.22]

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

    positions = {
        'SEAL-BFV': (0, 0),
        'E3-SEAL': (0, 1),
        'TFHE': (1, 0),
        'E3-TFHE': (1, 1),
        'Cingulata': (2, 0),
    }

    # plt.title('Runtime for Chi-Squared Test Benchmark', fontsize=10)
    plt.ylabel('Time [s]', labelpad=0)

    bar_width = 0.002
    spacer = 0.01
# {\fontsize{30pt}{3em}\selectfont{}{Mean WRFv3.5 LHF\r}{\fontsize{18pt}{3em}\selectfont{}(September 16 - October 30, 2012)}
    group_labels = [
        'SEAL\n{\\fontsize{7pt}{3em}\\selectfont{}(Native/E\\textsuperscript{3})}',
        'TFHE\n{\\fontsize{7pt}{3em}\\selectfont{}(Native/E\\textsuperscript{3})}',
        'Cingulata'
    ]

    # ['E3-SEAL', 'E3-TFHE', 'SEAL-BFV-Batched', 'SEAL-BFV', 'TFHE']

    def get_x_ticks_positions():
        group_widths = []
        x_pos_start = []
        x_pos_end = []
        for key, group in itertools.groupby(positions.values(), operator.itemgetter(0)):
            group_widths.append(len(list(group)) * bar_width)
        for w in group_widths:
            if not x_pos_start:
                x_pos_start.append(0 + spacer)
            else:
                x_pos_start.append(x_pos_end[-1] + spacer)
            x_pos_end.append(x_pos_start[-1] + w)
        result = list(map(add, x_pos_start, [w / 2 for w in group_widths]))
        return result, x_pos_start

    x_center, x_start = get_x_ticks_positions()

    def get_x_position(group_pos: tuple) -> int:
        return x_start[group_pos[0]] + (group_pos[1] * bar_width) + (bar_width / 2)

    plt.yscale('log')

    plt.xticks(x_center, group_labels, fontsize=9)  # rotation='35',)
    # adds a thousand separator
    # fig.axes[0].get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
    fig.axes[0].get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: human_format(x)))
    # add a grid
    ax = plt.gca()
    ax.grid(which='major', axis='y', linestyle=':')

    def ms_to_sec(num):
        return num / 1_000

    colors = ['0.1', '0.35', '0.5', '0.85']
    hatches = ['', '.', '///', '']

    # Plot Bars
    max_y_value = 0
    for i in range(len(labels)):
        if not labels[i] in positions:
            continue
        else:
            x_pos = get_x_position(positions[labels[i]])
        df = pandas_dataframes[i]
        d1 = ms_to_sec(df['t_keygen'].mean())
        d1_err = 0 if math.isnan(df['t_keygen'].std()) else df['t_keygen'].std()
        p1 = plt.bar(x_pos, d1, bar_width * 0.9, color=colors[0], hatch=hatches[0])
        d2 = ms_to_sec(df['t_input_encryption'].mean())
        d2_err = 0 if math.isnan(df['t_input_encryption'].std()) else df['t_input_encryption'].std()
        p2 = plt.bar(x_pos, d2, bar_width * 0.9, bottom=d1, color=colors[1], hatch=hatches[1])
        d3 = ms_to_sec(df['t_computation'].mean())
        d3_err = 0 if math.isnan(df['t_computation'].std()) else df['t_computation'].std()
        p3 = plt.bar(x_pos, d3, bar_width * 0.9, bottom=d1 + d2, color=colors[2], hatch=hatches[2])
        d4 = ms_to_sec(df['t_decryption'].mean())
        d4_err = 0 if math.isnan(df['t_decryption'].std()) else df['t_decryption'].std()
        total_err = ms_to_sec(d1_err + d2_err + d3_err + d4_err)
        max_y_value = d1 + d2 + d3 + d4 if (d1 + d2 + d3 + d4) > max_y_value else max_y_value
        p4 = plt.bar(x_pos, d4, bar_width * 0.9, yerr=total_err, ecolor='black', capsize=3, bottom=d1 + d2 + d3,
                     color=colors[3], hatch=hatches[3])
        print(labels[i].replace('\n', ' '), ": \n", d1, '\t', d2, '\t', d3, '\t', d4, '\t( total: ', d1 + d2 + d3 + d4,
              ')')

    max_y_rounded = (int(math.ceil(max_y_value / 10.0)) * 10) + 10
    # plt.yticks(np.arange(0, max_y_rounded, step=10))

    plt.yticks(fontsize=8)

    # Add Legend
    plt.legend((p4[0], p3[0], p2[0], p1[0]), ('Decryption', 'Computation', 'Encryption', 'Key Generation'))

    # Restore current figure
    plt.figure(previous_figure.number)

    return fig