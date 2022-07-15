import pickle
import numpy as np
import scipy
import pandas as pd
from matplotlib.widgets import Slider
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy import optimize
from sys import platform

import tkinter as tk
from tkinter import simpledialog

# import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.widgets import TextBox
import matplotlib.widgets as mwidgets
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from numba import jit, njit, float32, int32

# Constants:
sec_to_day = 1/(3600*24)


# Main function:
def interactive_search():
    pass


if __name__ == '__main__':

   userinput = 

    # Pick a TIC
    #TIC = 349907707
    TIC = 349832567

    # Read prewhitened frequencies
    pw_file = 'pw/freqs_comb_candicates_analysed.csv'
    pw = pd.read_csv(pw_file).query('tic == @TIC').reset_index(drop=True)

    # Read L-S periodogram
    pg_file = 'pg/pg_tess{TIC}_corrected_stitched.csv'
    pg = pd.read_csv(pg_file.format(TIC=TIC))
    # Add columns
    if not 'period' in pg:
        pg['period'] = 1/pg.freq
    if not 'e_period' in pw:
        pw['e_period'] = pw['e_frequency']*pw['period']**2

    # Define fundamental variables
    df_full = pd.DataFrame({})
    df_full['period'] = pw.sort_values(by=['period']).period.values
    df_full['amp'] = pw.sort_values(by=['period']).amp.values
    df_full['e_period'] = pw.sort_values(by=['period']).e_period.values
    df_full['selection'] = 1

    df = df_full.copy()

    # Shared variables of interaction

    ivar = {}
    ivar['keystroke'] = None
    ivar['pressed_button'] = None
    keystroke_i1 = 'control' if platform == 'linux' else 'i'
    keystroke_i2 = 'shift+control' if platform == 'linux' else 'I'

    plotted_lines = {}

    selection = {}

    comb_params = {}
    # np.median(pw_cluster.period.values)
    comb_params['P0'] = df.query('amp == amp.max()').period.values.item()
    comb_params['dP0'] = np.median(np.diff(df.period.values))
    comb_params['Sigma'] = 0
    comb_params['nr'] = 5
    comb_params['nl'] = 5

    # Generate template comb
    _P, _dP = pattern_period(comb_params['P0'],
                             comb_params['dP0'],
                             comb_params['Sigma'],
                             nr=comb_params['nr'],
                             nl=comb_params['nl'])

    comb_params['P'] = _P
    comb_params['dP'] = _dP

    comb_params['data'] = df

    # Periods in the data (observations after prewhitening) that agree with the comb within dP/4
   GoodnesFit
    comb_params['e_P0'] = 0
    comb_params['e_dP0'] = 0
    comb_params['e_Sigma'] = 0
    comb_params['mismatch'] = -1
    comb_params['residuals'] = -1
    comb_params['opt'] = None

    # Create figure
    # fig =  plt.figure(figszie=(10,14))
    fig = plt.figure(figsize=(18, 16))

    # Create axis grid

    big_vgrid = fig.add_gridspec(
        5, 1, height_ratios=[1.0, 0.5, 0.5, 0.2, 0.5], hspace=0.0)

    main_hgrid = big_vgrid[0].subgridspec(
        1, 2, width_ratios=[3, 1.2], wspace=0.2)
    # main_hgrid = fig.add_gridspec(1, 2, width_ratios=[3,1.2], wspace=0.2)

    sub_vgrid1 = main_hgrid[0].subgridspec(
        3, 1, height_ratios=[0.1, 1.0, 1.0], hspace=0.0)
    sub_vgrid2 = main_hgrid[1].subgridspec(
        3, 1, height_ratios=[0.1, 1.0, 1.0], hspace=0.0)

    ax_echelle = fig.add_subplot(sub_vgrid2[1:])
    ax_p = fig.add_subplot(sub_vgrid1[0])
    ax_pg = fig.add_subplot(sub_vgrid1[1])
    ax_dp = fig.add_subplot(sub_vgrid1[2])

    # Link x-axis for the periodogram and period-spacing pattern
    ax_pg.get_shared_x_axes().join(ax_dp, ax_pg)
    ax_pg.get_shared_x_axes().join(ax_p, ax_pg)

    ax_buttons = fig.add_subplot(big_vgrid[1])

    # Results plots before the last row (density plots)
    result1_hgrid = big_vgrid[2].subgridspec(
        1, 3, width_ratios=[1, 1, 1], wspace=0.3)

    result1_hgrid_1 = result1_hgrid[0].subgridspec(
        2, 1, height_ratios=[0.1, 1], hspace=0.1)
    result1_hgrid_2 = result1_hgrid[1].subgridspec(
        2, 1, height_ratios=[0.1, 1], hspace=0.1)
    result1_hgrid_3 = result1_hgrid[2].subgridspec(
        2, 1, height_ratios=[0.1, 1], hspace=0.1)

    ax_res_P0dP0 = fig.add_subplot(result1_hgrid_1[1])
    ax_res_dP0Sigma = fig.add_subplot(result1_hgrid_2[1])
    ax_res_SigmaP0 = fig.add_subplot(result1_hgrid_3[1])
    ax_res_P0dP0_cbar = fig.add_subplot(result1_hgrid_1[0])
    ax_res_dP0Sigma_cbar = fig.add_subplot(result1_hgrid_2[0])
    ax_res_SigmaP0_cbar = fig.add_subplot(result1_hgrid_3[0])

    # Results plots on the last row
    result2_hgrid = big_vgrid[4].subgridspec(
        1, 3, width_ratios=[1, 1, 1], wspace=0.3)

    ax_res_P0 = fig.add_subplot(result2_hgrid[0])
    ax_res_dP0 = fig.add_subplot(result2_hgrid[1])
    ax_res_Sigma = fig.add_subplot(result2_hgrid[2])

    ax_buttons.axis('off')

    # Plot available periods for the fit as triangles in the top panel
    plotted_lines['selection_p'] = []
    x = df.query('selection==1').period.values
    y = np.repeat(0.2, x.size)
    line, = ax_p.plot(x, y, color='k', marker=7, alpha=1,
                      zorder=2, picker=5, ls='None')
    l2 = line
    plotted_lines['selection_p'].append(line)

    # Plotted range
    ax_p.set_ylim(0, 1)
    ax_p.axis('off')

    # Plot the periodogram of the light curve
    x = pg.sort_values(by=['period']).period
    y = pg.sort_values(by=['period']).amp
    ax_pg.plot(x, y, lw=1, color='k', zorder=3)
    ax_pg.get_xaxis().set_visible(False)

    # Labels
    ax_pg.set_ylabel('amplitude')

    # Plotted range
    xlim = (df.period.min()-1/365, df.period.max()+1/365)
    ax_pg.set_xlim(xlim)

    # Overplot prewhitening results as vertical lines on the periodogram
    for p in df.period:
        ax_pg.axvline(p, color='k', ls='dotted', lw=1)

    # Plot prewhitening dP vs P
    plotted_lines['obs_dp'] = []
    _ = df.query('selection==1').period.values
    x = period_for_dP_plot(_, mode='middle')
    y = np.diff(_)
    line, = ax_dp.plot(x, y, lw=1, color='k', ls='dashed',
                       marker='.', zorder=2, picker=5)
    plotted_lines['obs_dp'].append(line)

    # Mark zero
    ax_dp.axhline(0, ls='dotted', lw=0.5, color='gray')

    # Labels
    ax_dp.set_xlabel('period (days)')
    ax_dp.set_ylabel('$\Delta P$ (days)')

    # Plot period echelle diagram
    colors = [choose_color(val) for val in df.selection.values]
    dp = np.median(np.diff(df.period.values))
    scatter1 = ax_echelle.scatter((df.period) % dp - dp, df.period,
                                  s=100.*(df.amp/df.amp.max()), color=colors, zorder=3, picker=5)
    scatter2 = ax_echelle.scatter((df.period) % dp + dp, df.period,
                                  s=100.*(df.amp/df.amp.max()), color=colors, zorder=3, picker=5)
    scatter3 = ax_echelle.scatter((df.period) % dp,      df.period,
                                  s=100.*(df.amp/df.amp.max()), color=colors, zorder=3, picker=5)
    plotted_lines['mod_dp'] = []
    line = ax_dp.axhline(dp, color='dodgerblue', lw=1, zorder=0, ls='dotted')
    plotted_lines['mod_dp'].append(line)

    # Plot period echelle diagram of the comb
    scatter11 = ax_echelle.scatter(
        (comb_params['P']) % dp - dp, comb_params['P'], s=30, color='r', zorder=4, alpha=0.3, marker='*')
    scatter22 = ax_echelle.scatter(
        (comb_params['P']) % dp + dp, comb_params['P'], s=30, color='r', zorder=4, alpha=0.3, marker='*')
    scatter33 = ax_echelle.scatter(
        (comb_params['P']) % dp,      comb_params['P'], s=30, color='r', zorder=4, alpha=0.3, marker='*')
    scatter44 = ax_echelle.scatter(
        (comb_params['P0']) % dp,      comb_params['P0'], s=30, color='gold', zorder=4, alpha=0.3, marker='*')

    # Plotted range
    ax_echelle.set_xlim(-dp, 2.*dp)
    ax_echelle.set_ylim(df.period.min(), df.period.max())

    # Separating the different echelles in the figure
    plotted_lines['echelle_vline'] = []
    ax_echelle.axvline(0,  ls='dashed', color='gray', lw=2, zorder=2)
    line = ax_echelle.axvline(dp, ls='dashed', color='gray', lw=2, zorder=2)
    plotted_lines['echelle_vline'].append(line)

    # Labels
    ax_echelle.set_ylabel('period (days)')
    ax_echelle.set_xlabel(f'period mod {dp:.5f} (days)')

    # Plot template comb on the periodogram
    plotted_lines['comb_pg'] = []
    for p in comb_params['P']:
        if p != comb_params['P0']:
            line = ax_pg.axvline(p, color='r', alpha=0.3, lw=2, zorder=0)
            plotted_lines['comb_pg'].append(line)
    line = ax_pg.axvline(
        comb_params['P0'], color='gold', alpha=0.9, lw=2, zorder=0)
    plotted_lines['comb_pg'].append(line)
    # plot_updated_comb(interactive=False)

    # Plot template dP vs P
    plotted_lines['comb_dp'] = []
    x = period_for_dP_plot(comb_params['P'], mode='middle')
    y = np.diff(comb_params['P'])
    line, = ax_dp.plot(x, y, lw=1, color='r', marker='*',
                       ls='solid', zorder=1, alpha=0.5)
    plotted_lines['comb_dp'].append(line)
    if comb_params['nr'] > 1:
        ind = np.abs(comb_params['P']-comb_params['P0']).argmin()
        _ = comb_params['P'][ind:ind+2]
        x = period_for_dP_plot(_, mode='middle')
        y = np.diff(_)
        line, = ax_dp.plot(x, y, lw=1, color='gold',
                           marker='*', ls='None', zorder=1, alpha=.5)
        plotted_lines['comb_dp'].append(line)

    # Plot observations that match the comb
    plotted_lines['observed_p'] = []

    # Sliders

    # Make space to place sliders
    fig.subplots_adjust(bottom=0.2, top=0.8)

    # Create axes for sliders
    height = 0.01
    width = 0.775
    x0 = 0.125
    y0 = 0.81
    hspacing = 0.015

    n_sliders = 5

    coords = [(x0, y0+i*hspacing, width, height) for i in range(n_sliders)]
    slider_axes = [fig.add_axes(
        c, facecolor='lightgoldenrodyellow') for c in coords]
    for ax in slider_axes:
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)

    ax_amp = slider_axes[0]
    ax_mod = slider_axes[1]
    ax_P0 = slider_axes[2]
    ax_dP0 = slider_axes[3]
    ax_Sigma = slider_axes[4]

    # Set range for each slider
    lower_limit_amp = 0
    upper_limit_amp = 1

    lower_limit_P0 = df.period.min()
    upper_limit_P0 = df.period.max()

    lower_limit_dP0 = 100*sec_to_day
    upper_limit_dP0 = 4000*sec_to_day

    _ = df.query('selection==1').period.values
    lower_limit_mod = min(lower_limit_dP0, np.diff(_).min())
    upper_limit_mod = max(upper_limit_dP0, np.diff(_).max())

    lower_limit_Sigma = -0.3  # -0.2 #Sigma + 0.5*Sigma
    upper_limit_Sigma = 0.3  # 0.2 #Sigma - 0.5*Sigma

    # Set resolution or step for each slider
    amp_resolution = 0.01
    mod_resolution = 0.0001
    P0_resolution = 0.0001
    dP0_resolution = 0.0001
    Sigma_resolution = 0.001

    # Set initial value for each slider
    amp_initial = 0
    mod_initial = comb_params['dP0']
    P0_initial = comb_params['P0']
    dP0_initial = comb_params['dP0']
    Sigma_initial = comb_params['Sigma']

    slider_amp = Slider(ax_amp,   'ampl',         lower_limit_amp,   upper_limit_amp,
                        valinit=amp_initial,   valfmt='%1.2f', facecolor='k', valstep=amp_resolution)
    slider_mod = Slider(ax_mod,   'mod',          lower_limit_mod,   upper_limit_mod,
                        valinit=mod_initial,   valfmt='%1.6f d', facecolor='dodgerblue', valstep=mod_resolution)
    slider_P0 = Slider(ax_P0,    '$P_0$',        lower_limit_P0,    upper_limit_P0,
                       valinit=P0_initial,    valfmt='%1.6f d', facecolor='gold',  valstep=P0_resolution)
    slider_dP0 = Slider(ax_dP0,   '$\Delta P_0$', lower_limit_dP0,   upper_limit_dP0,
                        valinit=dP0_initial,   valfmt='%1.6f d', facecolor='red',  valstep=dP0_resolution)
    slider_Sigma = Slider(ax_Sigma, '$\Sigma$',     lower_limit_Sigma, upper_limit_Sigma,
                          valinit=Sigma_initial, valfmt='%1.6f',   facecolor='red',  valstep=Sigma_resolution)

    slider_amp.on_changed(update_amplitude_tolerance)
    slider_mod.on_changed(update_echelle)
    slider_P0.on_changed(update_comb)
    slider_dP0.on_changed(update_comb)
    slider_Sigma.on_changed(update_comb)

    # Build interface: Buttons, checkbox and textbox

    height = 0.02
    width = 0.1
    dwidth = 0.01
    button_x0 = 0.10
    # button_y0 = 0.07

    xpos = button_x0+0*(width+dwidth)
    # ypos = 0.07
    # ypos = 0.4
    ypos = ax_buttons.get_position().y0+ax_buttons.get_position().height/2

    # <<--------------|
    button_fit = Button(plt.axes([xpos, ypos, width, height]), 'Fit')
    button_fit.on_clicked(do_fit3)

    xpos = button_x0+0*(width+dwidth)
    button_explore = Button(plt.axes(
        [xpos, ypos-height, width, height]), 'Explore result')  # <<--------------|
    button_explore.on_clicked(explore_results2)

    xpos = button_x0+1*(width+dwidth)
    button_save = Button(plt.axes([xpos, ypos, width, height]), 'Save')
    button_save.on_clicked(save)

    # Boxes

    xpos = button_x0+2*(width+dwidth)
    ax_box_mismatch = plt.axes([xpos, ypos, width, height])
    text_box_mismatch = TextBox(
        ax_box_mismatch, "", initial='0', color='lightgoldenrodyellow')
    ax_box_mismatch.set_title('Mismatch', y=-1.0)

    xpos = button_x0+3*(width+dwidth)
    ax_box_residuals = plt.axes([xpos, ypos, width, height])
    text_box_residuals = TextBox(
        ax_box_residuals, "", initial='0', color='lightgoldenrodyellow')
    ax_box_residuals.set_title('Residuals (d)', y=-1.0)

    xpos = button_x0+4*(width+dwidth)
    ax_box_P0 = fig.add_axes([xpos, ypos, width, height])
    text_box_P0 = TextBox(ax_box_P0, "", initial='0')
    cid_box_P0 = text_box_P0.on_submit(read_box_P0)
    ax_box_P0.set_title('$P_0$ (d)', y=-1.0)

    xpos = button_x0+5*(width+dwidth)
    ax_box_dP0 = fig.add_axes([xpos, ypos, width, height])
    text_box_dP0 = TextBox(ax_box_dP0, "", initial='0')
    cid_box_dP0 = text_box_dP0.on_submit(read_box_dP0)
    ax_box_dP0.set_title('$\Delta P_0$ (d)', y=-1.0)

    xpos = button_x0+6*(width+dwidth)
    ax_box_Sigma = fig.add_axes([xpos, ypos, width, height])
    text_box_Sigma = TextBox(ax_box_Sigma, "", initial='0')
    cid_box_Sigma = text_box_Sigma.on_submit(read_box_Sigma)
    ax_box_Sigma.set_title('$\Sigma$', y=-1.0)

    # Connect ID to the plot visualization
    cid_key = fig.canvas.mpl_connect('key_press_event', read_keystroke)
    cid_button = fig.canvas.mpl_connect('button_press_event', read_button)
    cid_pick = fig.canvas.mpl_connect('pick_event', update_selection)

    # Properties of the rectangle-span area-selector
    rect_props = dict(facecolor='grey', alpha=0.20)
    # Area selector
    span = mwidgets.SpanSelector(
        ax_pg, onselect, 'horizontal', rectprops=rect_props, useblit=False)

    ax_Sigma.set_title(f'TIC {TIC}')
    # fig.suptitle(f'TIC {TIC}')

    # mng = plt.get_current_fig_manager()
    # mng.resize(1800,800)

    # figManager = plt.get_current_fig_manager()
    # figManager.window.showMaximized()
    # ScrollableWindow(fig)
    plt.show()

    # Disconnect from the plot visuzlization
    for cid in [cid_key, cid_button, cid_pick, cid_box_P0, cid_box_dP0, cid_box_Sigma]:
        fig.canvas.mpl_disconnect(cid)
