import redshift_to_time     # redshift_to_time.cosmology_calculator(redshift, H0=69.6, WM=0.286, WV=0.714)
import data_Muzzin13        # data_Muzzin13.function_data(redshift, log(galaxy_stellar_mass))
import time
import random
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib.gridspec as gridspec


def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings
    Code copy from https://stackoverflow.com/a/50029441

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """
    ax = plt.gca() if ax is None else ax
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)
    lc.set_array(np.asarray(c))
    ax.add_collection(lc)
    ax.autoscale()
    return lc


def minimum_sft(galaxy_mass__):
    '''Thomas05 equaiton 5'''
    sft = math.exp(3.67-0.37*galaxy_mass__)
    return sft


def minimum_sft_new(galaxy_mass__):
    '''Yan 2019b'''
    sft = math.exp(8-0.7*galaxy_mass__)
    return sft


def generate_a_galaxy(age_of_the_universe, generated_galaxy_fraction):
    '''
    :param age_of_the_universe:
    :return: [star_formation_start_time, star_formation_stop_time, log_sfr]
    '''
    star_formation_start_time = random.random() * time_list[0]  # * generated_galaxy_fraction**0.25  # in [Gyr]
    # 10^-2.3 is the lowest  possible SFR
    # 10^5    is the highest possible SFR
    log_sfr = random.random() * (2.3 + 5) - 2.3  # log_10(star formation rate [solar mass per year])
    maximum_formation_time = age_of_the_universe - star_formation_start_time
    star_formation_stop_time = star_formation_start_time + random.random() * maximum_formation_time  # in [Gyr]
    if star_formation_stop_time > time_list[0]:
        star_formation_stop_time = time_list[0]*1.001
    log_galaxy_final_mass = 9 + math.log((star_formation_stop_time - star_formation_start_time), 10) + log_sfr
    # Here does not consider that
    # the final dynamical mass should be smaller than the total stellar mass formed
    star_formation_timescale = star_formation_stop_time - star_formation_start_time
    # minimum_formation_time = minimum_sft(log_galaxy_final_mass)
    minimum_formation_time = 0.3
    if 8 < log_galaxy_final_mass < 13 and star_formation_timescale > minimum_formation_time:
        output = [star_formation_start_time, star_formation_stop_time, log_sfr,
                  log_galaxy_final_mass, star_formation_timescale]
        return output
    else:
        return generate_a_galaxy(age_of_the_universe, generated_galaxy_fraction)


def data_prepare():
    '''
    There might be a different way which is to use the lower redshift ranges instead of middle redshifts.
    '''
    print('===\ndata_prepare...')
    global redshift_list, time_list, logMstar_list
    redshift_list_ = data_Muzzin13.redshift_list
    logMstar_list = data_Muzzin13.logMstar_list
    iii = 0
    redshift_list = []
    while iii < len(redshift_list_) - 1:
        redshift_ = (redshift_list_[iii] + redshift_list_[iii + 1]) / 2
        redshift_list.append(redshift_)
        (iii) = (iii + 1)
    time_list = []
    for redshift__ in redshift_list:
        time_ = redshift_to_time.cosmology_calculator(redshift__, H0, WM, WV)
        time_list.append(time_)
    print("redshifts:", redshift_list)
    print("time steps:", time_list)
    print("logMstar limits:", logMstar_list)
    print('===\n')
    return

def compute_galaxy_mass_at_each_time(a_galaxy):
    '''
    galaxy_sfh = [info_time_1, info_time_2...]
    info_time_i = [None_or_log(galaxy-mass), star-forming(1)_or_quiescent(0)]
    galaxy_sfh = compute_galaxy_mass_at_each_time(a_galaxy)
    '''
    galaxy_sfh = []
    for time__ in time_list:
        if a_galaxy[0] < time__:
            if a_galaxy[1] < time__:
                log_galaxy_final_mass = a_galaxy[3]
                None_or_log_galaxy_mass = log_galaxy_final_mass
                star_forming_or_quiescent = 0
            else:
                log_galaxy_mass_at_this_time = 9 + math.log((time__ - a_galaxy[0]), 10) + a_galaxy[2]
                None_or_log_galaxy_mass = log_galaxy_mass_at_this_time
                star_forming_or_quiescent = 1
        else:
            None_or_log_galaxy_mass = None
            star_forming_or_quiescent = 0
        galaxy_sfh.append([None_or_log_galaxy_mass, star_forming_or_quiescent])
    return galaxy_sfh


def compare_with_observation(galaxy_number_a__, galaxy_number_q__, galaxy_number_f__,
                                     galaxy_mass_a__, galaxy_mass_q__, galaxy_mass_f__,
                                     total_galaxy_number_at_low_z__, total_galaxy_mass_at_low_z__):
    '''
    Here calcualte the difference between model galaxies and observations
    :param galaxy_number_a__:
    :param galaxy_number_q__:
    :param galaxy_number_f__:
    :param galaxy_mass_a__:
    :param galaxy_mass_q__:
    :param galaxy_mass_f__:
    :param total_galaxy_number_at_low_z__:
    :param total_galaxy_mass_at_low_z__:
    :return:
    '''
    error__ = 0
    error_power = 2
    for ii in range(time_step_number):
        for jj in range(mass_range_number):
            if not ii > jj + 2:  # Note: At high redshift, use only the observations of massive galaxies
                # # if jj > 2 or ii > 4:
                # #     error_power = 1.5
                # # else:
                # #     error_power = 1
                # # error_power = redshift_list[ii] * logMstar_list[jj]**4
                # # error_power = redshift_list[ii] ** 1
                # # error__ += abs(galaxy_number_a__[ii][jj]/total_galaxy_number_at_low_z__ - \
                # #            data_Muzzin13.function_data(ii, jj)[0]/total_galaxy_number_at_low_z_obs)*error_power\
                # #            /data_Muzzin13.function_data(ii, jj)[0]
                # error__ += abs(galaxy_number_q__[ii][jj] / total_galaxy_number_at_low_z__ - \
                #                data_Muzzin13.function_data(ii, jj)[1] / total_galaxy_number_at_low_z_obs) * error_power \
                #            / data_Muzzin13.function_data(ii, jj)[1]
                # error__ += abs(galaxy_number_f__[ii][jj] / total_galaxy_number_at_low_z__ - \
                #                data_Muzzin13.function_data(ii, jj)[2] / total_galaxy_number_at_low_z_obs) * error_power \
                #            / data_Muzzin13.function_data(ii, jj)[2]
                # # error__ += abs(galaxy_mass_a__[ii][jj]/total_galaxy_mass_at_low_z__ - \
                # #            data_Muzzin13.function_data(ii, jj)[3]/total_galaxy_mass_at_low_z_obs)*error_power\
                # #            /data_Muzzin13.function_data(ii, jj)[0]
                # error__ += abs(galaxy_mass_q__[ii][jj] / total_galaxy_mass_at_low_z__ - \
                #                data_Muzzin13.function_data(ii, jj)[4] / total_galaxy_mass_at_low_z_obs) * error_power \
                #            / data_Muzzin13.function_data(ii, jj)[1]
                # error__ += abs(galaxy_mass_f__[ii][jj] / total_galaxy_mass_at_low_z__ - \
                #                data_Muzzin13.function_data(ii, jj)[5] / total_galaxy_mass_at_low_z_obs) * error_power \
                #            / data_Muzzin13.function_data(ii, jj)[2]

                model_number = galaxy_number_q__[ii][jj]/total_galaxy_number_at_low_z__
                if model_number == 0:
                    model_number = 1E-7
                error__ += abs(math.log(model_number, 10) - math.log(data_Muzzin13.function_data(ii, jj)[1]
                           /total_galaxy_number_at_low_z_obs, 10)) ** error_power
                model_number = galaxy_number_f__[ii][jj]/total_galaxy_number_at_low_z__
                if model_number == 0:
                    model_number = 1E-6
                error__ += abs(math.log(model_number, 10) - math.log(data_Muzzin13.function_data(ii, jj)[2]
                           /total_galaxy_number_at_low_z_obs, 10)) ** error_power
                model_mass = galaxy_mass_q__[ii][jj]/total_galaxy_mass_at_low_z__
                if model_mass == 0:
                    model_mass = 10**4.62
                error__ += abs(math.log(model_mass) - math.log(data_Muzzin13.function_data(ii, jj)[4]
                           /total_galaxy_mass_at_low_z_obs)) ** error_power
                model_mass = galaxy_mass_f__[ii][jj]/total_galaxy_mass_at_low_z__
                if model_mass == 0:
                    model_mass = 10**5.17
                error__ += abs(math.log(model_mass) - math.log(data_Muzzin13.function_data(ii, jj)[5]
                           /total_galaxy_mass_at_low_z_obs)) ** error_power
    return error__


def plot_compare_with_observation(galaxy_number_a__, galaxy_number_q__, galaxy_number_f__,
                                  total_galaxy_number_at_low_z__):
    lines = [[[], [], [], [], [], [], []], [[], [], [], [], [], [], []], [[], [], [], [], [], [], []]]
    lines_obs = [[[], [], [], [], [], [], []], [[], [], [], [], [], [], []], [[], [], [], [], [], [], []]]

    for ii in range(time_step_number):
        for jj in range(mass_range_number):
            if galaxy_number_a__[ii][jj] == 0:
                lines[0][ii].append(-9)
            else:
                lines[0][ii].append(math.log(galaxy_number_a__[ii][jj]/total_galaxy_number_at_low_z__*total_galaxy_number_at_low_z_obs, 10))
            lines_obs[0][ii].append(math.log(data_Muzzin13.function_data(ii, jj)[0], 10))
            if galaxy_number_q__[ii][jj] == 0:
                lines[1][ii].append(-9)
            else:
                lines[1][ii].append(math.log(galaxy_number_q__[ii][jj]/total_galaxy_number_at_low_z__*total_galaxy_number_at_low_z_obs, 10))
            lines_obs[1][ii].append(math.log(data_Muzzin13.function_data(ii, jj)[1], 10))
            if galaxy_number_f__[ii][jj] == 0:
                lines[2][ii].append(-9)
            else:
                lines[2][ii].append(math.log(galaxy_number_f__[ii][jj]/total_galaxy_number_at_low_z__*total_galaxy_number_at_low_z_obs, 10))
            lines_obs[2][ii].append(math.log(data_Muzzin13.function_data(ii, jj)[2], 10))
    return lines, lines_obs


def plot_all_galaxy_result():
    # print("galaxies:", galaxies)
    # print("galaxies_sfh_for_all_galaxy:", galaxies_sfh_for_all_galaxy)
    # print("galaxy_info_at_each_time:", galaxy_info_at_each_time)
    # print("galaxy_number_a:", galaxy_number_a)
    # print("galaxy_number_q:", galaxy_number_q)
    # print("galaxy_number_f:", galaxy_number_f)
    # print("galaxy_mass_a:", galaxy_mass_a)
    # print("galaxy_mass_q:", galaxy_mass_q)
    # print("galaxy_mass_f:", galaxy_mass_f)
    print("total_galaxy_number_at_low_z:", total_galaxy_number_at_low_z)
    # print("total_galaxy_mass_at_low_z:", total_galaxy_mass_at_low_z)
    print("log_average_galaxy_mass_at_low_z:",
          math.log(total_galaxy_mass_at_low_z/total_galaxy_number_at_low_z, 10))
    print('plotting results....................')

    ########## plot downsizing relation ##########

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(0, figsize=(6, 5))
    gap = int(len(galaxies)/the_number_of_galaxies_being_plotted)
    for index__ in range(gap, len(galaxies)+1, gap):
        a_galaxy = galaxies[index__]
        plot_galaxy_mass = a_galaxy[3]
        plot_galaxy_sft = a_galaxy[4]
        the_mean_age = age_of_the_universe - (a_galaxy[0] + a_galaxy[1]) / 2
        normalize_the_mean_age = (the_mean_age - age_of_the_universe)/time_list[0] + 1
        # normalize_the_mean_age = a_galaxy[1]/age_of_the_universe
        color = plt.cm.rainbow(normalize_the_mean_age)
        normalize_the_mass = (plot_galaxy_mass - 8) / 5
        normalize_the_mass_2 = normalize_the_mass ** 1.3 / 3 + 0.33
        plt.scatter(plot_galaxy_mass, plot_galaxy_sft, s=22, alpha=normalize_the_mass_2, c=color)
    age_ticks = [2, 4, 6, 8, 10, 12, 14]
    age_ticks_normalize = []
    for ticks in age_ticks:
        ticks_normalized = (ticks - age_of_the_universe)/time_list[0] + 1
        age_ticks_normalize.append(ticks_normalized)
    lc = multiline([8, 8], [0, 0], [0, 1], cmap='rainbow')
    # axcb = fig.colorbar(lc, ticks=age_ticks_normalize)
    axcb = fig.colorbar(lc)
    axcb.set_label('mean stellar age [Gyr]')
    axcb.ax.set_yticklabels(['2', '4', '6', '8', '10', '12', '14'])
    plt.xlabel(r'log$_{10}$(M$_{dyn}$ [M$_\odot$])')
    plt.ylabel(r'$t_{\rm sf}$ [Gyr]')
    plt.tight_layout()
    if save_figure == True:
        plt.savefig('downsizing.pdf', dpi=250)

    ########## plot SFH ##########

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(1, figsize=(6, 5))
    plt.xlabel(r'age of the Universe [Gyr]')
    plt.ylabel(r'log$_{10}$($SFR$ [M$_\odot$/yr])')
    for index__ in range(gap, len(galaxies)+1, gap):
        a_galaxy = galaxies[index__]
        star_formation_start_time = a_galaxy[0]
        star_formation_stop_time = a_galaxy[1]
        log_sfr = a_galaxy[2]
        galaxy_mass__ = a_galaxy[3]
        normalize_the_mass = (galaxy_mass__ - 8) / 5
        # color = plt.cm.hsv(normalize_the_mass)
        sft = a_galaxy[4]
        normalize_the_sft = sft / 14
        color = plt.cm.rainbow_r(normalize_the_sft)
        x = [star_formation_start_time, star_formation_stop_time, star_formation_stop_time, star_formation_start_time]
        y = [log_sfr, log_sfr, -2.3, -2.3]
        normalize_the_mass_2 = normalize_the_mass ** 2 / 5 + 0.05
        plt.fill(x, y, facecolor=color, alpha=normalize_the_mass_2, edgecolor='k')
    # mass_ticks = [8, 9, 10, 11, 12, 13]
    # mass_ticks_normalize = []
    # for ticks in mass_ticks:
    #     ticks_normalized = (ticks - 8) / 5
    #     mass_ticks_normalize.append(ticks_normalized)
    sft_ticks = [0, 2, 4, 6, 8, 10, 12, 14]
    sft_ticks_normalize = []
    for ticks in sft_ticks:
        ticks_normalized = ticks / 14
        sft_ticks_normalize.append(ticks_normalized)
    lc = multiline([0, 0], [0, 0], [0, 1], cmap='rainbow_r')
    # axcb = fig.colorbar(lc, ticks=mass_ticks_normalize)
    axcb = fig.colorbar(lc)
    # axcb.set_label(r'log$_{10}$($M_{dyn}$ [M$_\odot$])')
    # axcb.ax.set_yticklabels(['8', '9', '10', '11', '12', '13'])
    axcb.set_label(r't$_{sf}$ [Gyr]')
    axcb.ax.set_yticklabels(['0', '2', '4', '6', '8', '10', '12', '14'])
    plt.tight_layout()
    if save_figure == True:
        plt.savefig('SFH.pdf', dpi=250)

    ########## plot age--SFT relation ##########

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(2, figsize=(6, 5))
    plt.xlabel(r'mean stellar age [Gyr]')
    plt.ylabel(r'$t_{\rm sf}$ [Gyr]')

    plt.plot([0, age_of_the_universe/2], [0, age_of_the_universe], lw=0.8, ls='dashed', c='k',
             label='all galaxies end star formation at the same time')
    # the line with all galaxies still forming stars at the present day

    nolt1 = age_of_the_universe - time_list[0]*1.001  # nolt1: nearest obs lookback time
    plt.plot([nolt1, (age_of_the_universe+nolt1)/2], [0, age_of_the_universe-nolt1], lw=0.8, ls='dashed', c='k')
    # the line with all galaxies still forming stars at the nearest observation time

    plt.plot([age_of_the_universe/2, age_of_the_universe], [age_of_the_universe, 0], lw=0.8, c='k')
    # the line with all galaxies start forming star at the beginning of the time

    nolt = age_of_the_universe - time_list[0]  # nolt: nearest obs lookback time
    plt.plot([nolt/2, nolt], [nolt, 0], lw=0.8, c='k', label='all galaxies start forming star at the same time')
    # the line with all galaxies start forming star at the nearest observation time

    minimum_formation_time = minimum_sft(13)
    plt.plot([0, age_of_the_universe], [minimum_formation_time, minimum_formation_time], lw=0.8, ls='dotted',
             label=r'lowest possible $t_{sf}$')
    # the line with shortest possible sft

    for index__ in range(gap, len(galaxies)+1, gap):
        a_galaxy = galaxies[index__]
        plot_mean_age = age_of_the_universe-(a_galaxy[0]+a_galaxy[1])/2
        plot_galaxy_sft = a_galaxy[4]
        normalize_the_mass = (a_galaxy[3]-8)/5
        color = plt.cm.rainbow(normalize_the_mass)
        normalize_the_mass_2 = normalize_the_mass ** 2 / 2 + 0.1
        plt.scatter(plot_mean_age, plot_galaxy_sft, s=33, alpha=normalize_the_mass_2, c=color)

    mass_ticks = [8, 9, 10, 11, 12, 13, 14]
    mass_ticks_normalize = []
    for ticks in mass_ticks:
        ticks_normalized = (ticks-8) / 5
        mass_ticks_normalize.append(ticks_normalized)
    lc = multiline([8, 8], [0, 0], [0, 1], cmap='rainbow')
    # axcb = fig.colorbar(lc, ticks=mass_ticks_normalize)
    axcb = fig.colorbar(lc)
    axcb.set_label(r'log$_{10}$(M$_{dyn}$ [M$_\odot$])')
    axcb.ax.set_yticklabels(['8', '9', '10', '11', '12', '13', '14'])
    plt.legend(prop={'size': 7}, loc='best')
    plt.tight_layout()
    if save_figure == True:
        plt.savefig('triangle.pdf', dpi=250)

    ########## Madau plot ###########

    Madau_redshift = np.arange(0.5, 8, 0.2)
    Madau_time = []
    Madau_SFR_at_different_time = []
    for a_redshift in Madau_redshift:
        Madau_time.append(redshift_to_time.cosmology_calculator(a_redshift))
        Madau_SFR_at_different_time.append(0)
    time_index = len(Madau_SFR_at_different_time) - 1
    while time_index > 0:
        for a_galaxy in galaxies:
            if a_galaxy[0] < Madau_time[time_index] < a_galaxy[1]:
                Madau_SFR_at_different_time[time_index] += a_galaxy[2]
        (time_index) = (time_index - 1)
    for a_time in range(len(Madau_SFR_at_different_time)):
        if Madau_SFR_at_different_time[a_time] > 0:
            Madau_SFR_at_different_time[a_time] = math.log(Madau_SFR_at_different_time[a_time], 10)

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.figure(3, figsize=(6, 5))
    plt.xlabel(r'redshift')
    plt.ylabel(r'log$_{10}$($SFR$ [M$_\odot$/yr])')
    plt.loglog(Madau_redshift, Madau_SFR_at_different_time)
    plt.tight_layout()
    if save_figure == True:
        plt.savefig('Madau.pdf', dpi=250)

    ######## plot goodness of the fit ##########
    mass_list = [8.5, 9.5, 10.5, 11.25, 12]
    lines, lines_obs = plot_compare_with_observation(galaxy_number_a_new, galaxy_number_q_new, galaxy_number_f_new,
                                  total_galaxy_number_at_low_z_new)
    fig = plt.figure(4, figsize=(10, 5))

    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132, sharey=ax1)
    ax3 = plt.subplot(133, sharey=ax1)
    for time__ in range(time_step_number):
        normalize_the_time = time__/6.1+0.05
        color = plt.cm.rainbow(normalize_the_time)
        ax1.plot(mass_list, lines[0][time__], c=color, label='z={}'.format(redshift_list[time__]))
        ax1.plot(mass_list, lines_obs[0][time__], ls='dashed', c=color, alpha=0.11)
        ax2.plot(mass_list, lines[1][time__], c=color)
        ax2.plot(mass_list, lines_obs[1][time__], ls='dashed', c=color, alpha=0.11)
        ax3.plot(mass_list, lines[2][time__], c=color, label='Time={} Gyr'.format(round(time_list[time__], 1)))
        ax3.plot(mass_list, lines_obs[2][time__], ls='dashed', c=color, alpha=0.11)
        # select only galaxy obsvations that is used for fitting
        for mass__ in range(mass_range_number):
            if time__ > mass__ + 2:  # select only massive galaxy obsvations at high redshift
                lines_obs[0][time__][mass__] = None
                lines_obs[1][time__][mass__] = None
                lines_obs[2][time__][mass__] = None
        ax1.plot(mass_list, lines_obs[0][time__], ls='dashed', c=color)
        ax2.plot(mass_list, lines_obs[1][time__], ls='dashed', c=color)
        ax3.plot(mass_list, lines_obs[2][time__], ls='dashed', c=color)

    fig.subplots_adjust(wspace=0)  # make subplots close to each other
    plt.setp([a.get_yticklabels() for a in [fig.axes[1], fig.axes[2]]], visible=False)

    ax1.set_xlabel(r'log$_{10}$(M$_{dyn}$ [M$_\odot$])')
    ax2.set_xlabel(r'log$_{10}$(M$_{dyn}$ [M$_\odot$])')
    ax3.set_xlabel(r'log$_{10}$(M$_{dyn}$ [M$_\odot$])')
    ax1.set_ylabel(r'log$_{10}$($\Phi$ [# Mpc$^{-3}$])')

    # ax1.title.set_text('All')
    ax1.set_title('All')
    ax2.set_title(r'Number of galaxies with mass smaller than 10$^{13}$ M$_\odot$ and larger than the given value'
                       '\nQuiescent')
    ax3.set_title('Star Forming')
    ax1.plot([], [], c='k', label='Model')
    ax1.plot([], [], c='k', label='Observation', ls='dashed')
    ax1.legend(prop={'size': 7}, loc='lower left')
    ax3.legend(prop={'size': 7}, loc='lower left')
    if save_figure == True:
        plt.savefig('goodness_of_fit.pdf', dpi=250)

    plt.show()
    return


def output_text_file():
    file = open('generated_galaxies.txt', 'w')
    file.write("# galaxies:\n")
    file.write("{}\n".format(galaxies))
    file.write("# galaxies_sfh_for_all_galaxy:\n")
    file.write("{}\n".format(galaxies_sfh_for_all_galaxy))
    file.write("# galaxy_info_at_each_time:\n")
    file.write("{}\n".format(galaxy_info_at_each_time))
    file.write("# galaxy_number_a:\n")
    file.write("{}\n".format(galaxy_number_a))
    file.write("# galaxy_number_q:\n")
    file.write("{}\n".format(galaxy_number_q))
    file.write("# galaxy_number_f:\n")
    file.write("{}\n".format(galaxy_number_f))
    file.write("# galaxy_mass_a:\n")
    file.write("{}\n".format(galaxy_mass_a))
    file.write("# galaxy_mass_q:\n")
    file.write("{}\n".format(galaxy_mass_q))
    file.write("# galaxy_mass_f:\n")
    file.write("{}\n".format(galaxy_mass_f))
    file.write("# total_galaxy_number_at_low_z:\n")
    file.write("{}\n".format(total_galaxy_number_at_low_z))
    file.write("# total_galaxy_mass_at_low_z:\n")
    file.write("{}\n".format(total_galaxy_mass_at_low_z))
    file.write("# error:\n")
    file.write("{}\n".format(error))
    file.close()
    return


if __name__ == '__main__':
    start = time.time()
    H0 = 69.6
    WM = 0.286
    WV = 0.714
    # Change cosmology (and/or IMF assumption) here is not self-consistent
    # since the observational galaxy mass already assumed the standard model.
    # For a universe with longer lifetime, the real galaxy mass should be higher than the adopted observational value.
    # Thus this simulation result is biased towards lower galaxy mass.
    # Nevertheless, the influence of modifying the cosmology here seems not significant.
    age_of_the_universe = redshift_to_time.cosmology_calculator(0, H0, WM, WV)
    print("Modeled age of the universe:", age_of_the_universe)
    data_prepare()
    time_step_number = len(time_list)
    mass_range_number = len(logMstar_list)-1

    total_galaxy_number_at_low_z_obs = 0
    total_galaxy_mass_at_low_z_obs = 0
    for mass_range_ in range(mass_range_number):
        total_galaxy_number_at_low_z_obs += data_Muzzin13.function_data(0, mass_range_)[0]
        total_galaxy_mass_at_low_z_obs += data_Muzzin13.function_data(0, mass_range_)[3]
    # print('total_galaxy_number_at_low_z_obs', total_galaxy_number_at_low_z_obs)
    # print('total_galaxy_mass_at_low_z_obs', total_galaxy_mass_at_low_z_obs)
    print('log_average_galaxy_mass_at_low_z_obs',
          math.log(total_galaxy_mass_at_low_z_obs/total_galaxy_number_at_low_z_obs, 10))

    galaxies = []
    galaxies_sfh_for_all_galaxy = []
    galaxy_info_at_each_time = []
    galaxy_number_a = []
    galaxy_mass_a = []
    galaxy_number_q = []
    galaxy_mass_q = []
    galaxy_number_f = []
    galaxy_mass_f = []
    for a_time in range(time_step_number):
        galaxy_info_at_each_time.append([])
        galaxy_number_a.append([])
        galaxy_mass_a.append([])
        galaxy_number_q.append([])
        galaxy_mass_q.append([])
        galaxy_number_f.append([])
        galaxy_mass_f.append([])
        for a_mass in range(mass_range_number):
            galaxy_number_a[-1].append(0)
            galaxy_mass_a[-1].append(0)
            galaxy_number_q[-1].append(0)
            galaxy_mass_q[-1].append(0)
            galaxy_number_f[-1].append(0)
            galaxy_mass_f[-1].append(0)
    total_galaxy_number_at_low_z = 0
    total_galaxy_mass_at_low_z = 0
    
    # generate 1 or 2 random galaxy and calculate initial error:
    while total_galaxy_number_at_low_z == 0:
        new_galaxy = generate_a_galaxy(age_of_the_universe, 0)
        galaxy_sfh = compute_galaxy_mass_at_each_time(new_galaxy)
        galaxies.append(new_galaxy)
        galaxies_sfh_for_all_galaxy.append(galaxy_sfh)
        for at_time_i in range(time_step_number):
            log_galaxy_mass = galaxy_sfh[at_time_i][0]
            if log_galaxy_mass is not None:
                galaxy_info_at_each_time[at_time_i].append(galaxy_sfh[at_time_i])
                for a_mass_limit_index in range(mass_range_number):
                    if logMstar_list[a_mass_limit_index] < log_galaxy_mass:
                        galaxy_number_a[at_time_i][a_mass_limit_index] += 1
                        galaxy_mass = 10**log_galaxy_mass
                        galaxy_mass_a[at_time_i][a_mass_limit_index] += galaxy_mass
                        if galaxy_sfh[at_time_i][1] == 0:
                            galaxy_number_q[at_time_i][a_mass_limit_index] += 1
                            galaxy_mass_q[at_time_i][a_mass_limit_index] += galaxy_mass
                        if galaxy_sfh[at_time_i][1] == 1:
                            galaxy_number_f[at_time_i][a_mass_limit_index] += 1
                            galaxy_mass_f[at_time_i][a_mass_limit_index] += galaxy_mass
        for mass_range in range(mass_range_number):
            total_galaxy_number_at_low_z = total_galaxy_number_at_low_z + galaxy_number_a[0][mass_range]

    for mass_range in range(mass_range_number):
        total_galaxy_mass_at_low_z = total_galaxy_mass_at_low_z + galaxy_mass_a[0][mass_range]
    error = compare_with_observation(galaxy_number_a, galaxy_number_q, galaxy_number_f,
                                     galaxy_mass_a, galaxy_mass_q, galaxy_mass_f,
                                     total_galaxy_number_at_low_z, total_galaxy_mass_at_low_z)
    # output_text_file()

    # file = open('generated_galaxies.txt', 'r')
    # old_lines = file.readlines()
    # file.close()
    # # print(old_lines)
    #
    # galaxies = old_lines[1]
    # print('-----', galaxies)
    # mass = [float(x) for x in galaxies.split()]
    # print('mmm', mass)
    print('Initial error:', error)

    # # generate more galaxies but each one much improve the fit:
    # i = 111 111  # cost 1 min
    # i_total must be larger than 11111
    # Note that at least 10^6 galaxy need to be generated in order to have some massive galaxies in the sample
    i_total = 10000000
    i_total = 111111
    save_figure = False
    save_figure = True
    the_number_of_galaxies_being_plotted = min(i_total/10, 2000)
    i = i_total
    while i > 0:
    # while error > 3:
        generated_galaxy_fraction = 1 - i / i_total
        galaxy_number_a_new = []
        galaxy_number_q_new = []
        galaxy_number_f_new = []
        galaxy_mass_a_new = []
        galaxy_mass_q_new = []
        galaxy_mass_f_new = []
        for a_time in range(time_step_number):
            galaxy_number_a_new.append(galaxy_number_a[a_time].copy())
            galaxy_number_q_new.append(galaxy_number_q[a_time].copy())
            galaxy_number_f_new.append(galaxy_number_f[a_time].copy())
            galaxy_mass_a_new.append(galaxy_mass_a[a_time].copy())
            galaxy_mass_q_new.append(galaxy_mass_q[a_time].copy())
            galaxy_mass_f_new.append(galaxy_mass_f[a_time].copy())
        total_galaxy_number_at_low_z_new = 0
        total_galaxy_mass_at_low_z_new = 0
        new_galaxy = generate_a_galaxy(age_of_the_universe, generated_galaxy_fraction)
        galaxy_sfh = compute_galaxy_mass_at_each_time(new_galaxy)
        for at_time_i in range(time_step_number):
            log_galaxy_mass = galaxy_sfh[at_time_i][0]
            if log_galaxy_mass is not None:
                for a_mass_limit_index in range(mass_range_number):
                    if logMstar_list[a_mass_limit_index] < log_galaxy_mass:
                        galaxy_number_a_new[at_time_i][a_mass_limit_index] += 1
                        galaxy_mass = 10 ** log_galaxy_mass
                        galaxy_mass_a_new[at_time_i][a_mass_limit_index] += galaxy_mass
                        if galaxy_sfh[at_time_i][1] == 0:
                            galaxy_number_q_new[at_time_i][a_mass_limit_index] += 1
                            galaxy_mass_q_new[at_time_i][a_mass_limit_index] += galaxy_mass
                        if galaxy_sfh[at_time_i][1] == 1:
                            galaxy_number_f_new[at_time_i][a_mass_limit_index] += 1
                            galaxy_mass_f_new[at_time_i][a_mass_limit_index] += galaxy_mass
        for mass_range in range(mass_range_number):
            total_galaxy_number_at_low_z_new += galaxy_number_a_new[0][mass_range]
            total_galaxy_mass_at_low_z_new += galaxy_mass_a_new[0][mass_range]
        new_error = compare_with_observation(galaxy_number_a_new, galaxy_number_q_new, galaxy_number_f_new,
                                             galaxy_mass_a_new, galaxy_mass_q_new, galaxy_mass_f_new,
                                             total_galaxy_number_at_low_z_new, total_galaxy_mass_at_low_z_new)
        # add the new galaxy if fit is improved:
        if new_error < error:
            galaxies.append(new_galaxy)
            galaxies_sfh_for_all_galaxy.append(galaxy_sfh)
            galaxy_info_at_each_time[at_time_i].append(galaxy_sfh[at_time_i])
            galaxy_number_a = galaxy_number_a_new
            galaxy_number_q = galaxy_number_q_new
            galaxy_number_f = galaxy_number_f_new
            galaxy_mass_a = galaxy_mass_a_new
            galaxy_mass_q = galaxy_mass_q_new
            galaxy_mass_f = galaxy_mass_f_new
            total_galaxy_number_at_low_z = total_galaxy_number_at_low_z_new
            total_galaxy_mass_at_low_z = total_galaxy_mass_at_low_z_new
            error = new_error
        (i) = (i-1)
    print('Final error:', error)

    output_text_file()

    end = time.time()
    print("\nSimulation time:", end - start)

    plot_all_galaxy_result()