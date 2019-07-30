import redshift_to_time     # redshift_to_time.cosmology_calculator(redshift)
import data_Muzzin13        # data_Muzzin13.function_data(redshift, log(galaxy_stellar_mass))
import time
import random
import math
import matplotlib.pyplot as plt


def generate_a_galaxy(age_of_the_universe):
    '''
    :param age_of_the_universe:
    :return: [star_formation_start_time, star_formation_stop_time, log_sfr]
    '''
    star_formation_start_time = random.random() * time_list[0]  # in [Gyr]
    maximum_formation_time = age_of_the_universe - star_formation_start_time
    star_formation_stop_time = star_formation_start_time + random.random() * maximum_formation_time  # in [Gyr]
    log_sfr = random.random() * 10 - 5  # log_10(star formation rate [solar mass per year])
    log_galaxy_final_mass = 9 + math.log((star_formation_stop_time - star_formation_start_time), 10) + log_sfr
    # Here does not consider that
    # the final dynamical mass should be smaller than the total stellar mass formed
    if 8 < log_galaxy_final_mass < 13:
        star_formation_timescale = star_formation_stop_time - star_formation_start_time
        output = [star_formation_start_time, star_formation_stop_time, log_sfr,
                  log_galaxy_final_mass, star_formation_timescale]
        return output
    else:
        return generate_a_galaxy(age_of_the_universe)


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
        time_ = redshift_to_time.cosmology_calculator(redshift__)
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
    error__ = 0
    for ii in range(time_step_number):
        for jj in range(mass_range_number):
            error__ += abs(galaxy_number_a__[ii][jj]/total_galaxy_number_at_low_z__ - \
                       data_Muzzin13.function_data(ii, jj)[0]/total_galaxy_number_at_low_z_obs)
            error__ += abs(galaxy_number_q__[ii][jj]/total_galaxy_number_at_low_z__ - \
                       data_Muzzin13.function_data(ii, jj)[1]/total_galaxy_number_at_low_z_obs)
            error__ += abs(galaxy_number_f__[ii][jj]/total_galaxy_number_at_low_z__ - \
                       data_Muzzin13.function_data(ii, jj)[2]/total_galaxy_number_at_low_z_obs)
            error__ += abs(galaxy_mass_a__[ii][jj]/total_galaxy_mass_at_low_z__ - \
                       data_Muzzin13.function_data(ii, jj)[3]/total_galaxy_mass_at_low_z_obs)
            error__ += abs(galaxy_mass_q__[ii][jj]/total_galaxy_mass_at_low_z__ - \
                       data_Muzzin13.function_data(ii, jj)[4]/total_galaxy_mass_at_low_z_obs)
            error__ += abs(galaxy_mass_f__[ii][jj]/total_galaxy_mass_at_low_z__ - \
                       data_Muzzin13.function_data(ii, jj)[5]/total_galaxy_mass_at_low_z_obs)
    return error__


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

    plot_galaxy_mass_list = []
    plot_galaxy_sft_list = []
    for a_galaxy in galaxies:
        plot_galaxy_mass_list.append(a_galaxy[3])
        plot_galaxy_sft_list.append(a_galaxy[4])

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.figure(0, figsize=(6, 5))
    plt.scatter(plot_galaxy_mass_list, plot_galaxy_sft_list, s=3)
    plt.xlabel(r'log$_{10}$(M$_{dyn}$ [M$_\odot$])')
    plt.ylabel(r't$_{\rm sf}$')

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.figure(1, figsize=(6, 5))
    plt.xlabel(r'Time')
    plt.ylabel(r'SFR')
    for a_galaxy in galaxies:
        star_formation_start_time = a_galaxy[0]
        star_formation_stop_time = a_galaxy[1]
        log_sfr = a_galaxy[2]
        x = [star_formation_start_time, star_formation_stop_time, star_formation_stop_time, star_formation_start_time]
        y = [log_sfr, log_sfr, -5, -5]
        plt.fill(x, y, 'b', alpha=0.0333, edgecolor='r')
    plt.tight_layout()
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

    age_of_the_universe = redshift_to_time.cosmology_calculator(0)
    print("age_of_the_universe:", age_of_the_universe)
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
    
    # generate 1 random galaxy and calculate initial error:
    i = 1
    while i > 0:
        new_galaxy = generate_a_galaxy(age_of_the_universe)
        galaxy_sfh = compute_galaxy_mass_at_each_time(new_galaxy)
        galaxies.append(new_galaxy)
        galaxies_sfh_for_all_galaxy.append(galaxy_sfh)
        for at_time_i in range(time_step_number):
            log_galaxy_mass = galaxy_sfh[at_time_i][0]
            if log_galaxy_mass is not None:
                galaxy_info_at_each_time[at_time_i].append(galaxy_sfh[at_time_i])
                for a_mass_limit_index in range(mass_range_number):
                    if logMstar_list[a_mass_limit_index] < log_galaxy_mass < logMstar_list[a_mass_limit_index + 1]:
                        galaxy_number_a[at_time_i][a_mass_limit_index] += 1
                        galaxy_mass = 10**log_galaxy_mass
                        galaxy_mass_a[at_time_i][a_mass_limit_index] += galaxy_mass
                        if galaxy_sfh[at_time_i][1] == 0:
                            galaxy_number_q[at_time_i][a_mass_limit_index] += 1
                            galaxy_mass_q[at_time_i][a_mass_limit_index] += galaxy_mass
                        if galaxy_sfh[at_time_i][1] == 1:
                            galaxy_number_f[at_time_i][a_mass_limit_index] += 1
                            galaxy_mass_f[at_time_i][a_mass_limit_index] += galaxy_mass
        (i) = (i - 1)

    for mass_range in range(mass_range_number):
        total_galaxy_number_at_low_z = total_galaxy_number_at_low_z + galaxy_number_a[0][mass_range]
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
    i = 11111
    while i > 0:
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
        new_galaxy = generate_a_galaxy(age_of_the_universe)
        galaxy_sfh = compute_galaxy_mass_at_each_time(new_galaxy)
        for at_time_i in range(time_step_number):
            log_galaxy_mass = galaxy_sfh[at_time_i][0]
            if log_galaxy_mass is not None:
                for a_mass_limit_index in range(mass_range_number):
                    if logMstar_list[a_mass_limit_index] < log_galaxy_mass < logMstar_list[a_mass_limit_index + 1]:
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