import redshift_to_time     # redshift_to_time.cosmology_calculator(redshift)
import data_Muzzin13        # data_Muzzin13.function_data(redshift, log(galaxy_stellar_mass))
import time
import random
import math


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
    if 8 < log_galaxy_final_mass < 13:
        output = [star_formation_start_time, star_formation_stop_time, log_sfr, log_galaxy_final_mass]
        return output
    else:
        return generate_a_galaxy(age_of_the_universe)


def data_prepare():
    '''
    There might be a different way which is to use the lower redshift ranges instead of middle redshifts.
    '''
    print('===\ndata_prepare...')
    global redshift_list, time_list, logMstar_list_
    redshift_list_ = data_Muzzin13.redshift_list
    logMstar_list_ = data_Muzzin13.logMstar_list
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
    print("times:", time_list)
    print("logMstar:", logMstar_list_)
    print('===\n')
    return

def compute_galaxy_mass_at_each_time(a_galaxy):
    '''
    galaxy_info = [info_time_1, info_time_2...]
    info_time_i = [None_or_log(galaxy-mass), star-forming(1)_or_quiescent(0)]
    galaxy_info = compute_galaxy_mass_at_each_time(a_galaxy)
    '''
    galaxy_info = []
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
        galaxy_info.append([None_or_log_galaxy_mass, star_forming_or_quiescent])
    return galaxy_info

def compare_with_observation(galaxy_info_at_each_time__):
    # galaxy_list = []
    # for ii in range(len(time_list)):
    #     galaxy_list.append([])
    #     for jj in range(len(redshift_list)):
    #         galaxy_list[ii].append([])

    # error_list = []
    # for ii in range(len(time_list)):
    #     time__ = time_list[ii]
    #     error_list.append([])
    #     for jj in range(len(redshift_list)):
    #         redshift__ = redshift_list[jj]
    #         error_list[ii].append([])

    # galaxy_number = len(galaxies)
    # while galaxy_number > 0:
    #     a_galaxy = galaxies[galaxy_number]
    #     for time__ in time_list:
    #         if a_galaxy[0] < time__:
    #             if a_galaxy[1] < time__:
    #
    #     (galaxy_number) = (galaxy_number-1)
    error__ = 0

    return error__


def plot_all_galaxy_result():
    print("galaxies:", galaxies)
    print("galaxies_info_list:", galaxies_info_list)
    print("galaxy_info_at_each_time:", galaxy_info_at_each_time)
    return


if __name__ == '__main__':
    start = time.time()

    age_of_the_universe = redshift_to_time.cosmology_calculator(0)
    print("age_of_the_universe:", age_of_the_universe)
    data_prepare()

    galaxies = []
    galaxies_info_list = []
    galaxy_info_at_each_time = []
    galaxy_number = []
    galaxy_total_mass = []
    for a_time in range(len(time_list)):
        galaxy_info_at_each_time.append([])
        galaxy_number.append([])
        galaxy_total_mass.append([])
        for a_mass in range(len(logMstar_list_)-1):
            galaxy_number[-1].append([])
            galaxy_total_mass[-1].append([])

    print("galaxy_number", galaxy_number)
    print("galaxy_total_mass", galaxy_total_mass)

    # generate 10 random galaxy and calculate initial error:
    i = 2
    while i > 0:
        new_galaxy = generate_a_galaxy(age_of_the_universe)
        galaxy_info = compute_galaxy_mass_at_each_time(new_galaxy)
        galaxies.append(new_galaxy)
        galaxies_info_list.append(galaxy_info)
        for info_at_time_i in range(len(galaxy_info)):
            if galaxy_info[info_at_time_i][0] is not None:
                galaxy_info_at_each_time[info_at_time_i].append(galaxy_info[info_at_time_i])

        (i) = (i - 1)

    error = compare_with_observation(galaxy_info_at_each_time)

    # # generate more galaxies but each one much improve the fit:
    # i = 2
    # while i > 0:
    #     new_galaxy = generate_a_galaxy(age_of_the_universe)
    #     galaxy_info = compute_galaxy_mass_at_each_time(new_galaxy)
    #     galaxy_info_at_each_time_new = galaxy_info_at_each_time
    #     for info_at_time_i in range(len(galaxy_info)):
    #         if galaxy_info[info_at_time_i][0] is not None:
    #             galaxy_info_at_each_time_new[info_at_time_i].append(galaxy_info[info_at_time_i])
    #     new_error = compare_with_observation(galaxy_info_at_each_time_new)
    #     # add the new galaxy if fit is improved:
    #     if new_error < error:
    #         galaxies.append(new_galaxy)
    #         galaxies_info_list.append(galaxy_info)
    #         galaxy_info_at_each_time = galaxy_info_at_each_time_new
    #         error = new_error
    #     (i) = (i-1)



    plot_all_galaxy_result()

    end = time.time()
    print("Run time:", end - start)
