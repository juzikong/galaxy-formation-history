import redshift_to_time     # redshift_to_time.cosmology_calculator(redshift)
import data_Muzzin13        # data_Muzzin13.function_data(redshift, log(galaxy_stellar_mass))
import time
import random


def generate_a_galaxy(age_of_the_universe):
    star_formation_start_time = random.random() * age_of_the_universe
    maximum_formation_time = age_of_the_universe - star_formation_start_time
    star_formation_stop_time = star_formation_start_time + random.random() * maximum_formation_time
    log_sfr = random.random() * 10 - 5  # log_10(star formation rate [solar mass per year])
    output = [star_formation_start_time, star_formation_stop_time, log_sfr]
    return output


def compare_with_observation():
    error__ = 0
    return error__


def plot_all_galaxy_result():
    print(galaxies)
    return


galaxies = [[0, 0, 0]]

if __name__ == '__main__':
    start = time.time()

    age_of_the_universe = redshift_to_time.cosmology_calculator(0)
    print(age_of_the_universe)

    error = compare_with_observation()

    i = 3
    while i > 0:
        new_galaxy = generate_a_galaxy(age_of_the_universe)
        galaxies.append(new_galaxy)
        new_error = compare_with_observation()
        if new_error < error:
            del galaxies[-1]
        (i) = (i-1)

    plot_all_galaxy_result()

    end = time.time()
    print("Run time:", end - start)
