def function_data(Z, logMstar):
    '''
        This function returns the number (and mass) of galaxies (stellar mass) within given mass and redshift range
        redshift range  Z = 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 4
        mass range      logMstar = 8, 9, 10, 11, 11.5, 13
        Data from Muzzin et al. (2013) their Table. 2
    '''
    log_eta_a = [[-1.49, 1.89, 2.35, 3.25, 4.52],
                 [-1.69, 2.00, 2.45, 3.49, 4.88],
                 [-2.04, 2.26, 2.66, 3.77, 5.45],
                 [-2.41, 2.56, 2.88, 3.97, 5.76],
                 [-3.06, 3.11, 3.30, 4.20, 5.89],
                 [-2.88, 3.09, 3.43, 4.32, 5.57],
                 [-3.14, 3.65, 4.20, 5.04, 5.80]]
    log_eta_q = [[-2.19, 2.45, 2.66, 3.33, 4.60], 
                 [-2.36, 2.70, 2.85, 3.60, 5.10], 
                 [-3.01, 3.01, 3.12, 3.93, 5.74], 
                 [-3.39, 3.40, 3.48, 4.29, 6.27], 
                 [-3.58, 3.67, 3.91, 4.82, 6.37], 
                 [-4.02, 4.05, 4.19, 5.00, 6.63], 
                 [-5.44, 5.44, 5.46, 5.83, 6.97]]
    log_eta_f = [[-1.58, 2.03, 2.65, 4.08, 6.07],
                 [-1.70, 2.09, 2.66, 4.10, 6.18],
                 [-1.94, 2.30, 2.85, 4.27, 6.39],
                 [-2.21, 2.52, 3.01, 4.25, 6.05],
                 [-3.18, 3.23, 3.42, 4.32, 6.04],
                 [-2.99, 3.20, 3.55, 4.42, 5.61],
                 [-3.11, 3.67, 4.27, 5.11, 5.85]]
    log_rho_a = [[8.41, 8.40, 8.35, 7.98, 7.08],
                 [8.26, 8.25, 8.19, 7.73, 6.72],
                 [8.02, 8.01, 7.95, 7.41, 6.13],
                 [7.79, 7.79, 7.74, 7.20, 5.81],
                 [7.43, 7.43, 7.41, 6.98, 5.69],
                 [7.32, 7.32, 7.28, 6.91, 6.04],
                 [6.64, 6.63, 6.57, 6.31, 5.89]]
    log_rho_q = [[8.19, 8.19, 8.18, 7.91, 7.00],
                 [7.96, 7.96, 7.94, 7.61, 6.48],
                 [7.65, 7.65, 7.64, 7.25, 5.83],
                 [7.29, 7.30, 7.29, 6.88, 5.29],
                 [6.83, 6.82, 6.80, 6.38, 5.21],
                 [6.58, 6.58, 6.57, 6.19, 4.95],
                 [5.58, 5.58, 5.58, 5.44, 4.62]]
    log_rho_f = [[7.99, 7.97, 7.85, 7.08, 5.50],
                 [7.97, 7.95, 7.84, 7.06, 5.39],
                 [7.78, 7.76, 7.66, 6.88, 5.17],
                 [7.65, 7.64, 7.56, 6.93, 5.53],
                 [7.31, 7.31, 7.29, 6.86, 5.53],
                 [7.21, 7.21, 7.17, 6.82, 6.00],
                 [6.59, 6.58, 6.50, 6.25, 5.86]]
    '''
        The units are:  eta [Mpc^-3]    rho [solar mass Mpc^-3]
        a, q, and f meaning: all, quiescent, and star-forming galaxies are taken into account
    '''
    output = [log_eta_a[Z][logMstar], log_eta_q[Z][logMstar], log_eta_f[Z][logMstar],
              log_rho_a[Z][logMstar], log_rho_q[Z][logMstar], log_rho_f[Z][logMstar]]
    return output


print(function_data(1, 1))
