import numpy as np
from pyXSteam.XSteam import XSteam
import matplotlib.pyplot as plt
import scipy.special as sp

steam_table = XSteam(XSteam.UNIT_SYSTEM_BARE) # m, kg, second, Kelvin, MPa, W, KJ

# Input parameters

N_fa = 157  # number of fuel assemblies
N_fr = 264  # number of fuel rods per assembly
N_cell = 101    # number of cells per fuel rod
N_spacer = 10   # number of spacers
Q_tot_actual = 3400e6 * 1.05  # total thermal power, MW
H = 4.267   # height of the fuel rod, m
l_cell = H / N_cell # length of each cell
d_r = 0.0095    # diameter of the fuel rod, m
r_co = d_r/2    # cladding outer radius
r_ci = r_co - 0.000572 # cladding inner radius
r_fo = 0.0081915 / 2    # fuel outer radius
r_f2 = (2/3) * r_fo
r_f1 = (1/3) * r_fo
pitch = 0.0126  # pitch of the fuel assemblies
P_operating = 15.513 # operating pressure, MPa
T_inlet = 280.7 + 273.15 # inlet temperature, Kelvin
T_outlet = 324.4 + 273.15 # outlet temperature, Kelvin
p_inlet = P_operating # inlet pressure, MPa
dp_inletorifice = - 0.25 * p_inlet # the pressure drop due to orifice
W_actual = 14300    # mass flow rate, kg/s
W_eff = W_actual * 0.95 # mass flow rate, kg/s
W_channel = W_eff / (N_fa * N_fr)    # mass flow rate through a single channel
H_by_H_e = 5 / 6
H_e = H / H_by_H_e  # extrapolated height
R_by_R_e = 5 / 6    # ... extrapolated radius
clad_roughness = 0.0025e-3   # cladding roughness
i_inlet = steam_table.h_pt(P_operating, T_inlet) # inlet enthalpy, KJ/kg

f_z = (np.pi * H_by_H_e) / (2 * np.sin((np.pi / 2) * H_by_H_e))
f_r = (2.405 * R_by_R_e) / (2 * sp.jv(1, (2.405 * R_by_R_e)))


# cell axial positions, m

cell_length = H / N_cell
z_cells = np.linspace(- H/2 + cell_length/2, H/2 - cell_length/2, N_cell)


# local loss coefficients

xi = np.zeros(N_cell)

# inlet + orifices
xi[0] = 0.5 + 6.096725342605476

# outlet
xi[-1] = 1.0

# spacer cells
spacer_cells = np.zeros(N_spacer)
for j in range(1, N_spacer + 1):
    spacer_cells[j - 1] = np.ceil(N_cell * j / (N_spacer + 1)) - 1
    
# spacer pressure loss coefficient
xi_spacer = 0.7

for j in spacer_cells:
    xi[int(j)] = xi_spacer


# cell power, both axial and radial

q_rodAve = Q_tot_actual / (N_fa * N_fr)
q_cellAve = q_rodAve / N_cell

# average channel
def average_q_cell_z(z):
    return q_cellAve * f_z * np.cos(np.pi * z / H_e)

# hot channel
def hot_q_cell_z(z):
    return q_cellAve * f_r * f_z * np.cos(np.pi * z / H_e)


# common calculations for the channel
D_h = d_r * ((4 / np.pi) * (pitch / d_r)**2 - 1)    # hydraulic diameter
P_w = np.pi * d_r     # m, wetted perimeter
A_h = D_h * P_w / 4     # m^2, flow area


# power and flow rate percentages
power_percs = np.array([84])
flow_percs = np.arange(100, 100 + 1, 1)


# array to hold total pressure drop data at different power levels
total_pressure_drop = np.zeros((len(power_percs), len(flow_percs)))

# array to hold quality data at different power levels, flow rate is the last one
x_eq_plot = np.zeros((len(power_percs), len(z_cells)))


# # iterate over power percentages (average channel)
# for q_perc_index in range(0, len(power_percs)):
    
#     print(f'{power_percs[q_perc_index]}% power')
    
#     # iterate over mass flow rate percentages
#     for w_perc_index in range(0, len(flow_percs)):
#         # print(f'{flow_percs[w_perc_index]}% mass flow rate')
        
#         # calculate mass flow rate through channel
#         W = W_channel * flow_percs[w_perc_index] / 100
        
#         G = W / A_h   # mass flux through channel
#         G_R = G / (1356.23)
        
#         # power distribution along the cells
#         q_cells = average_q_cell_z(z_cells) * power_percs[q_perc_index] / 100
#         q_cells = q_cells * q_rodAve / np.sum(q_cells)

#         # enthalpy distribution along the cells
#         i = np.ones(N_cell) * i_inlet
#         for j in range(1, N_cell):
#             i[j] = i[j - 1] + (q_cells[j - 1] / 1e3) / W
            
#         # arrays to hold different pressure drops along the cells
#         dpf = np.zeros(N_cell)
#         dpg = np.zeros(N_cell)
#         dpl = np.zeros(N_cell)
#         dpa = np.zeros(N_cell)
#         dp  = np.zeros(N_cell)
        
#         # array to hold pressure values along the cells
#         p = np.zeros(N_cell)
#         p[0] = p_inlet
        
#         # array to hold equilibrium quality
#         x_eq = np.zeros(N_cell)
        
#         # array to hold temperature
#         T = np.zeros(N_cell)    # T_lb
#         T_co = np.zeros(N_cell)
#         T_ci = np.zeros(N_cell)
#         T_fo = np.zeros(N_cell)
#         T_f2 = np.zeros(N_cell)
#         T_f1 = np.zeros(N_cell)
#         T_fc = np.zeros(N_cell)
#         T[0] = T_inlet
        
#         # array to hold heat flux
#         q2 = np.zeros(N_cell)
        
#         # array to hold critical heat flux
#         q2_cr = np.zeros(N_cell)
        
#         # array to hold DNBR values
#         DNBR = np.zeros(N_cell)
        
#         # calculating parameters at the inlet
#         rho_inlet = steam_table.rho_ph(P_operating, i[0])
#         dpl[0] = xi[0] * G**2 / (2 * rho_inlet)
        
#         i_l_inlet = steam_table.hL_p(p[0])
#         i_g_inlet = steam_table.hV_p(p[0])
#         x_eq[0] = (i[0] - i_l_inlet) / (i_g_inlet - i_l_inlet)
        
#         a_1 = 0.5328
#         a_2 = 0.1212
#         a_3 = -0.3040
#         a_4 = 0.3285
        
#         c_1 = 1.6151
#         c_2 = 1.4066
#         c_3 = 0.4843
#         c_4 = -2.0749
        
#         P_H = np.pi * d_r # heated perimeter
#         # D_H = d_r * ((4 / np.pi) * (pitch / d_r)**2 - 1)    # heated diameter
#         A_H = P_H * cell_length # D_H * P_H / 4 # heated area
        
#         # CHF calculation at the first cell
#         p_R = p[0] / 22 # steam_table.psat_t(T[0])
        
#         A = a_1 * p_R**a_2 * G_R**(a_3 + a_4 * p_R)
#         B = 3.1544e6
#         C = c_1 * p_R**c_2 * G_R**(c_3 + c_4 * p_R)
        
#         q1 = q_cells / cell_length
#         q2 = q_cells / A_H
#         q3 = q_cells / (np.pi * r_fo**2 * cell_length)
#         q2_R_0 = q2[0] / B
            
#         q2_cr[0] = B * (A - x_eq[0]) / (C)
        
#         DNBR[0] = q2_cr[0] / q2[0]
        
#         # temperature calculation for the first cell
        
#         # liquid to cladding
        
#         # Reynold's number
#         mu = np.zeros(N_cell)
#         Re = np.zeros(N_cell)
#         rho_0 = steam_table.rho_ph(p[0], i[0]) # kg/m^3
#         mu[0] = steam_table.my_ph(p[0], i[0]) # Pa.s
#         Re[0] = G * D_h / mu[0]
        
#         # array to hold Prandtl's number
#         Pr = np.zeros(N_cell)
#         Cp = steam_table.Cp_ph(p[0], i[0]) * 1e3    # joules per kg per kelvin
#         lambda_l = steam_table.tc_ph(p[0], i[0])
#         Pr[0] = Cp * mu[0] / lambda_l
        
#         # array to hold Nusselt's number
#         Nu = np.zeros(N_cell)
#         Nu[0] = 0.023 * Re[0]**0.8 * Pr[0]**0.4
        
#         # array to hold heat transfer coefficient
#         h = np.zeros(N_cell)
#         h[0] = Nu[0] * lambda_l / D_h
        
#         T_co[0] = T[0] + q2[0] / h[0]
        
#         # cladding
#         T_ci[0] = T_co[0] + 1
#         for k in range(0, 3):
#             T_cAvg = (T_ci[0] + T_co[0]) / 2 - 273.15 # celcius
#             lambda_clad = 12.6 + 0.0118 * T_cAvg
#             delta_T_c = (q3[0] * r_fo**2 / (2 * lambda_clad)) * np.log(r_co / r_ci)
#             T_ci[0] = T_co[0] + delta_T_c
        
#         # gas gap    
#         lambda_gas_gap = 0.6
#         # delta_T_gas_gap = (q1[0] / 4 * np.pi) * (2 * np.log(r_ci / r_fo) / lambda_gas_gap)
#         delta_T_gas_gap = (q3[0] * r_fo**2 / (2 * lambda_gas_gap)) * np.log(r_ci / r_fo)
#         T_fo[0] = T_ci[0] + delta_T_gas_gap
        
#         # fuel
#         T_f2[0] = T_fo[0] + 1
#         for k in range(0, 3):
#             T_f2Avg = (T_f2[0] + T_fo[0]) / 2 # kelvin
#             t_f2 = T_f2Avg / 1000
#             lambda_f2 = 100 / (7.5408 + 17.692 * t_f2 + 3.6142 * t_f2**2) + 6400 * np.exp(-16.35/t_f2) / t_f2**(5/2)
#             delta_T_f2 = (q3[0] * (r_fo**2 - r_f2**2)) / (4 * lambda_f2)
#             T_f2[0] = T_fo[0] + delta_T_f2
        
#         T_f1[0] = T_f2[0] + 1
#         for k in range(0, 3):
#             T_f1Avg = (T_f1[0] + T_f2[0]) / 2 # kelvin
#             t_f1 = T_f1Avg / 1000
#             lambda_f1 = 100 / (7.5408 + 17.692 * t_f1 + 3.6142 * t_f1**2) + 6400 * np.exp(-16.35/t_f1) / t_f1**(5/2)
#             delta_T_f1 = (q3[0] * (r_f2**2 - r_f1**2)) / (4 * lambda_f1)
#             T_f1[0] = T_f2[0] + delta_T_f1
            
#         T_fc[0] = T_fc[0] + 1
#         for k in range(0, 3):
#             T_fcAvg = (T_fc[0] + T_f1[0]) / 2 # kelvin
#             t_fc = T_fcAvg / 1000
#             lambda_fc = 100 / (7.5408 + 17.692 * t_fc + 3.6142 * t_fc**2) + 6400 * np.exp(-16.35/t_fc) / t_fc**(5/2)
#             delta_T_fc = (q3[0] * r_f1**2) / (4 * lambda_fc)
#             T_fc[0] = T_f1[0] + delta_T_fc
        
#         # iterate over all the cells
#         for j in range(1, N_cell):
#             # calculate pressure at the cell
#             p[j] = p_inlet - dp[j - 1] / 1e6   # dp in Pa, p in MPa
            
#             # cell temperatures
#             T[j] = steam_table.t_ph(p[j], i[j])
            
#             # q2[j] = q_cells[j] / A_H
            
#             # calculate quality at the cell
#             i_g = steam_table.hV_p(p[j])
#             i_l = steam_table.hL_p(p[j])
#             x_eq[j] = (i[j] - i_l) / (i_g - i_l)
            
#             if x_eq[j] <= 0:
#                 rho = steam_table.rho_ph(p[j], i[j]) # kg/m^3
#                 mu[j] = steam_table.my_ph(p[j], i[j]) # Pa.s
#                 Re[j] = G * D_h / mu[j]
#                 C_f = 1 / (- 3.6 * np.log10(((clad_roughness / D_h) / 3.7)**1.11 + 6.9 / Re[j]))**2

#                 # pressure drops from the inlet to the cell
#                 dpf[j] = dpf[j - 1] + 4 * C_f * cell_length * G**2 / (D_h * 2 * rho) # Pa
#                 dpg[j] = dpg[j - 1] + cell_length * rho * 9.81 * 1 # Pa
#                 dpl[j] = dpl[j - 1] + xi[j] * G**2 / (2 * rho) # Pa
                
#                 Cp = steam_table.Cp_ph(p[j], i[j]) * 1e3    # joules per kg per kelvin
#                 lambda_l = steam_table.tc_ph(p[j], i[j])
#                 Pr[j] = Cp * mu[j] / lambda_l
#                 Nu[j] = 0.023 * Re[j]**0.8 * Pr[j]**0.4
#                 h[j] = Nu[j] * lambda_l / D_h
#                 T_co[j] = T[j] + q2[j] / h[j]
                            
#             elif x_eq[j] > 0:
#                 if x_eq[j] > 1:
#                     x_eq[j] = 1
#                 rho_l = steam_table.rhoL_p(p[j])
#                 rho_g = steam_table.rhoV_p(p[j])
#                 T_sat = steam_table.tsat_p(p[j])
#                 mu_l = steam_table.my_pt(p[j], T_sat - 0.01)
#                 mu_g = steam_table.my_pt(p[j], T_sat + 0.01)
#                 rho_mix = rho_l / (x_eq[j] * (rho_l / rho_g - 1) + 1)
#                 # mu_mix = x_eq[j] * mu_g + (1 - x_eq[j - 1]) * mu_l
#                 mu_mix = 1 / (x_eq[j] / mu_g + (1 - x_eq[j]) / mu_l)
#                 mu[j] = mu_mix
#                 Re[j] = G * D_h / mu[j]
#                 C_f = 1 / (- 3.6 * np.log10(((clad_roughness / D_h) / 3.7)**1.11 + 6.9 / Re[j]))**2
                
#                 # pressure drops from the inlet to the cell
#                 dpf[j] = dpf[j - 1] + 4 * C_f * cell_length * G**2 / (D_h * 2 * rho) # Pa
#                 dpg[j] = dpg[j - 1] + cell_length * rho_mix * 9.81 * 1 # Pa
#                 dpl[j] = dpl[j - 1] + xi[j] * G**2 / (2 * rho_mix) # Pa
#                 # dpa[j] = dpa[j - 1] + G**2 / (rho)  # Pa
#                 dpa[j] = dpa[j - 1] + G**2 * (1/rho_g - 1/rho_l) * (x_eq[j] - x_eq[j - 1])
                
#                 # chen correlation
#                 X_tt = ((1 - x_eq[j]) / x_eq[j])**0.9 * (rho_g / rho_l)**0.5 * (mu_l / mu_g)**0.1
#                 Re_chen = G * (1 - x_eq[j]) * D_h / mu_l
#                 if 1/X_tt <= 0.1:
#                     F = 1
#                 else:
#                     F = 2.35 * (0.213 + 1/X_tt)**0.736
#                 S = 1 / (1 + 2.56 * 1e-6 * F**1.463 * Re_chen**1.17)
#                 lambda_satl = steam_table.tcL_p(p[j])
#                 Cp_satl = steam_table.CpL_p(p[j]) * 1e3
#                 Pr_f = Cp_satl * mu_l / lambda_satl
#                 h_mac = 0.023 * (lambda_satl / D_h) * Re_chen**0.8 * Pr_f**0.4 * F
                
#                 delta_T_sup = 5
#                 T_w = T_sat + delta_T_sup
#                 p_s = steam_table.psat_t(T_w)
#                 sigma = steam_table.st_p(p[j])
                
#                 for k in range(0, 3):
#                     h_mic = 0.00122 * ((lambda_satl**0.79 * Cp_satl**0.45 * rho_l**0.49) / (sigma**0.5 * mu_l**0.29 * (i_g - i_l)**0.24 * rho_g**0.24)) * delta_T_sup**0.24 * (p_s - p[j])**0.75 * S
#                     h[j] = h_mic + h_mac
#                     delta_T_sup_new = q2[j] / h[j]
#                     delta_T_sup = delta_T_sup_new
                
#                 delta_T_lb = delta_T_sup
#                 T_co[j] = T[j] + delta_T_lb
                    
#             # cladding temperature calculation
#             T_ci[j] = T_co[j] + 1
#             for k in range(0, 3):
#                 T_cAvg = (T_ci[j] + T_co[j]) / 2 - 273.15 # celcius
#                 lambda_clad = 12.6 + 0.0118 * T_cAvg
#                 delta_T_c = (q3[j] * r_fo**2 / (2 * lambda_clad)) * np.log(r_co / r_ci)
#                 T_ci[j] = T_co[j] + delta_T_c
            
#             # gas gap temperature calculation
#             lambda_gas_gap = 0.6
#             delta_T_gas_gap = (q3[j] * r_fo**2 / (2 * lambda_gas_gap)) * np.log(r_ci / r_fo)
#             T_fo[j] = T_ci[j] + delta_T_gas_gap
            
#             # fuel temperature calculation
#             T_f2[j] = T_fo[j] + 1
#             for k in range(0, 3):
#                 T_f2Avg = (T_f2[j] + T_fo[j]) / 2 # kelvin
#                 t_f2 = T_f2Avg / 1000
#                 lambda_f2 = 100 / (7.5408 + 17.692 * t_f2 + 3.6142 * t_f2**2) + 6400 * np.exp(-16.35/t_f2) / t_f2**(5/2)
#                 delta_T_f2 = (q3[j] * (r_fo**2 - r_f2**2)) / (4 * lambda_f2)
#                 T_f2[j] = T_fo[j] + delta_T_f2
            
#             T_f1[j] = T_f2[j] + 1
#             for k in range(0, 3):
#                 T_f1Avg = (T_f1[j] + T_f2[j]) / 2 # kelvin
#                 t_f1 = T_f1Avg / 1000
#                 lambda_f1 = 100 / (7.5408 + 17.692 * t_f1 + 3.6142 * t_f1**2) + 6400 * np.exp(-16.35/t_f1) / t_f1**(5/2)
#                 delta_T_f1 = (q3[j] * (r_f2**2 - r_f1**2)) / (4 * lambda_f1)
#                 T_f1[j] = T_f2[j] + delta_T_f1
                
#             T_fc[j] = T_fc[j] + 1
#             for k in range(0, 3):
#                 T_fcAvg = (T_fc[j] + T_f1[j]) / 2 # kelvin
#                 t_fc = T_fcAvg / 1000
#                 lambda_fc = 100 / (7.5408 + 17.692 * t_fc + 3.6142 * t_fc**2) + 6400 * np.exp(-16.35/t_fc) / t_fc**(5/2)
#                 delta_T_fc = (q3[j] * r_f1**2) / (4 * lambda_fc)
#                 T_fc[j] = T_f1[j] + delta_T_fc
                        
#             # total pressure drop from the inlet to the cell
#             dp[j] = dpf[j] + dpg[j] + dpl[j] + dpa[j]
            
#             # CHF calculation
#             p_R = p[j] / 22 # steam_table.psat_t(T_inlet)
            
#             A = a_1 * p_R**a_2 * G_R**(a_3 + a_4 * p_R)
#             B = 3.1544e6
#             C = c_1 * p_R**c_2 * G_R**(c_3 + c_4 * p_R)
            
#             q2_R = q2[j] / B
            
#             q2_cr[j] = B * (A - x_eq[0]) / (C + (x_eq[j] - x_eq[0]) / q2_R)
            
#             DNBR[j] = q2_cr[j] / q2[j]
            
#         total_pressure_drop[q_perc_index, w_perc_index] = dp[N_cell - 1]
        
#     x_eq_plot[q_perc_index,] = x_eq

# average = i

# iterate over power percentages (hot channel)
for q_perc_index in range(0, len(power_percs)):
    
    print(f'{power_percs[q_perc_index]}% power')
    
    # iterate over mass flow rate percentages
    for w_perc_index in range(0, len(flow_percs)):
        print(f'{flow_percs[w_perc_index]}% mass flow rate')
        
        # calculate mass flow rate through channel
        W = W_channel * flow_percs[w_perc_index] / 100
        
        G = W / A_h   # mass flux through channel
        G_R = G / (1356.23)
        
        # power distribution along the cells
        q_cells = hot_q_cell_z(z_cells)     # calculate cell powers
        q_cells = q_cells * q_rodAve * f_r / np.sum(q_cells)    # correct it
        q_cells = hot_q_cell_z(z_cells) * power_percs[q_perc_index] / 100   # adjust for raised/lowered power percs

        # enthalpy distribution along the cells
        i = np.ones(N_cell) * i_inlet
        for j in range(1, N_cell):
            i[j] = i[j - 1] + (q_cells[j - 1] / 1e3) / W
            
        # arrays to hold different pressure drops along the cells
        dpf = np.zeros(N_cell)
        dpg = np.zeros(N_cell)
        dpl = np.zeros(N_cell)
        dpa = np.zeros(N_cell)
        dp  = np.zeros(N_cell)
        
        # array to hold pressure values along the cells
        p = np.zeros(N_cell)
        p[0] = p_inlet
        
        # array to hold equilibrium quality
        x_eq = np.zeros(N_cell)
        
        # array to hold temperature
        T = np.zeros(N_cell)    # T_lb
        T_co = np.zeros(N_cell)
        T_ci = np.zeros(N_cell)
        T_fo = np.zeros(N_cell)
        T_f2 = np.zeros(N_cell)
        T_f1 = np.zeros(N_cell)
        T_fc = np.zeros(N_cell)
        T[0] = T_inlet
        
        # array to hold heat flux
        q2 = np.zeros(N_cell)
        
        # array to hold critical heat flux
        q2_cr = np.zeros(N_cell)
        
        # array to hold DNBR values
        DNBR = np.zeros(N_cell)
        
        # calculating parameters at the inlet
        rho_inlet = steam_table.rho_ph(P_operating, i[0])
        dpl[0] = xi[0] * G**2 / (2 * rho_inlet)
        
        i_l_inlet = steam_table.hL_p(p[0])
        i_g_inlet = steam_table.hV_p(p[0])
        x_eq[0] = (i[0] - i_l_inlet) / (i_g_inlet - i_l_inlet)
        
        a_1 = 0.5328
        a_2 = 0.1212
        a_3 = -0.3040
        a_4 = 0.3285
        
        c_1 = 1.6151
        c_2 = 1.4066
        c_3 = 0.4843
        c_4 = -2.0749
        
        P_H = np.pi * d_r # heated perimeter
        # D_H = d_r * ((4 / np.pi) * (pitch / d_r)**2 - 1)    # heated diameter
        A_H = P_H * cell_length # D_H * P_H / 4 # heated area
        
        # CHF calculation at the first cell
        p_R = p[0] / 22 # steam_table.psat_t(T[0])
        
        A = a_1 * p_R**a_2 * G_R**(a_3 + a_4 * p_R)
        B = 3.1544e6
        C = c_1 * p_R**c_2 * G_R**(c_3 + c_4 * p_R)
        
        q1 = q_cells / cell_length
        q2 = q_cells / A_H
        q3 = q_cells / (np.pi * r_fo**2 * cell_length)
        q2_R_0 = q2[0] / B
            
        q2_cr[0] = B * (A - x_eq[0]) / (C)
        
        DNBR[0] = q2_cr[0] / q2[0]
        
        # temperature calculation for the first cell
        
        # liquid to cladding
        
        # Reynold's number
        mu = np.zeros(N_cell)
        Re = np.zeros(N_cell)
        rho_0 = steam_table.rho_ph(p[0], i[0]) # kg/m^3
        mu[0] = steam_table.my_ph(p[0], i[0]) # Pa.s
        Re[0] = G * D_h / mu[0]
        
        # array to hold Prandtl's number
        Pr = np.zeros(N_cell)
        Cp = steam_table.Cp_ph(p[0], i[0]) * 1e3    # joules per kg per kelvin
        lambda_l = steam_table.tc_ph(p[0], i[0])
        Pr[0] = Cp * mu[0] / lambda_l
        
        # array to hold Nusselt's number
        Nu = np.zeros(N_cell)
        Nu[0] = 0.023 * Re[0]**0.8 * Pr[0]**0.4
        
        # array to hold heat transfer coefficient
        h = np.zeros(N_cell)
        h[0] = Nu[0] * lambda_l / D_h
        
        T_co[0] = T[0] + q2[0] / h[0]
        
        # cladding
        T_ci[0] = T_co[0] + 1
        for k in range(0, 3):
            T_cAvg = (T_ci[0] + T_co[0]) / 2 - 273.15 # celcius
            lambda_clad = 12.6 + 0.0118 * T_cAvg
            delta_T_c = (q3[0] * r_fo**2 / (2 * lambda_clad)) * np.log(r_co / r_ci)
            T_ci[0] = T_co[0] + delta_T_c
        
        # gas gap    
        lambda_gas_gap = 0.6
        # delta_T_gas_gap = (q1[0] / 4 * np.pi) * (2 * np.log(r_ci / r_fo) / lambda_gas_gap)
        delta_T_gas_gap = (q3[0] * r_fo**2 / (2 * lambda_gas_gap)) * np.log(r_ci / r_fo)
        T_fo[0] = T_ci[0] + delta_T_gas_gap
        
        # fuel
        T_f2[0] = T_fo[0] + 1
        for k in range(0, 3):
            T_f2Avg = (T_f2[0] + T_fo[0]) / 2 # kelvin
            t_f2 = T_f2Avg / 1000
            lambda_f2 = 100 / (7.5408 + 17.692 * t_f2 + 3.6142 * t_f2**2) + 6400 * np.exp(-16.35/t_f2) / t_f2**(5/2)
            delta_T_f2 = (q3[0] * (r_fo**2 - r_f2**2)) / (4 * lambda_f2)
            T_f2[0] = T_fo[0] + delta_T_f2
        
        T_f1[0] = T_f2[0] + 1
        for k in range(0, 3):
            T_f1Avg = (T_f1[0] + T_f2[0]) / 2 # kelvin
            t_f1 = T_f1Avg / 1000
            lambda_f1 = 100 / (7.5408 + 17.692 * t_f1 + 3.6142 * t_f1**2) + 6400 * np.exp(-16.35/t_f1) / t_f1**(5/2)
            delta_T_f1 = (q3[0] * (r_f2**2 - r_f1**2)) / (4 * lambda_f1)
            T_f1[0] = T_f2[0] + delta_T_f1
            
        T_fc[0] = T_fc[0] + 1
        for k in range(0, 3):
            T_fcAvg = (T_fc[0] + T_f1[0]) / 2 # kelvin
            t_fc = T_fcAvg / 1000
            lambda_fc = 100 / (7.5408 + 17.692 * t_fc + 3.6142 * t_fc**2) + 6400 * np.exp(-16.35/t_fc) / t_fc**(5/2)
            delta_T_fc = (q3[0] * r_f1**2) / (4 * lambda_fc)
            T_fc[0] = T_f1[0] + delta_T_fc
        
        # iterate over all the cells
        for j in range(1, N_cell):
            # calculate pressure at the cell
            p[j] = p_inlet - dp[j - 1] / 1e6   # dp in Pa, p in MPa
            
            # cell temperatures
            T[j] = steam_table.t_ph(p[j], i[j])
            
            # q2[j] = q_cells[j] / A_H
            
            # calculate quality at the cell
            i_g = steam_table.hV_p(p[j])
            i_l = steam_table.hL_p(p[j])
            x_eq[j] = (i[j] - i_l) / (i_g - i_l)
            
            if x_eq[j] <= 0:
                rho = steam_table.rho_ph(p[j], i[j]) # kg/m^3
                mu[j] = steam_table.my_ph(p[j], i[j]) # Pa.s
                Re[j] = G * D_h / mu[j]
                C_f = 1 / (- 3.6 * np.log10(((clad_roughness / D_h) / 3.7)**1.11 + 6.9 / Re[j]))**2

                # pressure drops from the inlet to the cell
                dpf[j] = dpf[j - 1] + 4 * C_f * cell_length * G**2 / (D_h * 2 * rho) # Pa
                dpg[j] = dpg[j - 1] + cell_length * rho * 9.81 * 1 # Pa
                dpl[j] = dpl[j - 1] + xi[j] * G**2 / (2 * rho) # Pa
                
                Cp = steam_table.Cp_ph(p[j], i[j]) * 1e3    # joules per kg per kelvin
                lambda_l = steam_table.tc_ph(p[j], i[j])
                Pr[j] = Cp * mu[j] / lambda_l
                Nu[j] = 0.023 * Re[j]**0.8 * Pr[j]**0.4
                h[j] = Nu[j] * lambda_l / D_h
                T_co[j] = T[j] + q2[j] / h[j]
                            
            elif x_eq[j] > 0:
                if x_eq[j] > 1:
                    x_eq[j] = 1
                rho_l = steam_table.rhoL_p(p[j])
                rho_g = steam_table.rhoV_p(p[j])
                T_sat = steam_table.tsat_p(p[j])
                mu_l = steam_table.my_pt(p[j], T_sat - 0.01)
                mu_g = steam_table.my_pt(p[j], T_sat + 0.01)
                rho_mix = rho_l / (x_eq[j] * (rho_l / rho_g - 1) + 1)
                # mu_mix = x_eq[j] * mu_g + (1 - x_eq[j - 1]) * mu_l
                mu_mix = 1 / (x_eq[j] / mu_g + (1 - x_eq[j]) / mu_l)
                mu[j] = mu_mix
                Re[j] = G * D_h / mu[j]
                C_f = 1 / (- 3.6 * np.log10(((clad_roughness / D_h) / 3.7)**1.11 + 6.9 / Re[j]))**2
                
                # pressure drops from the inlet to the cell
                dpf[j] = dpf[j - 1] + 4 * C_f * cell_length * G**2 / (D_h * 2 * rho) # Pa
                dpg[j] = dpg[j - 1] + cell_length * rho_mix * 9.81 * 1 # Pa
                dpl[j] = dpl[j - 1] + xi[j] * G**2 / (2 * rho_mix) # Pa
                # dpa[j] = dpa[j - 1] + G**2 / (rho)  # Pa
                dpa[j] = dpa[j - 1] + G**2 * (1/rho_g - 1/rho_l) * (x_eq[j] - x_eq[j - 1])
                
                # chen correlation
                X_tt = ((1 - x_eq[j]) / x_eq[j])**0.9 * (rho_g / rho_l)**0.5 * (mu_l / mu_g)**0.1
                Re_chen = G * (1 - x_eq[j]) * D_h / mu_l
                if 1/X_tt <= 0.1:
                    F = 1
                else:
                    F = 2.35 * (0.213 + 1/X_tt)**0.736
                S = 1 / (1 + 2.56 * 1e-6 * F**1.463 * Re_chen**1.17)
                lambda_satl = steam_table.tcL_p(p[j])
                Cp_satl = steam_table.CpL_p(p[j]) * 1e3
                Pr_f = Cp_satl * mu_l / lambda_satl
                h_mac = 0.023 * (lambda_satl / D_h) * Re_chen**0.8 * Pr_f**0.4 * F
                
                delta_T_sup = 5
                T_w = T_sat + delta_T_sup
                p_s = steam_table.psat_t(T_w)
                sigma = steam_table.st_p(p[j])
                
                for k in range(0, 3):
                    h_mic = 0.00122 * ((lambda_satl**0.79 * Cp_satl**0.45 * rho_l**0.49) / (sigma**0.5 * mu_l**0.29 * (i_g - i_l)**0.24 * rho_g**0.24)) * delta_T_sup**0.24 * (p_s - p[j])**0.75 * S
                    h[j] = h_mic + h_mac
                    delta_T_sup_new = q2[j] / h[j]
                    delta_T_sup = delta_T_sup_new
                
                delta_T_lb = delta_T_sup
                T_co[j] = T[j] + delta_T_lb
                    
            # cladding temperature calculation
            T_ci[j] = T_co[j] + 1
            for k in range(0, 3):
                T_cAvg = (T_ci[j] + T_co[j]) / 2 - 273.15 # celcius
                lambda_clad = 12.6 + 0.0118 * T_cAvg
                delta_T_c = (q3[j] * r_fo**2 / (2 * lambda_clad)) * np.log(r_co / r_ci)
                T_ci[j] = T_co[j] + delta_T_c
            
            # gas gap temperature calculation
            lambda_gas_gap = 0.6
            delta_T_gas_gap = (q3[j] * r_fo**2 / (2 * lambda_gas_gap)) * np.log(r_ci / r_fo)
            T_fo[j] = T_ci[j] + delta_T_gas_gap
            
            # fuel temperature calculation
            T_f2[j] = T_fo[j] + 1
            for k in range(0, 3):
                T_f2Avg = (T_f2[j] + T_fo[j]) / 2 # kelvin
                t_f2 = T_f2Avg / 1000
                lambda_f2 = 100 / (7.5408 + 17.692 * t_f2 + 3.6142 * t_f2**2) + 6400 * np.exp(-16.35/t_f2) / t_f2**(5/2)
                delta_T_f2 = (q3[j] * (r_fo**2 - r_f2**2)) / (4 * lambda_f2)
                T_f2[j] = T_fo[j] + delta_T_f2
            
            T_f1[j] = T_f2[j] + 1
            for k in range(0, 3):
                T_f1Avg = (T_f1[j] + T_f2[j]) / 2 # kelvin
                t_f1 = T_f1Avg / 1000
                lambda_f1 = 100 / (7.5408 + 17.692 * t_f1 + 3.6142 * t_f1**2) + 6400 * np.exp(-16.35/t_f1) / t_f1**(5/2)
                delta_T_f1 = (q3[j] * (r_f2**2 - r_f1**2)) / (4 * lambda_f1)
                T_f1[j] = T_f2[j] + delta_T_f1
                
            T_fc[j] = T_fc[j] + 1
            for k in range(0, 3):
                T_fcAvg = (T_fc[j] + T_f1[j]) / 2 # kelvin
                t_fc = T_fcAvg / 1000
                lambda_fc = 100 / (7.5408 + 17.692 * t_fc + 3.6142 * t_fc**2) + 6400 * np.exp(-16.35/t_fc) / t_fc**(5/2)
                delta_T_fc = (q3[j] * r_f1**2) / (4 * lambda_fc)
                T_fc[j] = T_f1[j] + delta_T_fc
                        
            # total pressure drop from the inlet to the cell
            dp[j] = dpf[j] + dpg[j] + dpl[j] + dpa[j]
            
            # CHF calculation
            p_R = p[j] / 22 # steam_table.psat_t(T_inlet)
            
            A = a_1 * p_R**a_2 * G_R**(a_3 + a_4 * p_R)
            B = 3.1544e6
            C = c_1 * p_R**c_2 * G_R**(c_3 + c_4 * p_R)
            
            q2_R = q2[j] / B
            
            q2_cr[j] = B * (A - x_eq[0]) / (C + (x_eq[j] - x_eq[0]) / q2_R)
            
            DNBR[j] = q2_cr[j] / q2[j]
            
        total_pressure_drop[q_perc_index, w_perc_index] = dp[N_cell - 1]
        
    x_eq_plot[q_perc_index,] = x_eq
    
hot = i
print(np.amin(DNBR))

# plt.plot(z_cells, hot, label='hot channel')
# plt.plot(z_cells, average, label='average channel')

plt.plot(z_cells, T, label='liquid bulk')
# plt.plot(z_cells, T_co, label='cladding outer')
# plt.plot(z_cells, T_ci, label='cladding inner')
# plt.plot(z_cells, T_fo, label='fuel outer')
# plt.plot(z_cells, T_f2, label='fuel cell 2')
# plt.plot(z_cells, T_f1, label='fuel cell 1')
# plt.plot(z_cells, T_fc, label='fuel center')

# plt.plot(z_cells, q2, label='heat flux')
# plt.plot(z_cells, q2_cr, label='critical heat flux')
# plt.plot(z_cells, DNBR, label='DNBR')
plt.xlabel('Axial height, (m)', fontsize='18')
plt.ylabel('Temperature (kelvin)', fontsize='18')
plt.legend(fontsize='18')
plt.show()

# fig, ax1 = plt.subplots()

# ax1.set_xlabel('Axial distance, z (m)', fontsize='18')
# ax1.set_ylabel('heat flux (W/m^2)', fontsize='18')
# ax1.plot(z_cells, q2, label="heat flux, q''(z)", color='g')
# ax1.plot(z_cells, q2_cr, label="critical heat flux, q''_cr(z)", color='r')
# ax1.tick_params(axis='y')
# ax1.legend(loc='upper left', fontsize='18')

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# ax2.set_ylabel('DNBR', fontsize='18')  # we already handled the x-label with ax1
# ax2.set_ylim(bottom=0.8, top=1.2)
# ax2.plot(z_cells[70:100], DNBR[70:100], label="DNBR", color='b')
# ax2.vlines(z_cells[np.argmin(DNBR)], 0, np.amin(DNBR))
# ax2.hlines(np.amin(DNBR), z_cells[np.argmin(DNBR)], np.amax(z_cells))
# ax2.tick_params(axis='y')
# ax2.legend(loc='upper right', fontsize='18')

# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# plt.show()

# for q_perc_index in range(0, len(power_percs)):
#     plt.plot(flow_percs, total_pressure_drop[q_perc_index, :], label=f'{power_percs[q_perc_index]}% power')

# plt.xlabel('Mass flux percentages (%)', fontsize='18')
# plt.ylabel('Total pressure drop over the channel (Pa)', fontsize='18')
# plt.legend(loc='upper left', fontsize='18')
# # plt.yscale('log')
# plt.show()

# print(np.amax(T_fc), np.amax(T_fo), np.amax(T_ci), np.amax(T_co), np.amax(T))
# print(np.argmax(T_fc), np.argmax(T_fo), np.argmax(T_ci), np.argmax(T_co), np.argmax(T))
# print(z_cells[51], z_cells[53], z_cells[61], z_cells[71], z_cells[79])

# plt.plot([0, r_f1, r_f2, r_fo, r_ci, r_co, 0.0089], [T_fc[51], T_f1[51], T_f2[51], T_fo[51], T_ci[51], T_co[51], T[51]])
# plt.xlabel('radial distance (m)', fontsize='18')
# plt.ylabel('temperature (k)', fontsize='18')
# plt.show()