import numpy as np
import matplotlib as mplot

import scipy
import sympy as sp

from utils import *





def solveTurbine(turbine):
    gamma = 1.1253
    MW = 11.498
    C_p = 3893    
    c_z = turbine.flowCoeff*turbine.U ##assume change in c_z is low across stator (since change in density is low 1-2%)
    rho_o1 = solve_ideal_gas_density(turbine.p_o1, turbine.T_o1, turbine.R)
    rho_1 = solve_static_rho(gamma, turbine.M_in, rho_o1)
    A_zN = solve_area(turbine.mDot, c_z, rho_1)
    r_h = solve_rhub(turbine.r_m, A_zN)
    r_tip = solve_rtip(turbine.r_m, r_h)
    h_N = r_tip-r_h
    T_2 = float(turbine.T_o1*(1+(gamma-1)/2*turbine.M_2**2)**-1)
    
    alpha_2 = 72 ##general upper bound for alpha2
    c_theta_2 = c_z*tand(alpha_2)
    c_2 = np.sqrt(c_z**2+c_theta_2**2)
    sigma_zN = solve_axial_solidity(0, alpha_2, c_z, turbine.M_2, gamma)
    sigma_N = sigma_zN/np.cos(np.atan(c_theta_2/2/c_z))
    nu = 0.00003621 ##millipoise to pa*s from CEA
    rho_2_guess = rho_1*.998
    mu0 = rho_2_guess*nu
    # mu_N = float((sutherlandsLaw(mu0, T_2, T_o1)))
    mu_N = nu*rho_2_guess
    nu_N = mu_N/rho_2_guess
    # o_N = Re_t_N*mu_N/c_2/rho_2_guess
    o_N = (7.19e-7)*turbine.Re_t_N/c_2/rho_2_guess
    s_N = o_N/cosd(alpha_2)
    b_zN = sigma_zN*s_N
    b_N = sigma_N*s_N
    D_hN =sigma_N*s_N*h_N*cosd(alpha_2)/(s_N*cosd(alpha_2)+h_N)
    Re_e_N = c_2*D_hN*rho_2_guess/mu_N
    zeta_star = 1.04+0.06*(alpha_2/100)**2
    C_jN = 0.993+0.021*(b_zN/h_N)
    zeta_N = (zeta_star*C_jN-1)*(10**5/Re_e_N)**0.25
    exponent = gamma/(gamma-1)
    c_2 = np.sqrt(c_z**2+c_theta_2**2)
    debug_p_o1 = turbine.p_o1
    p_o2 = -(1-((1-c_2**2/2/C_p/turbine.T_o1/(1-zeta_N))/(1-c_2**2/2/C_p/turbine.T_o1))**exponent)*turbine.p_o1+turbine.p_o1
    pressure_loss_ratio_nozzle = p_o2/turbine.p_o1
    print("The blade height is: " + str(h_N)) 
    T_o3 = abs(np.divide(turbine.loadingCoeff*turbine.U**2,C_p)+turbine.T_o1)
    rho_2o_real = solve_ideal_gas_density(p_o2, turbine.T_o1, turbine.R)
    rho_2_real = solve_static_rho(gamma, turbine.M_2, rho_2o_real)
    A_zR = solve_area(turbine.mDot,c_z,rho_2_real)
    beta_2 = np.atan((c_z*tand(alpha_2)-turbine.U)/c_z)*180/np.pi 
    beta_3 = -beta_2 ## full impulse turbine
    alpha_3 = np.atan((2+turbine.loadingCoeff)/2/turbine.flowCoeff)*180/np.pi
    w_theta_2 = c_z*tand(beta_2)
    w_theta_3 = c_z*tand(beta_3)
    T_3 = T_o3-c_z**2/2/C_p/cosd(alpha_3)
    M_3 = c_z/cosd(alpha_3)/np.sqrt(gamma*turbine.R*T_3)
    M_3rel = c_z/cosd(beta_3)/np.sqrt(gamma*turbine.R*T_3)
    M_2rel = c_z/cosd(beta_2)/np.sqrt(gamma*turbine.R*T_2)
    rho_3_guess = solve_static_rho(gamma, M_3, rho_2_real)
    sigma_zR = solve_axial_solidity(beta_2, beta_3, c_z, M_3rel, gamma)
    sigma_R = sigma_zR/np.cos(np.atan((w_theta_2+w_theta_3)/2/c_z))
    w_3 = np.sqrt(w_theta_3**2+c_z**2)
    # mu_R = sutherlandsLaw(mu0, T_3, T_o1)
    # nu_R = mu_N/rho_3_guess
    o_R = turbine.Re_t_R*7.19e-7/w_3/rho_3_guess
    s_R = o_R/cosd(alpha_3)
    b_zR = sigma_zR*s_R
    b_R = sigma_R*s_R
    h_R = h_N
    D_hR =sigma_R*s_R*h_R*cosd(beta_3)/(s_R*cosd(beta_3)+h_R)
    T_o2_rel = T_2*(1+(gamma-1)/2*M_2rel**2)
    Re_e_R = c_2*D_hN*rho_3_guess/mu0
    zeta_star = 1.04+0.06*(beta_3/100)**2   
    C_jR = 0.993+0.021*(b_zR/h_R)
    zeta_R = (zeta_star*C_jR-1)*(10**5/Re_e_R)**0.25
    
    p_2 = p_o2*(1+(gamma-1)/2*turbine.M_2**2)**(-exponent)
    p_o2_rel = p_2*(1+(gamma-1)/2*M_2rel**2)**exponent
    p_o3_rel = np.multiply(-(1-((1-w_3**2/2/C_p/T_o2_rel*(1/(1-zeta_R)))/(1-(w_3**2)/2/C_p/T_o2_rel))**exponent), p_o2_rel)+p_o2_rel
    p_3 = p_o3_rel*(1+(gamma-1)/2*M_3rel**2)**(-exponent)
    p_o3 = p_3*(1+(gamma-1)/2*M_3**2)**exponent


    
    numStatorBlades = 2*np.pi*turbine.r_m/s_N
    print("The stator has: " + str(numStatorBlades) + " blades")
   
   
    numRotorBlades = 2*np.pi*turbine.r_m/s_R
    print("The rotor has: " + str(numRotorBlades) + " blades")
  

    return (1-T_o3/turbine.T_o1)/(1-p_o3/turbine.p_o1)


def balance_optimizer(system):
    gamma = 1.1253
    MW = 11.498
    C_p = 3893
    
    ##constants

    pressureRise = 300*6894.76 #psi to pa
    MR_gasgen = .3
    T_o1 = 1100
    MR_chamber = 2

    gasGenAverageRho = (800+1100*MR_gasgen)/(1+MR_gasgen)
    
    rotor_Re = 1e5
    stator_Re = 1e5
    rhoOx = 1100
    rhoFuel = 800

    # a = -loadingCoeff*(U**2)*eta_turbine-pressureRise/gasGenAverageRho/eta_pump_fuel
    
    
    # b = Base_ox_mdot/rhoOx*pressureRise/eta_pump_ox+Base_fuel_mdot/rhoFuel*pressureRise/eta_pump_fuel
    
    # mDot_gasgen = b/a
    
    
    def balance(g):
        return system.base_ox_mdot/rhoOx*pressureRise/system.eta_pump_ox+system.base_fuel_mdot/rhoFuel*pressureRise/system.eta_pump_fuel+g/gasGenAverageRho/system.eta_pump_fuel*pressureRise+system.U**2*g*system.loadingCoeff
    
        
    mDot_gasgen = scipy.optimize.newton(balance, .1)
    pumpPower = system.base_ox_mdot/rhoOx*pressureRise/system.eta_pump_ox+system.base_fuel_mdot/rhoFuel*pressureRise/system.eta_pump_fuel+mDot_gasgen/gasGenAverageRho/system.eta_pump_fuel
    r_m = system.U*30/np.pi/system.N_rpm ##answer in meters
    print("The turbine mean radius is: "+ str(r_m*39.37))
    turbine = turbineParameters(mDot_gasgen, T_o1, pressureRise*0.8, gamma, system.loadingCoeff, r_m, system.U, system.N_rpm, 8314.5/MW, system.flowCoeff, .3, 1, stator_Re, rotor_Re, pumpPower*1.05)
    efficiency = solveTurbine(turbine)
    print("The turbine efficiency is: " + str(efficiency))

    return mDot_gasgen


class turbineParameters:
    def __init__(self,mDot,T_o1, p_o1, gamma, loadingCoeff, r_m, U_r, N_rpm, R, flowCoeff, M_in, M_2, Re_t_N, Re_t_R, Power):
        self.mDot = float(mDot)
        self.T_o1 = float(T_o1)
        self.p_o1 = float(p_o1)
        self.gamma = float(gamma)
        self.loadingCoeff = float(loadingCoeff)
        self.r_m = float(r_m)
        self.U = float(U_r)
        self.N_rpm = float(N_rpm)
        self.R = float(R)
        self.flowCoeff = float(flowCoeff)
        self.M_in = float(M_in)
        self.M_2 = float(M_2)
        self.Re_t_N = float(Re_t_N)
        self.Re_t_R = float(Re_t_R)
        self.Power = float(Power)











    

    
    




    
     


    


    


    


    
