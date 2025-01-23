import numpy as np


def solve_U(cz, flow_coeff):
    return cz/flow_coeff


def solve_beta(cz, U, alpha1):
    return np.arctan((cz*np.tan(alpha1*np.pi/180)-U)/cz)*180/np.pi

def solve_T1(To1, cz, cp, alpha1):
    return To1-cz**2/2/cp/(np.cos(alpha1*np.pi/180)**2)


def cosd(angle):
    return np.cos(angle*np.pi/180)

def tand(angle):
    return np.tan(angle*np.pi/180)

def sind(angle):
    return np.sin(angle*np.pi/180)

def carters_rule(camber, solidity):
    return .25*np.abs(camber)/np.sqrt(solidity)

def solve_ideal_gas_density(p,T, R):
    return p/T/R

def solve_static_rho(gamma, M, rho_o):
    return rho_o*((1+(gamma-1)/2*M**2)**(-1/(gamma-1)))

def solve_area(mdot, cz, rho):
    return mdot/rho/cz

def solve_rhub(rm, area):
    return np.sqrt(rm**2-area/2/np.pi)

def solve_rtip(rm, rh):
    return np.sqrt(2*rm**2-rh**2)

def solve_rm(rt, rh):
    return np.sqrt((rt**2-rh**2)/2)

def solve_pressure_isentropic(po, gamma, M):
    exponent = -gamma/(gamma-1)
    return po*((1+(gamma-1)/2*(M**2))**(exponent))

def solve_rm(rt, rh):
    return np.sqrt((rt**2+rh**2)/2)

def solve_eta_stage(gamma, Pr, U, loading_coeff, cp, To1):
    exponent = (gamma-1)/gamma + 0j
    return (Pr**exponent-1)/((U**2*loading_coeff)/(cp*To1))

def solve_pressure_coeff(angle1, angle2, pressure_loss_coeff):
    return 1-cosd(angle1)/cosd(angle2)-pressure_loss_coeff

def solve_To_relation(pressure_ratio, To2, gamma):
    exponent = gamma/(gamma-1) + 0j
    return To2*pressure_ratio**exponent

def solve_axial_solidity(angle_in, angle_out, c_z, M_out, gamma):
    exponent = gamma/(gamma-1)
    a = (c_z*tand(angle_in/(c_z*tand(angle_out)))-1)
    b = ((gamma*M_out**2/2)/(1+(gamma-1)/2*M_out**2)**exponent-1)
    return a*sind(2*angle_out)*b/.8

def sutherlandsLaw(mu0, T, T_normal):
    return float(mu0*(T/T_normal)**0.7)