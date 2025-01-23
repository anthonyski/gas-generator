from powerBalance import *



class systemParameters:
    def __init__(self,loadingCoeff,N_rpm,pressureRise, MR_gasgen, T_o1, MR_chamber, base_ox_mdot, base_fuel_mdot, eta_pump_ox, eta_pump_fuel, eta_turbine, flowCoeff, U):
        self.loadingCoeff = float(loadingCoeff)
        self.N_rpm = float(N_rpm)
        self.pressureRise = float(pressureRise)
        self.MR_gasgen = float(MR_gasgen)
        self.T_o1 = float(T_o1)
        self.MR_chamber = float(MR_chamber)
        self.base_ox_mdot = float(base_ox_mdot)
        self.base_fuel_mdot = float(base_fuel_mdot)
        self.eta_pump_ox = float(eta_pump_ox)
        self.eta_pump_fuel = float(eta_pump_fuel)
        self.eta_turbine = float(eta_turbine)
        self.flowCoeff = float(flowCoeff)
        self.U = float(U)

gasGenSystem = systemParameters(-2, 25000, 300*6894.76,.3,1100,2,3.2,1.6,.6,.6,.8,.4,150)

mdots = balance_optimizer(gasGenSystem)

print("The gas generator mdot is: " + str(mdots))