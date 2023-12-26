import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

"""Script to generate plots for batch kinetics in a few conditions.
Considers well-mixed, batch reaction conditions, with a constant surface area and volume.
Assumes all phases are solid or aqueous, and there is no head space.

Runs as script.
 """


def precipitationRateCalculation(ionActivityProduct, equilibriumK, rateConstantk, reactionOrdern=1):
    """Calculates the reaction rate for a rate law given as:
    r = k*(1-IAP/K)^n
    rate, r: gives the surface area normalized reaction rate for the reactant in mol/s
    rateConstantk: k is the rate constant times the (implied) surface area, should be in mol/s
    ionActivityProduct, IAP: is the ion activity product
    reactionOrdern, n: The reaction order, which is assumed 1 default
    """
    rate = rateConstantk*(1-ionActivityProduct/equilibriumK)**reactionOrdern
    return rate


def simulatePrecipNoCO3(totCa, totC, k_react, Keq_solid, n_rxnOrder, dt, t_end):
    # Initialize parameters and calculation arrays
    times = np.arange(0, t_end+dt, dt)
    NPoints = np.shape(times)
    Ca = np.zeros(NPoints)
    Ca[0] = totCa
    CO3 = np.zeros(NPoints)
    CO3[0] = totC
    rates = np.zeros(NPoints)

    # Begin calculation
    for timeIndex in range(0, NPoints[0]-1, 1):
        IAP = Ca[timeIndex]*CO3[timeIndex]
        rate = precipitationRateCalculation(IAP, Keq_solid, k_react, n_rxnOrder)
        rates[timeIndex] = rate
        Ca[timeIndex+1] = Ca[timeIndex]+rate*dt
        CO3[timeIndex+1] = CO3[timeIndex]+rate*dt
    rates[-1] = precipitationRateCalculation(Ca[-1]*CO3[-1], Keq_solid, k_react, n_rxnOrder)
    plt.figure(1)
    plt.plot(times, Ca, label='Calcium, no equilibrium calc')
    plt.plot(times, CO3, label='Carbonate, no equilibrium calc')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (M)')

    plt.figure(2)
    plt.plot(times, rates)
    # Save results to dictionary
    resultDict = {'CO3': CO3, 'Ca': Ca, 'rates': rates, 'times': times}
    return resultDict


def HRootFinder(K1,K2,Ca_total,C_total):
        #TODO: TEST ME
        Kw = 10**-14
        functionH = lambda x: K1*K2*x**4+(K2+2*Ca_total*K2*K1)*x**3+(1+2*Ca_total-Kw*K1*K2-K2*C_tot)*x**2+(2*Ca_total-2*C_total-Kw*K2)*x-Kw
        roots = root(functionH, [-1, -10**-7, 10**-7, 1])
        return roots


def simulatePrecipCO3Equil(totCa, totC, k_react, Keq_solid, n_rxnOrder, K1, K2, dt, t_end):
    # Initialize parameters
    K2 = 10**-14
    times = np.arange(0, t_end+dt, dt)
    NPoints = np.shape(times)
    Ca = np.zeros(NPoints)
    CO3 = np.zeros(NPoints)
    C_total = np.zeros(NPoints)
    C_total[0] = totC
    Ca_total = np.zeros(NPoints)
    Ca_total[0] = totCa
    H = np.zeros(NPoints)
    rates = np.zeros(NPoints)

    # Begin calculation
    for timeIndex in range(0, NPoints[0]-1, 1):
        # Calculate equlibrium state, starting with H+
        # Need funtion to actually calculate roots and then give the physically meaningful result
        time = timeIndex
    result = {'totC': C_total, 'totCa': Ca_total, 'H+': H, 'Ca': Ca, 'CO3': CO3, 'rates':rates}
    return result


def main():
    # User inputs/flags
    Ca_Initial = 1E-3 # Molar
    CO3_Initial = 1E-3 # Molar
    k = 1E-6 # Mol/s
    Keq = 10**-8.54
    dt = 1E-3 # s
    t_end = 100
    reactionOrder = 1

    nonEquilData = simulatePrecipNoCO3(Ca_Initial, CO3_Initial, k, Keq, reactionOrder, dt, t_end)

    return 0


main()