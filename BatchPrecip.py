import numpy as np
import matplotlib.pyplot as plt

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


def main():
    # User inputs/flags
    Ca_Initial = 1E-3 # Molar
    CO3_Initial = 1E-3 # Molar
    k = 1E-6 # Mol/s
    Keq = 10**-8.54
    dt = 1E-3 # s
    t_end = 100
    k = 1E-6 # Mol/s
    Keq = 10**-8.54
    dt = 1E-3 # s
    reactionOrder = 1

    # Initialize parameters and calculation arrays
    times = np.arange(0, t_end+dt, dt)
    NPoints = np.shape(times)
    Ca = np.zeros(NPoints)
    Ca[0] = Ca_Initial
    CO3 = np.zeros(NPoints)
    CO3[0] = CO3_Initial
    rates = np.zeros(NPoints)

    # Begin calculation
    for timeIndex in range(0, NPoints[0]-1, 1):
        IAP = Ca[timeIndex]*CO3[timeIndex]
        rate = precipitationRateCalculation(IAP, Keq, k, reactionOrder)
        rates[timeIndex] = rate
        Ca[timeIndex+1] = Ca[timeIndex]+rate*dt
        CO3[timeIndex+1] = CO3[timeIndex]+rate*dt
    rates[-1] = precipitationRateCalculation(Ca[-1]*CO3[-1], Keq, k, reactionOrder)
    plt.figure(1)
    plt.plot(times, Ca, label='Calcium')
    plt.plot(times, CO3, label='Carbonate')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (M)')

    plt.figure(2)
    plt.plot(times, rates)
    plt.show()

    return 0


main()