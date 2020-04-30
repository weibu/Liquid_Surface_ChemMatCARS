from periodictable import Sm, Gd
class NeutronScatteringFactors(object):
    def __init__(self, isotope, Pr, Pi, Eo, E_lambda, gamma_lambda, Gamma_lambda, Delta_E = 0):
         # Paper gives coefficients Ar, Br, Cr, ... from lowest to highest
         self.isotope = isotope
         self.Pr = Pr[::-1]
         self.Pi = Pi[::-1]
         self.Eo = Eo
         self.E_lambda = E_lambda
         self.gamma_lambda = gamma_lambda
         self.Gamma_lambda = Gamma_lambda
         self.Delta_E = Delta_E
         self.a = 1.16 * isotope.isotope**(1/3.) + 0.6 # fm
         self.a *= 0.1 # 1e-12 cm

    def sf(self, E): # eV
        # Eq. 1 of Lynn
        k = 0.0021968*numpy.sqrt(E) # 1e12/cm

        # Eq. 11, 12 of Lynn, with coefficients reversed
        Rex = numpy.polyval(self.Pr,E-Eo)
        Sex = numpy.polyval(self.Pi,E-Eo)

        # Eq. 10 of Lynn, with common subexpressions extracted for readability
        Ep = self.E_lambda + Delta_lambda - E
        Gp = self.Gamma_lambda/2
        gp = self.a*self.gamma_lambda**2/(Ep + Gp**2)
        Rp = self.a*(1-Rex)
        Sp = self.a*Sex
        real = -Rp + 2*k*Rp*Sp + gp*(Ep + 2*k*Gp*Rp)
        imag = Sp + k*(Rp**2 - Sp**2) + gp*(Gp - 2*k*Ep*Rp)
        return real,imag

Sm149 = NeutronScatteringFactors(Sm[149],
            Pr=[0.0,0.75,0.441], Pi=[0.013, 0.0, 0.025],
            Eo=0.025, E_lambda=0.0973, 
            gamma_lambda =0.573, Gamma_lambda = 0.0656)

