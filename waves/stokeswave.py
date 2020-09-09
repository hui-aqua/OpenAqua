"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
"""
import numpy as np
from numpy import pi


class Stocks2wave():
    def __str__(self):
        s0 = 'The environment is stokes second order wave theory condition and the specific parameters are:\n'
        s1 = 'water Depth= ' + str(self.waterDepth) + ' m\n'
        s2 = 'wave Period= ' + str(self.wavePeriod) + ' s\n'
        s3 = 'wave Length= ' + str(self.waveLength) + ' m\n'
        s4 = 'wave Height= ' + str(self.waveHeight) + ' m\n'
        S = s0 + s1 + s2 + s3 + s4
        return S

    def get_velocity(self, posi, time):
        zeta = self.wavek * posi[0] - 2 * pi / self.wavePeriod * time
        yita = self.get_surface(posi, time)
        if posi[2] < yita:
            horizonvelocity = self.piht * np.cosh(self.wavek * (posi[2] + self.waterDepth)) * np.cos(zeta) / np.sinh(
                self.wavek * self.waterDepth) + 0.75 * self.piht * self.pihl * np.cosh(
                2 * self.wavek * (posi[2] + self.waterDepth)) * np.cos(2 * zeta) / pow(4, np.sinh(
                self.wavek * self.waterDepth))
            vericalvelocity = self.piht * np.sinh(self.wavek * (posi[2] + self.waterDepth)) * np.sin(zeta) / np.sinh(
                self.wavek * self.waterDepth) + 0.75 * self.piht * self.pihl * np.sinh(
                2 * self.wavek * (posi[2] + self.waterDepth)) * np.sin(2 * zeta) / pow(4, np.sinh(
                self.wavek * self.waterDepth))
            if self.waterDepth > self.waveLength * 0.5:  # if it is deep water
                horizonvelocity = self.piht * np.exp(self.wavek * posi[2]) * np.cos(
                    zeta) + 0.75 * self.piht * self.pihl * np.cosh(
                    2 * self.wavek * (posi[2] + self.waterDepth)) * np.cos(2 * zeta) / pow(4, np.sinh(
                    self.wavek * self.waterDepth))
                vericalvelocity = self.piht * np.exp(self.wavek * posi[2]) * np.sin(
                    zeta) + 0.75 * self.piht * self.pihl * np.sinh(
                    2 * self.wavek * (posi[2] + self.waterDepth)) * np.sin(2 * zeta) / pow(4, np.sinh(
                    self.wavek * self.waterDepth))
        else:
            horizonvelocity = 0.0
            vericalvelocity = 0.0
        velo = np.array([0.0, 0.0, 0.0])
        velo[0] = horizonvelocity
        velo[1] = 0.0
        velo[2] = vericalvelocity
        return velo

    def get_acceleration(self, posi, time):
        zeta = self.wavek * posi[0] - 2 * pi / self.wavePeriod * time
        yita = self.get_surface(posi, time)
        if posi[2] < yita:
            horizontalacceleration = self.piht2 * np.cosh(self.wavek * (posi[2] + self.waterDepth)) * np.sin(
                zeta) / np.sinh(self.wavek * self.waterDepth) + 1.5 * self.piht2 * self.pihl * np.cosh(
                2 * self.wavek * (posi[2] + self.waterDepth)) * np.sin(2 * zeta) / pow(4, np.sinh(
                self.wavek * self.waterDepth))
            vericalaccelateration = -self.piht2 * np.sinh(self.wavek * (posi[2] + self.waterDepth)) * np.cos(
                zeta) / np.sinh(self.wavek * self.waterDepth) - 1.5 * self.piht2 * self.pihl * np.sinh(
                2 * self.wavek * (posi[2] + self.waterDepth)) * np.cos(2 * zeta) / pow(4, np.sinh(
                self.wavek * self.waterDepth))
            if self.waterDepth > self.waveLength * 0.5:  # if it is deep water
                horizontalacceleration = self.piht2 * np.exp(self.wavek * posi[2]) * np.sin(
                    zeta) + 1.5 * self.piht2 * self.pihl * np.cosh(
                    2 * self.wavek * (posi[2] + self.waterDepth)) * np.sin(2 * zeta) / pow(4, np.sinh(
                    self.wavek * self.waterDepth))
                vericalaccelateration = -self.piht2 * np.exp(self.wavek * posi[2]) * np.cos(
                    zeta) - 1.5 * self.piht2 * self.pihl * np.sinh(
                    2 * self.wavek * (posi[2] + self.waterDepth)) * np.cos(2 * zeta) / pow(4, np.sinh(
                    self.wavek * self.waterDepth))
        else:
            horizontalacceleration = 0.0
            vericalaccelateration = 0.0
        acce = np.array([0.0, 0.0, 0.0])
        acce[0] = horizontalacceleration
        acce[1] = 0.0
        acce[2] = vericalaccelateration
        return acce

    def get_surface(self, posi, time):
        zeta = self.wavek * posi[0] - 2 * pi / self.wavePeriod * time
        yita = self.waveHeight / 2 * np.cos(
            zeta) + pi * self.waveHeight * self.waveHeight / 8.0 / self.waveLength * np.cosh(
            self.wavek * self.waterDepth) * (2 + np.cosh(2 * self.wavek * self.waterDepth)) * np.cos(2 * zeta) / pow(
            np.sinh(self.wavek * self.waterDepth), 3)
        return yita
