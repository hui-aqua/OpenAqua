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
import Airywave as aw

class Stocks2wave(aw.Airywave):
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


if __name__ == "__main__":
    # validation of code, ref to Figure3-3 on page 47 from DNV GL-RP205 Ver. 2008
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    g1=1
    g2=2
    gs = gridspec.GridSpec(g1, g2)           # Create 1x2 sub plots
    
    water_d=[10,20,30,40,50,60,80,100,1000]
    waves_length=[]
    waves_phasevelocity=[]
    wave_period=np.linspace(1,20,50)
    for item in water_d:
        waves_length.append([Stocks2wave(1,i,item,0).wave_Length for i in wave_period])
        waves_phasevelocity.append([Stocks2wave(1,i,item,0).wave_phase_velocity for i in wave_period])
    
    plt.figure(figsize=(6.3, 4.0))
    
    ax = plt.subplot(gs[0, 0])
    for item in waves_length:
        plt.plot(wave_period,item,label= "Depth "+str(water_d[waves_length.index(item)]))
    plt.xlabel("Wave period (s)")
    plt.ylabel("Wave length (m)")
    plt.xlim(0, 20)
    plt.ylim(0,700)
    plt.grid(True)
    plt.legend()
            
    ax = plt.subplot(gs[0, 1])
    for item in waves_phasevelocity:
        plt.plot(wave_period,item,label= "Depth "+str(water_d[waves_phasevelocity.index(item)]))
    plt.xlabel("Wave period (s)")
    plt.ylabel("Phase velocity (m/s)")
    plt.xlim(0, 20)
    plt.ylim(0,35)
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()
        
        
        
    # # validation 2 shows the wave elevation according to time and space
    # import matplotlib.pyplot as plt
    # import matplotlib.gridspec as gridspec
    # g1=2
    # g2=1
    # gs = gridspec.GridSpec(g1, g2)           # Create 1x2 sub plots
    
    # water_d=[10,20,30,40,50,60,80,100,1000]
    
    # time_slice=np.linspace(0,100,1000)
    # wave_elevation_with_time=[]
    
    # space_slice=np.ones((1000,3))
    # x_axis=[]
    # for posi in range(1000):
    #     space_slice[posi]=[posi/10,0,0]
    #     x_axis.append(posi/10)
    # wave_elevation_with_x=[]
    # wave_height=1.5
    # wave_period=6
   
    # for item in water_d:
    #     wave_elevation_with_time.append([Stocks2wave(wave_height,wave_period,item,0).get_elevation([0,0,0],i) for i in time_slice])
    #     wave_elevation_with_x.append(Stocks2wave(wave_height,wave_period,item,0).get_elevation_at_nodes(space_slice,0))
    
    # plt.figure(figsize=(6.3, 4.0))
    
    # ax = plt.subplot(gs[0, 0])
    # for item in water_d:
    #     plt.plot(time_slice,wave_elevation_with_time[water_d.index(item)],label= "Depth "+str(item))
    # plt.xlabel("Time (s)")
    # plt.ylabel("Wave elevation (m)")
    # plt.xlim(0, 100)
    # plt.ylim(-3,3)
    # plt.grid(True)
    # plt.legend()
            
    # ax = plt.subplot(gs[1,0])
    # for item in water_d:
    #     plt.plot(x_axis,wave_elevation_with_x[water_d.index(item)],label= "Depth "+str(item))
    # plt.xlabel("X (m)")
    # plt.ylabel("Wave elevation (m)")
    # plt.xlim(0, 100)
    # plt.ylim(-3,3)
    # plt.grid(True)
    # plt.legend()
    
    # plt.tight_layout()
    # plt.show()
    
    
    
    
    