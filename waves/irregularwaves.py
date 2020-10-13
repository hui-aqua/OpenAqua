import numpy as np
import Airywave as wave
import wave_spectrum as ws



class irregular_sea:
    """
    The default wave spectra is JONSWAP
    ...
    ----------
    *NOTE: The maximum applied simulation time should be less than 3h.
    """
    def __init__(self, significant_wave_height, peak_period, gamma, water_depth, wave_direction):
        """
        Parameters
        ------------
        sigcnificant_wave_height：significant wave height | float |Uinit [m]
        peak_period：peak period | float | Unit [s]
        gamma：gamma | float | Unit [-]
        water_depth： water depth of the sea, assume flat sea floor. A position number | float | Unit [m]
        wave_direction: direction of wave propagation. | float | Unit [degree]
        """
        self.tp=peak_period
        self.hs=significant_wave_height
        self.water_depth=water_depth
        self.list_of_waves=[]

        time_max=3600*2  # 3h, We assume the simulations will not exceed 3h.
        fre_max=3        # we assume the highest eigenfrequency of studied structure is below 3 Hz.
        d_fre=2 * np.pi / time_max                # get the resolution for frequence
        fre_range=np.arange(d_fre,fre_max,d_fre)
        design_wave_spectra=ws.jonswap_spectra(fre_range, significant_wave_height, peak_period, gamma)
        list_xi=np.sqrt(2*d_fre*design_wave_spectra)
        for index, item in enumerate(list_xi):
            wave_period=2*np.pi/fre_range[index]
            self.list_of_waves.append(wave.Airywave(item*2, wave_period, water_depth, wave_direction, np.random.uniform(0,360)))

    def __str__(self):
        """ Print the information of the present object. """
        s0 = 'The environment is irregular waves condition and the specific parameters are:\n'
        s1 = 'significant wave height = ' + str(self.hs) + ' m\n'
        s2 = 'peak period= ' + str(self.tp) + ' s\n'
        s3 = 'Numver of wave components is '+str(len(self.list_of_waves))
        S = s0 + s1 + s2+s3
        return S

    def get_elevations_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for geting the elevations \n
        :return: Get a list of elevations at one position with a time squence \n
        """
        wave_elevations=np.zeros((len(self.list_of_waves),len(time_list)))
        for index, each_wave in enumerate(self.list_of_waves):
            wave_elevations[index]=each_wave.get_elevation(position,time_list)
        return np.sum(wave_elevations,axis=0)

    def get_velocity_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for geting the elevations \n
        :return: Get a list of elevations at one position with a time squence \n
        """
        waves_velocities=np.zeros((len(self.list_of_waves),len(time_list),3))
        for index, each_wave in enumerate(self.list_of_waves):
            waves_velocities[index]=each_wave.get_velocity_with_time(position,time_list)
        return np.sum(waves_velocities,axis=0)

    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        wave_elevations=np.zeros((len(self.list_of_waves),len(list_of_point)))
        for index, each_wave in enumerate(self.list_of_waves):
            wave_elevations[index]=each_wave.get_elevation_at_nodes(list_of_point,global_time)
        return np.sum(wave_elevations,axis=0)

    def get_velocity_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_velocity=np.zeros((len(self.list_of_waves),len(list_of_point),3))
        list_elevation=self.get_elevation_at_nodes(list_of_point, global_time)
        # wheeler streching method
        list_of_point[:,2] = (list_of_point[:, 2] - list_elevation)*self.water_depth / (self.water_depth + list_elevation)

        for index, each_wave in enumerate(self.list_of_waves):
            node_velocity[index]=each_wave.get_velocity_at_nodes(list_of_point,global_time)
        velo=np.sum(node_velocity,axis=0)
        # ensure the velocity is zero above the wave elevation
        for index in range(len(list_elevation)):
            if list_of_point[index,2]>list_elevation[index]:
                velo[index]=0
        return velo

    def get_acceleration_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_acceleration=np.zeros((len(self.list_of_waves),len(list_of_point),3))
        list_elevation=self.get_elevation_at_nodes(list_of_point, global_time)
        # wheeler streching method
        list_of_point[:,2] = (list_of_point[:, 2] - list_elevation)*self.water_depth / (self.water_depth + list_elevation)

        for index, each_wave in enumerate(self.list_of_waves):
            node_acceleration[index]=each_wave.get_velocity_at_nodes(list_of_point,global_time)
        acce=np.sum(node_acceleration,axis=0)
        for index in range(len(list_elevation)):
            if list_of_point[index,2]>list_elevation[index]:
                acce[index]=0
        return acce

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    plt.rcParams['image.cmap'] = 'summer'

    sea_state=irregular_sea(5,8,3,60,0)
    time_max = 3600 # [s]
    dt=1
    time_frame=np.arange(0, time_max, dt)
    print(sea_state)


    g1 = 2
    g2 = 1
    gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots
    elevations=sea_state.get_elevations_with_time(np.array([0,0,0]),time_frame)
    position=np.zeros((600,3))
    for i in range(20):
        for j in range(30):
            position[i+20*j]=[5*i,0,5.0-j]
    position2=np.zeros((100,3))
    for i in range(100):
        position2[i]=[i,0,0]
    plt.figure(figsize=(6.7, 4.5))
    # for t in time_frame:
    t=0
    plt.clf()
    ax = plt.subplot(gs[0, 0])
    plt.title("surface elevation with time at position x=0,y=0")
    plt.plot(time_frame,elevations ,color='k',linewidth=0.5)
    plt.xlabel("time (s)")
    plt.ylabel("surface elevation (m)")
    plt.xlim(0, 3600)
    plt.ylim(-5, 5)
    ax = plt.subplot(gs[1, 0])
    plt.title("serface elevation and velocity at "+str(t)+" s")
    velocity=sea_state.get_velocity_at_nodes(position,t)
    plt.text(70,-35,"Maximum velocity "+str(round(max(np.linalg.norm(velocity,axis=1)),2))+" m/s")
    plt.quiver(position[:,0],
                position[:,2],
                velocity[:,0],
                velocity[:,2],
                units='xy',
                scale=0.5,
                )
    plt.plot([i for i in range(100)],sea_state.get_elevation_at_nodes(position2,t))
    plt.plot([0,100],[0,0],color='r',linewidth=0.5)
    plt.xlabel("X (m)")
    plt.ylabel("elevation (m)")
    plt.xlim(0, 100)
    plt.ylim(-26, 6)
    # q.set_array(np.linspace(0,2,10))
    # fig.colorbar(q)
    plt.tight_layout()
    plt.savefig('./figures/irrugularwaves_'+str(t)+'.png', dpi=600)


    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # # time_frame=np.arange(1000,4600,dt)
    # # yita=sea_state.get_elevations_with_time([0,0,0],time_frame)
    # # print(np.std(yita))
    # # plt.plot(time_frame,yita)
    # x_axis=np.array(np.arange(-10,10,1.1)).tolist()
    # position=np.zeros((10*len(x_axis),3))
    # for i in range(10):
    #     for index, item in enumerate(x_axis):
    #         position[i*len(x_axis)+index]=[6*i,6*i,-float(item)/2]
    # # print(position.shape)
    # # print(position)
    # velocity=sea_state.get_velocity_at_nodes(position, 0)
    # velocity_mag=[]
    # for each in velocity:
    #     velocity_mag.append(np.linalg.norm(each))
    # print("max velocity is " +str(max(velocity_mag)) + " m/s")
    # print("and the position is "+str(position[velocity_mag.index(max(velocity_mag))]))
    # # Flatten and normalize
    # velocity_mag=np.array(velocity_mag)
    # normal_velo = (velocity_mag.ravel() - velocity_mag.min()) / velocity_mag.ptp()
    # # print(normal_velo)
    # # Repeat for each body line and two head lines
    # c = np.concatenate((normal_velo, np.repeat(normal_velo, 2)))
    # # Colormap
    # c = plt.cm.summer(c)
    # q=ax.quiver(position[:,0],
    #             position[:,1],
    #             position[:,2],
    #             velocity[:,0],
    #             velocity[:,1],
    #             velocity[:,2],
    #             colors=c,
    #             normalize=False,
    #             )
    # # q.set_array(np.linspace(0,2,10))
    # fig.colorbar(q)

    # for i in range(60):
    #     position=np.zeros((len(x_axis),3))
    #     for index, item in enumerate(x_axis):
    #         position[index]=[item,i,0]
    #     yita=sea_state.get_elevation_at_nodes(position, 0)
    #     ax.plot(x_axis,[i]*len(x_axis),yita,color="b")

    # # plt.plot(x_axis,yita)
    # ax.set_title("JONSWAP sea condition Hs=4m, Tp=8s r=3, t="+str(0)+"s")
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # ax.set_xlim(0, 60)
    # ax.set_ylim(0, 60)
    # ax.set_zlim(-30, 5)
    # plt.savefig('./figures/waves.png', dpi=300)
    plt.show()