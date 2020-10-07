import numpy as np
import Airywave as wave
import wave_spectrum as ws



class irregular_sea:
    """
    The default wave spectra is JONSWAP 
    """

    def __init__(self, significant_wave_height, peak_period, gamma, water_depth, wave_direction):    
        self.tp=peak_period
        self.hs=significant_wave_height
        self.list_of_waves=[]

        time_max=3600*3  # 3h, We assume the simulations will not exceed 3h.
        fre_max=3        # we assume the highest eigenfrequency of studied structure is below 3 Hz.
        
        d_fre=2 * np.pi / time_max                # get the resolution for frequence
        fre_range=np.arange(d_fre,fre_max,d_fre)  
        design_wave_spectra=ws.jonswap_spectra(fre_range, significant_wave_height, peak_period, gamma)
        list_xi=np.sqrt(2*d_fre*design_wave_spectra)
        for index, item in enumerate(list_xi):
            wave_period=2*np.pi/fre_range[index]
            self.list_of_waves.append(wave.Airywave(item*2, wave_period, water_depth, wave_direction, np.random.uniform(0,180)))
        

    def __str__(self):
        s0 = 'The environment is irregular waves condition and the specific parameters are:\n'
        s1 = 'significant wave height = ' + str(self.hs) + ' m\n'
        s2 = 'peak period= ' + str(self.tp) + ' s\n'
        S = s0 + s1 + s2 
        return S       

    def get_elevations_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for geting the elevations \n
        :return: Get a list of elevations at one position with a time squence \n
        """
        wave_elevations=np.zeros((len(self.list_of_waves),len(time_list)))
        for index, wach_wave in enumerate(self.list_of_waves):
            wave_elevations[index]=wach_wave.get_elevation(position,time_list)
        return np.sum(wave_elevations,axis=0)
    
    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        wave_elevations=np.zeros((len(self.list_of_waves),len(list_of_point)))
        for index, wach_wave in enumerate(self.list_of_waves):
            wave_elevations[index]=wach_wave.get_elevation_at_nodes(list_of_point,global_time)
        return np.sum(wave_elevations,axis=0)

    def get_velocity_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_velocity=np.zeros((len(self.list_of_waves),len(list_of_point),3))
        for index, wach_wave in enumerate(self.list_of_waves):
            node_velocity[index]=wach_wave.get_velocity_at_nodes(list_of_point,global_time)
        # print(np.sum(node_velocity,axis=0))
        # print(np.sum(node_velocity,axis=0).shape)
        # return np.sum(node_velocity,axis=0)
        velo=np.sum(node_velocity,axis=0)
        # np.where(velo<10,velo,0)
        print(np.where(velo<10, velo, 0)) 
        # replace the value larger than 10 with 0
        # print(np.where(velo<10, velo, 0).shape)
        return np.where(velo<10, velo,0)

    def get_acceleration_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        node_acceleration=np.zeros((len(self.list_of_waves),len(list_of_point),3))
        for index, wach_wave in enumerate(self.list_of_waves):
            node_acceleration[index]=wach_wave.get_velocity_at_nodes(list_of_point,global_time)
        return np.sum(node_acceleration,axis=0)



if __name__ == "__main__":   
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    plt.rcParams['image.cmap'] = 'summer'

    sea_state=irregular_sea(4,8,3,60,45)
    print(sea_state)
    time_max = 3600 # [s]
    dt=0.1

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # time_frame=np.arange(1000,4600,dt)
    # yita=sea_state.get_elevations_with_time([0,0,0],time_frame)
    # print(np.std(yita))  
    # plt.plot(time_frame,yita)
    
    x_axis=np.array(np.arange(0,60,2)).tolist()
    position=np.zeros((5*len(x_axis),3))
    for i in range(5):
        for index, item in enumerate(x_axis):
            position[i*len(x_axis)+index]=[6*i,6*i,-item/2]
    print(position.shape)   
    print(position)     
    velocity=sea_state.get_velocity_at_nodes(position, 20000)

    velocity_mag=[]
    for each in velocity:
        velocity_mag.append(np.linalg.norm(each))
    print(velocity_mag)
    # Flatten and normalize
    velocity_mag=np.array(velocity_mag)
    normal_velo = (velocity_mag.ravel() - velocity_mag.min()) / velocity_mag.ptp()
    # print(normal_velo)
    # Repeat for each body line and two head lines
    c = np.concatenate((normal_velo, np.repeat(normal_velo, 2)))
    # Colormap
    c = plt.cm.summer(c)
    q=ax.quiver(position[:,0], 
                position[:,1],
                position[:,2],
                velocity[:,0],
                velocity[:,1],
                velocity[:,2],
                colors=c,
                normalize=False,
                )
    # q.set_array(np.linspace(0,2,10))
    fig.colorbar(q)

    for i in range(60):
        position=np.zeros((len(x_axis),3))
        for index, item in enumerate(x_axis):
            position[index]=[item,i,0]    
        yita=sea_state.get_elevation_at_nodes(position, 20000)
        ax.plot(x_axis,[i]*len(x_axis),yita,color="b")
    

    # plt.plot(x_axis,yita)
    ax.set_title("JONSWAP sea condition Hs=4m, Tp=8s r=3 ")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 60)
    ax.set_zlim(-30, 5)
    plt.savefig('./figures/waves.png', dpi=600)
    plt.show()