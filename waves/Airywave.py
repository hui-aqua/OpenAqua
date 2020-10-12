"""
-------------------------------------------\n
-         University of Stavanger          \n
-         Hui Cheng (PhD student)          \n
-          Lin Li (Medveileder)            \n
-      Prof. Muk Chen Ong (Supervisor)     \n
-------------------------------------------\n
Any questions about this code,
please email: hui.cheng@uis.no \n
wave direction is x+
"""
import numpy as np
from numpy import pi


class Airywave:
    """
    Using Airy wave theory      \n
    Ref. DNV GL-RP205 Ver. 2008:P45
    Linear wave theory.
    """

    def __init__(self, wave_height=1.0, wave_period=10.0, water_depth=60.0, direction=0.0, phase=0.0):
        """
        :param wave_height: [float] Unit: [m]. wave height.
        :param wave_period: [float] Unit: [s]. wave period.
        :param water_depth: [float] Unit: [m]. wave depth.
        :param direction： [float] Unit: [degree]. direction of propagation, measured from the positive x-axis.
        :param phase: [float] Unit: [degree]. phase.
        """
        self.gravity = 9.81
        self.wave_Height = wave_height
        self.wave_Period = wave_period
        self.water_Depth = water_depth
        self.wave_beta = pi * direction / 180.0
        self.wave_phase = pi * phase / 180.0

        # 1 Calculation
        # wave length
        alpha = [1, 0.666, 0.445, -0.105, 0.272]
        omega_ba = 4.0 * pow(pi, 2) * water_depth / self.gravity / pow(wave_period, 2)
        f_omega = 0.0
        for index, item in enumerate(alpha):
            f_omega += item * pow(omega_ba, index)
        self.wave_Length = self.wave_Period * pow(self.gravity * self.water_Depth, 0.5) * pow(f_omega / (1 + omega_ba * f_omega), 0.5)
        # wave number
        self.wave_k = 2 * pi / self.wave_Length
        # angular frequency
        self.omega = pow(self.gravity * self.wave_k * np.tanh(self.wave_k * self.water_Depth), 0.5)
        # phase velocity
        self.wave_phase_velocity = pow(self.gravity / self.wave_k * np.tanh(self.wave_k * self.water_Depth), 0.5)
        # self.wavep2=pow(9.81/2/pi/self.wave_Length*np.tanh(2*pi*water_depth/self.wave_Length),-0.5)
        # print("waveperiod is " +str())

        # 2 for easy calculation
        self.pi_h_t = pi * wave_height / wave_period
        self.pi_h_t_2 = 2 * wave_height * pow(pi / wave_period, 2)

    def __str__(self):
        """
        Print the information of object
        """
        s0 = 'The environment is airy wave condition and the specific parameters are:\n'
        s1 = 'water Depth= ' + str(self.water_Depth) + ' m\n'
        s2 = 'wave Period= ' + str(self.wave_Period) + ' s\n'
        s2_1 = 'wave number= ' + str(self.wave_k) + ' \n'
        s3 = 'wave Length= ' + str(self.wave_Length) + ' m\n'
        s4 = 'wave Height= ' + str(self.wave_Height) + ' m\n'
        s5 = 'wave phase velocity= ' + str(self.wave_phase_velocity) + ' m/s\n'
        s6 = 'wave direction：= ' + str(self.wave_beta) + ' degrdd\n'
        S = s0 + s1 + s2 + s2_1+ s3 + s4 + s5 + s6
        return S

    def calc_theta(self, position, global_time):
        """
        A private function. \n
        :param position: [np.array].shape=(n,3) or [np.array].shape=(n,2) coordinates Unit: [m]. \n The coordinate of the point which you want to know the wave surface elevation. can be [x,y] or [x,y,z]
        :param global_time: [float] Unit: [s].
        :return: [float] Unit: [m]. The sea surface level in Z direction. At the targeted position.
        """
        if len(position.shape) == 1:
            # only one point
            return self.wave_k * (position[0] * np.cos(self.wave_beta) + position[1] * np.sin(self.wave_beta)) - self.omega * global_time + self.wave_phase
        elif len(position.shape) == 2:
            # a list of point
            return self.wave_k * (position[:, 0] * np.cos(self.wave_beta) + position[:, 1] * np.sin(self.wave_beta)) - self.omega * global_time + self.wave_phase

    def get_elevation(self, position, global_time):
        """
        A private function. \n
        :param position: [np.array].shape=(n,3) coordinates Unit: [m]. The position of the point which you want to know the wave surface elevation.
        :param global_time: [float] Unit: [s].
        :return: scale or [np.array].shape=(n,) Unit: [m]. The sea surface level in Z direction.
        """
        return self.wave_Height / 2 * np.cos(self.calc_theta(position, global_time))

    def get_velocity_with_time(self, position, global_time):
        """
        :param position: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave velocity
        :param global_time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return:  [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the velocity at the targeted point.
        """
        theta = self.calc_theta(position, global_time)
        eta = self.get_elevation(position, global_time)
        # wheeler's method
        # z_wheeler = (position[2] - self.wave_Height / 2) / (1 + self.wave_Height / 2 / self.water_Depth)
        z_wheeler = position[2]
        velocity=np.zeros((len(global_time),3))
        velocity[:, 0] = np.cos(self.wave_beta) * self.pi_h_t * np.cosh(self.wave_k * (z_wheeler + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        velocity[:, 1] = np.sin(self.wave_beta) * self.pi_h_t * np.cosh(self.wave_k * (z_wheeler + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        velocity[:, 2] =                          self.pi_h_t * np.sinh(self.wave_k * (z_wheeler + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        for i in range(len(eta)):
            if position[2] > eta[i]:
                velocity[i] = 0.0

        # if position[2] < eta:
        #     horizontal_velocity = self.pi_h_t * np.cosh(self.wave_k * (z_wheeler + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        #     vertical_velocity = self.pi_h_t * np.sinh(self.wave_k * (z_wheeler + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        # else:
        #     horizontal_velocity = 0.0
        #     vertical_velocity = 0.0
        # velocity = np.array([0.0, 0.0, 0.0])
        # velocity[0] = horizontal_velocity * np.cos(self.wave_beta)
        # velocity[1] = horizontal_velocity * np.sin(self.wave_beta)
        # velocity[2] = vertical_velocity
        # if np.linalg.norm(velocity) > 5:
        #     print("Warning!! Velocity at " + str(position) + " is very large as " + str(velocity))
        return velocity

    def get_acceleration_with_time(self, position, global_time):
        """
        :param position: [np.array].shape=(1,3) | Unit: [m]. The position of the point which you want to know the wave acceleration.
        :param global_time: [float] Unit: [s]. The time which you want to know the wave velocity
        :return: [np.array].shape=(1,3) Unit: [m/s]. A numpy array of the acceleration at the targeted point.
        """
        theta = self.calc_theta(position, global_time)
        eta = self.get_elevation(position, global_time)
        # wheeler's method
        z_streched = (position[2] - self.wave_Height / 2) / (1 + self.wave_Height / 2 / self.water_Depth)
        # z_streched=position[2]
        acceleration = np.zeros((len(global_time), 3))
        acceleration[:, 0] = np.cos(self.wave_beta) * self.pi_h_t_2 * np.cosh(
            self.wave_k * (z_streched + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        acceleration[:, 1] = np.sin(self.wave_beta) * self.pi_h_t_2 * np.cosh(
            self.wave_k * (z_streched + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        acceleration[:, 2] = -self.pi_h_t_2 * np.sinh(self.wave_k * (z_streched + self.water_Depth)) * np.cos(theta) / np.sinh(
            self.wave_k * self.water_Depth)
        for i in range(len(global_time)):
            if position[2] > eta[i]:
                acceleration[i] = 0.0
        return acceleration

    def get_elevations_with_time(self, position, time_list):
        """
        Public function.\n
        :param position: [np.array].shape=(n,3) Unit: [m]. The position of one node
        :param time_list: [np.array].shape=(n,1) | Uint: [s]. The time sequence for geting the elevations
         \n
        :return: Get a list of elevations at one position with a time squence \n
        """
        return self.get_elevation(position, time_list)

    def get_elevation_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of elevation at a list of point \n
        """
        return self.get_elevation(list_of_point, global_time)

    def get_velocity_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point:  [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: [float] Unit: [s]. Physical time.
        :return: Get a list of velocity at a list of point\n
        """
        theta = self.calc_theta(list_of_point, global_time)
        eta = self.get_elevation(list_of_point, global_time)
        # wheeler streching method
        # z_wheeler = (list_of_point[:, 2] - self.wave_Height / 2) / (1 + self.wave_Height / 2 / self.water_Depth)
        z_wheeler=list_of_point[:,2]
        # print("Here the z is "+str(z_wheeler))
        velocity = np.zeros((len(list_of_point), 3))
        velocity[:, 0] = np.cos(self.wave_beta) * self.pi_h_t * np.cosh(self.wave_k * (z_wheeler + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        velocity[:, 1] = np.sin(self.wave_beta) * self.pi_h_t * np.cosh(self.wave_k * (z_wheeler + self.water_Depth)) * np.cos(theta) / np.sinh(self.wave_k * self.water_Depth)
        velocity[:, 2] =                          self.pi_h_t * np.sinh(self.wave_k * (z_wheeler + self.water_Depth)) * np.sin(theta) / np.sinh(self.wave_k * self.water_Depth)
        # for i in range(len(list_of_point)):
        #     if list_of_point[i, 2] > eta[i]:
        #         velocity[i] = 0.0
            # elif np.linalg.norm(velocity[i]) >0.1:
            #     print(np.linalg.norm(velocity[i]))
            #     print(list_of_point[i])
            #     print(self.pi_h_t )
            #     print(np.cosh(self.wave_k * (z_wheeler[i] + self.water_Depth)) * np.cos(theta[i]) / np.sinh(self.wave_k * self.water_Depth))
            #     print(self.wave_k * (z_wheeler[i] + self.water_Depth))
        return velocity

    def get_acceleration_at_nodes(self, list_of_point, global_time):
        """
        Public function.\n
        :param list_of_point: [np.array].shape=(n,3) Unit: [m]. A list of points's positions
        :param global_time: time [s] \n
        :return: Get a list of acceleration at a list of point \n
        """
        zeta = self.calc_theta(list_of_point, global_time)
        yita = self.get_elevation(list_of_point, global_time)
        # wheeler streching method
        z_streched = (list_of_point[:, 2] - self.wave_Height / 2) / (1 + self.wave_Height / 2 / self.water_Depth)
        # z_streched=position[2]
        acce = np.zeros((len(list_of_point), 3))
        acce[:, 0] = np.cos(self.wave_beta) * self.pi_h_t_2 * np.cosh(
            self.wave_k * (z_streched + self.water_Depth)) * np.sin(zeta) / np.sinh(self.wave_k * self.water_Depth)
        acce[:, 1] = np.sin(self.wave_beta) * self.pi_h_t_2 * np.cosh(
            self.wave_k * (z_streched + self.water_Depth)) * np.sin(zeta) / np.sinh(self.wave_k * self.water_Depth)
        acce[:, 2] = -self.pi_h_t_2 * np.sinh(self.wave_k * (z_streched + self.water_Depth)) * np.cos(zeta) / np.sinh(
            self.wave_k * self.water_Depth)
        for i in range(len(list_of_point)):
            if list_of_point[i, 2] > yita[i]:
                acce[i] = 0.0
        return acce

    # def elements_velocities(self, position_nodes, elements, global_time):
    #     """
    #     :param position_nodes: a numpy list of position \n
    #     :param elements: a python list of element \n
    #     :param global_time: time [s] \n
    #     :return: Get a numpy array of velocity at a list of elements \n
    #     """
    #     velocity_list = []
    #     for element in elements:
    #         element_center = np.array([0, 0, 0])
    #         for node in element:
    #             element_center += position_nodes[node] / len(element)
    #         velocity_on_element = self.get_velocity(element_center, global_time)
    #         velocity_list.append(velocity_on_element)
    #     return np.array(velocity_list)

    # def elements_accelerations(self, position_nodes, elements, global_time):
    #     """
    #     :param position_nodes: a numpy list of position \n
    #     :param elements: a python list of element \n
    #     :param global_time: time [s] \n
    #     :return: Get a numpy array of acceleration at a list of elements \n
    #     """
    #     acceleration_list = []
    #     for element in elements:
    #         element_center = np.array([0, 0, 0])
    #         for index in element:
    #             element_center += position_nodes[index] / len(element)
    #         acceleration_list.append(self.get_acceleration(element_center, global_time))
    #     return np.array(acceleration_list)


if __name__ == "__main__":

    ## validation 1, ref to Figure3-3 on page 47 from DNV GL-RP205 Ver. 2008
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    g1 = 1
    g2 = 2
    gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots

    water_d = [10, 20, 30, 40, 50, 60, 80, 100, 1000]
    waves_length = []
    waves_phasevelocity = []
    wave_period = np.linspace(1, 20, 50)
    for item in water_d:
        waves_length.append([Airywave(1, i, item, 0).wave_Length for i in wave_period])
        waves_phasevelocity.append([Airywave(1, i, item, 0).wave_phase_velocity for i in wave_period])

    plt.figure(figsize=(6.3, 4.0))

    ax = plt.subplot(gs[0, 0])
    for item in waves_length:
        plt.plot(wave_period, item, label="Depth " + str(water_d[waves_length.index(item)]))
    plt.xlabel("Wave period (s)")
    plt.ylabel("Wave length (m)")
    plt.xlim(0, 20)
    plt.ylim(0, 700)
    plt.grid(True)
    plt.legend()

    ax = plt.subplot(gs[0, 1])
    for item in waves_phasevelocity:
        plt.plot(wave_period, item, label="Depth " + str(water_d[waves_phasevelocity.index(item)]))
    plt.xlabel("Wave period (s)")
    plt.ylabel("Phase velocity (m/s)")
    plt.xlim(0, 20)
    plt.ylim(0, 35)
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig('./figures/waveperiod_vs_wavelengthAndphasevelocity.png', dpi=600)
    # plt.show()

    ## validation 2 shows the wave elevation according to time and space
    g1 = 2
    g2 = 1
    gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots

    water_d = [10, 20, 30, 40, 50, 60, 80, 100, 1000]

    time_slice = np.linspace(0, 100, 1000)
    velocities_with_time = []

    space_slice = np.ones((1000, 3))
    x_axis = []
    for posi in range(1000):
        space_slice[posi] = [posi / 2, 0, 0]
        x_axis.append(posi / 2)
    elevation_with_time = []
    wave_height = 1.5
    wave_period = 10

    for item in water_d:
        velocities_with_time.append([Airywave(wave_height, wave_period, item, 0).get_elevation(np.array([0, 0, 0]), i) for i in time_slice])
        elevation_with_time.append(Airywave(wave_height, wave_period, item, 0).get_elevation_at_nodes(space_slice, 0))
    plt.figure()

    ax = plt.subplot(gs[0, 0])
    plt.title("wave elevation at x=0,y=0")
    for item in water_d:
        plt.plot(time_slice, velocities_with_time[water_d.index(item)], label="Depth " + str(item))
    plt.xlabel("Time (s)")
    plt.ylabel("Wave elevation (m)")
    plt.xlim(0, 100)
    plt.ylim(-3, 3)
    plt.grid(True)
    plt.legend()

    ax = plt.subplot(gs[1, 0])
    plt.title("wave elevation when t=0")
    for item in water_d:
        plt.plot(x_axis, elevation_with_time[water_d.index(item)], label="Depth " + str(item))
    plt.xlabel("X (m)")
    plt.ylabel("Wave elevation (m)")
    plt.xlim(0, 500)
    plt.ylim(-3, 3)
    plt.grid(True)
    # plt.legend()

    plt.tight_layout()
    plt.savefig('./figures/wave_shape.png', dpi=600)
    # plt.show()


    ## validation 3 shows the wave velocity  according to time and space
    #TODO: error place
    g1 = 2
    g2 = 1
    gs = gridspec.GridSpec(g1, g2)  # Create 1x2 sub plots

    water_d = [10, 20, 30, 40, 50, 60, 80, 100, 1000]

    time_slice = np.linspace(0, 100, 1000)
    velocities_with_time = []
    elevation_with_time = []
    space_slice = np.ones((1000, 3))
    x_axis = []
    for posi in range(1000):
        space_slice[posi] = [posi / 2, 0, 0]
        x_axis.append(posi / 2)

    wave_height = 1.5
    wave_period = 10

    for item in water_d:
        velocities_with_time.append(np.linalg.norm(Airywave(wave_height, wave_period, item, 0).get_velocity_with_time(np.array([0, 0, 0]), time_slice) ,axis=1))
        elevation_with_time.append( np.linalg.norm(Airywave(wave_height, wave_period, item, 0).get_acceleration_with_time(np.array([0, 0, 0]), time_slice),axis=1))
    plt.figure()

    ax = plt.subplot(gs[0, 0])
    plt.title("wave velocity at x=0,y=0")
    for item in water_d:
        plt.plot(time_slice, velocities_with_time[water_d.index(item)], label="Depth " + str(item))
    plt.xlabel("Time (s)")
    plt.ylabel("Wave velocity (m)")
    plt.xlim(0, 100)
    plt.ylim(-3, 3)
    plt.grid(True)
    plt.legend()

    ax = plt.subplot(gs[1, 0])
    plt.title("wave velocity when x=0,y=0")
    for item in water_d:
        plt.plot(x_axis, elevation_with_time[water_d.index(item)], label="Depth " + str(item))
    plt.xlabel("X (m)")
    plt.ylabel("Time (s)")
    plt.xlim(0, 500)
    plt.ylim(-3, 3)
    plt.grid(True)
    # plt.legend()

    plt.tight_layout()
    plt.savefig('./figures/wave_velocity.png', dpi=600)
    # plt.show()





    ## validation 4 shows the wave velocity and acceleration
    plt.figure()

    wave1 = Airywave(5, 8, 600, 0, 0)
    print(wave1)
    x_list = np.linspace(0, 90, 10)
    z_list = np.linspace(5, -60, 20)

    yita_list = []
    for x in x_axis:
        yita_list.append(wave1.get_elevation(np.array([x, 0, 0]), 0))

    posi = []
    # velo=[]
    # acce=[]
    for x in x_list:
        for z in z_list:
            posi.append([x, 0, z])
            # velo.append(wave1.get_velocity(np.array([x,0,z]),0).tolist())
            # acce.append(wave1.get_acceleration(np.array([x,0,z]),0).tolist())
    posi = np.array(posi)
    # print(posi)
    velo = wave1.get_velocity_at_nodes(posi, 0)
    acce = wave1.get_acceleration_at_nodes(posi, 0)

    print("maximum velocity mag is " + str(np.max(np.linalg.norm(velo, axis=1))))
    # print("acceleration mag is"+str(np.linalg.norm(acce,axis=1)))

    ax = plt.subplot(gs[0, 0])
    plt.title("velocity")
    plt.plot(x_axis, yita_list, color="b")
    plt.quiver(posi[:, 0], posi[:, 2], velo[:, 0], velo[:, 2])

    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.xlim(-10, 100)
    plt.ylim(-60, 10)
    plt.grid(True)

    ax = plt.subplot(gs[1, 0])
    plt.title("Acceleration")
    plt.plot(x_axis, yita_list, color="b")
    plt.quiver(posi[:, 0], posi[:, 2], acce[:, 0], acce[:, 2])

    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.xlim(-10, 100)
    plt.ylim(-60, 10)
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('./figures/velocityandacceleration.png', dpi=600)
    plt.show()
