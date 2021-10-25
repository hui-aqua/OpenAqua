from vpython import *
import numpy as np

# GlowScript 3.1 VPython

scene.caption = "None"
scene.width = 900
scene.height = 900

num_ball = 5000  # total number of knots
knot = []
for i in range(num_ball):
    knot.append(sphere(pos=vector(0, 0, 0),
                       color=color.green,
                       radius=0.05,
                       make_trail=True,
                       retain=5)
                )
    knot[i].velocity_base=vector(0.5,0,0)

wallbottom = box(pos=vector(0, 0, 0), size=vector(0.05, 10, 10), color=color.gray(0.7))
walltop = box(pos=vector(10, 0, 0), size=vector(0.05, 10, 10), color=color.gray(0.7))

t_end = 600  # simulation duration is 10 s
dt = 0.1
t = 0  # initial time 0s
# velocity_base = vector(0.5,0,0)

while t < t_end:
    rate(60)
    scene.caption = "To rotate 'camera', drag with right button or Ctrl-drag.\n"\
    "To zotexture=textures.rugom, drag with middle button or Alt/Option depressed, or use scroll wheel.\n"\
    "On a two-button mouse, middle is left + right.\n"\
    "To pan left/right and up/down, Shift-drag.\n"\
    "Touch screen: pinch/extend to zoom, swipe or two-finger rotate.\n"\
    + "Time = " + str(round(t, 4)) + " s"
    for i in range(num_ball):
        if not(10>knot[i].pos.x >= 0):
            knot[i].velocity_base.x=-knot[i].velocity_base.x
        knot[i].pos += (knot[i].velocity_base+vector(0,np.random.random()-0.5,np.random.random()-0.5)) * dt
    t += dt
