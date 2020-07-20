# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 11:41:26 2020

@author: tannuri
"""
import numpy as np
import control as co_general
import matplotlib.pyplot as plt
import control.matlab as co
from matplotlib import animation
plt.close('all')

# exercicio de atraso / PI
#G = co.tf(1,[1,3,2]) * co.tf(1,[1,10])
#Gc = 64 
##Gc2 = co.tf(64*np.asarray([1,0.1]),[1,0])
#Gc2 = co.tf(64*10*np.asarray([10,1]),[100,1])

# exercicio de avanço / PD
G = co.tf(1,[1,10,24,0])
Gc = 63 
Gc2 = co.tf(100*np.asarray([1/3,1]),[1/12,1]) #avanco
#Gc2 = co.tf(191*np.asarray([1/6,1]),[1]) #PD


T = np.linspace(0,20,3000)
co.damp(co.feedback(G*Gc,1))
y,t = co.step(co.feedback(G*Gc,1),T)
u,t = co.step(co.feedback(Gc,G),T)

co.damp(co.feedback(G*Gc2,1))
y2,t = co.step(co.feedback(G*Gc2,1),T)
u2,t = co.step(co.feedback(Gc2,G),T)

plt.subplot(2,1,1)
plt.plot(t,y,t,y2)
plt.grid()
plt.legend(['y1','y2'])
plt.subplot(2,1,2)
plt.plot(t,u,t,u2)
plt.grid()
plt.legend(['u1','u2'])





# set up figure and animation
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-1, 5), ylim=(-1, 5))
ax2 = fig.add_axes([0.88, 0.11, 0.07, 0.8])
ax2.set_ylim([-5,5])
ax2.set_xlim([-1,3])

ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
line2, = ax.plot([], [], 'o-', lw=2)
energia1, = ax2.plot([], [], 'o-', lw=2)
energia2, = ax2.plot([], [], 'o-', lw=2)

time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def init():
    """initialize animation"""
    line.set_data([], [])
    line2.set_data([], [])
    energia1.set_data([], [])
    energia2.set_data([], [])
    time_text.set_text('')
    return line, line2, energia1, energia2, time_text

def animate(i):
    """perform animation step"""
    line.set_data([0,5*np.cos(y[i]*60*np.pi/180)],[0,5*np.sin(y[i]*60*np.pi/180)])
    line2.set_data([0,5*np.cos(y2[i]*60*np.pi/180)],[0,5*np.sin(y2[i]*60*np.pi/180)])
    energia1.set_data([1,1],[0,u[i]])
    energia2.set_data([2,2],[0,u2[i]])
    time_text.set_text('time = %.1f' % t[i])
    return line, line2, energia1, energia2, time_text

ani = animation.FuncAnimation(fig, animate, frames=3000,
                              interval=0.01, blit=True, init_func=init)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#ani.save('compara.mp4', fps=30)

plt.show()