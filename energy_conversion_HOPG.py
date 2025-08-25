import pprint
import numpy as np
from time import sleep
from scipy.optimize import differential_evolution
from scipy import integrate
from scipy.interpolate import interp1d
from numpy.linalg import solve
import matplotlib.pyplot as plt
import csv

d = 3.357 

def I(t):
	return 0.297 * np.exp(-t/57.6) + 0.723

def phi(x):
	m = 0.17
	c = 2.27
	return m * x + c

filter1 = np.loadtxt('Be_500um_transmittance.dat', skiprows=2, unpack=True)
filter1 = interp1d(filter1[0,:]*1.0E-03, filter1[1,:], kind='cubic')

filter2 = np.loadtxt('CH2_50um_transmittance.dat', skiprows=2, unpack=True)
filter2 = interp1d(filter2[0,:]*1.0E-03, filter2[1,:], kind='cubic')

t = input("Input the scaning time delay : ")

t = float(t)

shot_num = 'G43801'
file = 'profile/200911_G43801_HOPG_1212'
Energy = 1.406 #kJ

data = np.genfromtxt('{}.csv'.format(file), delimiter=',', skip_header= 1, skip_footer= 0, autostrip=True)

E_def = 8.048															#基準の線のエネルギー
theta_def = np.arcsin(12.3984/(2.0 * d * E_def))						#基準の線のブラッグ角
# place = data[np.argmax(data,axis=0)[1],0]															#基準の線のIP上での位置[cm]
place = 1.005

s = [-0.01988027, -0.79709931]

# s[0] = 0.65319243
# s[1] =  -3.4443449

diff = 13 + s[0] - place - (50 - s[1])*np.tan(np.arcsin(12.3984/(2*d*E_def)))

s[0] = s[0] - diff

theta_rad = np.arctan((13 + s[0] - data[:,0])/(50 - s[1]))
thera_deg = theta_rad * 180 / np.pi

E = 12.3984 / (2*d*np.sin(theta_rad))

def dEdx(x):
	a = 12.3984
	b = 13 + s[0]
	c = 50 - s[1]
	theta = np.arctan((b-x)/c)
	def dtdx(x):
		return -c/(c**2 + (b - x)**2)
	return (a/(2*d)) * (2 * np.cos(theta))/(-1 + np.cos(2*theta)) * dtdx(x)

Photon = data[:,1] * 1000 / (I(t) * phi(E) * 0.0025 * dEdx(data[:,0])*filter1(E)*filter2(E))

out = np.array([E,Photon]).T

np.savetxt('data/HOPG_{}.csv'.format(shot_num), out, delimiter=',')

out = np.array([E,Photon/Energy]).T

np.savetxt('data_normalize/HOPG_energy_normalize_{}.csv'.format(shot_num), out, delimiter=',')

plt.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(19.20, 10.80))
ax = fig.add_subplot(111)
# ax.text(np.min(E), (np.max(Photon) - np.min(Photon))*0.8 + np.min(Photon), "{}".format(shot_num), size=35, color='black')
ax.text(np.min(E), 4.5E+5, "{}\n{} kJ".format(shot_num,Energy), size=35, color='black')
#plt.imshow(img)
ax.set_xlabel("Energy keV", size=30)
ax.set_ylabel("Photon/keV", size=30)
ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
plt.ylim([0,6E+5])
# ax.tick_params(axis='x', labelsize=18)
# ax.tick_params(axis='y', labelsize=18)
ax.plot(E,Photon)

plt.savefig("data/HOPG_{}.png".format(shot_num), bbox_inches='tight', pad_inches=0)

plt.figure()
plt.rcParams.update({'font.size': 22})
fig = plt.figure(figsize=(19.20, 10.80))
ax = fig.add_subplot(111)
# ax.text(np.min(E), (np.max(Photon) - np.min(Photon))*0.8 + np.min(Photon), "{}".format(shot_num), size=35, color='black')
ax.text(np.min(E), 4.5E+5, "{}".format(shot_num), size=35, color='black')
#plt.imshow(img)
ax.set_xlabel("Energy keV", size=30)
ax.set_ylabel("Photon/keV/kJ", size=30)
ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
plt.ylim([0,6E+5])
# ax.tick_params(axis='x', labelsize=18)
# ax.tick_params(axis='y', labelsize=18)

ax.plot(E,Photon/Energy)

plt.savefig("data_normalize/HOPG_energy_normalize_{}.png".format(shot_num), bbox_inches='tight', pad_inches=0)
