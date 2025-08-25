# HOPG X線分光器データから Ka線収率を計算するプログラム
# 指定されたエネルギー範囲でのフォトン数を積分し、立体角と反射率で補正

import pprint
import numpy as np
from time import sleep
from scipy.optimize import differential_evolution
from scipy import integrate
from scipy.interpolate import interp1d
import math as m
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as tick

def getNearestValue(list, num):
    """
    概要: リストからある値に最も近い値を返却する関数
    @param list: データ配列
    @param num: 対象値
    @return 対象値に最も近い値
    """

    # リスト要素と対象値の差分を計算し最小値のインデックスを取得
    idx = np.abs(np.asarray(list) - num).argmin()
    return idx

d = 3.357 

# solidangle = (250/m.cos(13.27)*2*m.pi*0.1/360)*(250/m.cos(13.27)*2*m.pi*0.3985/360)/(4*m.pi*250**2)
solidangle = 4*m.asin(m.sin(2.23/2*m.pi/180)*m.sin(0.4426/2*m.pi/180))

# reflectivity = 0.7   #ZYA
reflectivity = 0.3019
#reflectivity = 0.712 #ZYB
# reflectivity2 = 0.0113
reflectivity2 = 0.025/(500/np.cos(13.27*np.pi/180)*2*np.tan(2.23/2*np.pi/180))

min_E = 8.0
max_E = 8.1

# 指定されたショット番号範囲でKa線収率を計算
for i in range(43775, 43803):
	shotnum = i
	# データファイルの存在確認
	if os.path.exists('data/HOPG_G{}.csv'.format(shotnum)):
		# HOPGスペクトラムデータの読み込み
		data = np.genfromtxt('data/HOPG_G{}.csv'.format(shotnum), delimiter=',', autostrip=True)
		a1, b1 = np.hsplit(data, 2)
		E = np.ravel(a1)      # エネルギー [keV]
		Photon = np.ravel(b1) # フォトン数密度 [photons/keV]

		# バックグラウンド減算（7.9-7.99 keVの平均値を使用）
		Photon = (Photon - np.mean(Photon[(E > 7.9) & (E < 7.99)]))

		# Ka線エネルギー範囲のインデックスを取得
		min_E_idx = getNearestValue(E, min_E)
		max_E_idx = getNearestValue(E, max_E)

		# エネルギー幅を計算してフォトン数を積分
		delta_E = E[min_E_idx+1:max_E_idx+1] - E[min_E_idx:max_E_idx]
		tot = np.sum(Photon[min_E_idx:max_E_idx]*delta_E)

		# 立体角と反射率で補正してKa線収率を計算
		result = tot/reflectivity/reflectivity2/solidangle

		print('Ka yeild of {} is {:.2E} [Photon/str]'.format(shotnum, result))
