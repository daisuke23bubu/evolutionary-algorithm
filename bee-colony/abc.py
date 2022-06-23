# ABC で重回帰分析

import numpy as np
import pandas as pd
import random as rd
from sklearn.preprocessing import StandardScaler

data_x = np.array([[7.54,23.23],[3.12,23.23],[3.44,23.77],[12.52,27.00],[8.20,25.87],[5.36,26.33]])
data_y = np.array([[21.90],[23.12],[23.43],[23.77],[24.00],[25.40]])


# 説明変数の数
N = 2

# データ数
D = 6

# 食料源の初期位置
pos = np.zeros(N)
for i in range(N):
    pos[i] = rd.uniform(-1, 1)

pos_tmp = np.zeros(N)


# 標準化
sc = StandardScaler()
exSData = sc.fit_transform(data_x)
resSData = (data_y - data_y.mean()) / data_y.std()


# 評価値
def valueF(pos):
    f = 0
    g = 0
    for i in range(D):
        for j in range(N):
            g += pos[j] * exSData[i][j]
        
        f += pow((resSData[i] - g), 2)
            
    return f


# 収穫バチ
e_bee = 100
e_bee_matrix = np.zeros((e_bee, N))

# 初期の食料源 (解候補の行列)
X = np.zeros((e_bee, N))
for i in range(e_bee):
    for j in range(N):
        X[i][j] = rd.uniform(-1, 1)

# 生成し置換する
def change(matrix, bee):
    for i in range(bee):
        for j in range(N):
            m = rd.randint(0, bee-1)
            pos_tmp[j] = X[m][j]

        k = rd.randint(0, N-1)   
        pos_tmp[k] = pos_tmp[k] + rd.uniform(-1, 1) * (pos_tmp[k] - X[rd.randint(0, bee-1)][k])
        
        # 小さくなったら更新
        if valueF(pos) > valueF(pos_tmp):
            matrix[i] = pos_tmp

        else:
            matrix[i] = X[m]

    return matrix


# 追従バチ
o_bee = 10
o_bee_matrix = np.zeros((o_bee, N))

# 追従する収穫バチを選ぶ確率
def prob_b(matrix, k):
    listf = np.zeros(e_bee)

    for i in range(e_bee):
        pos_damy = np.zeros(N)

        for j in range(N):
            pos_damy[j] = matrix[i][j]

        listf[i] = valueF(pos_damy)

    vmax = listf.max()
    vmin = listf.min()

    trVal = np.zeros(e_bee)

    for i in range(e_bee):
        trVal[i] = (vmax - listf[i]) / (vmax - vmin)

    prob = trVal[k] / trVal.sum()

    return prob

def onlooker_bee(o_bee_matrix):
    # ハチ用のリスト
    bee_list = []
    for i in range(e_bee):
        bee_list.append(i)
        
    # 重み用のリスト
    weight_list = []
    for i in range(e_bee):
        weight_list.append(prob_b(e_bee_matrix, i))
        
    index_list = rd.choices(bee_list, k=o_bee, weights=weight_list)
    
    for i in range(o_bee):
        for j in range(N):
            o_bee_matrix[i][j] = e_bee_matrix[index_list[i]][j]

    return o_bee_matrix



# 偵察バチ (scout bee)
def scout_bee(X):
    for i in range(e_bee):
        if (np.all(X[i] == e_bee_matrix[i])):
            for j in range(N):
                X[i][j] = rd.uniform(-1, 1)

    return X

# 最良の食料源を更新する
def best_sol(pos):
    best_pos = np.zeros(N)
    for i in range(N):
        best_pos[i] = o_bee_matrix[0][i]
    
    if valueF(pos) > valueF(best_pos):
        pos = best_pos

    return pos

    


if __name__ == "__main__":

    # 終了条件
    esp = pow(10, -5)

    # 最大反復回数
    max_generation = 1000
    
    for generation in range(max_generation):

        # 収穫バチフェーズ
        e_bee_matrix = change(e_bee_matrix, e_bee)

        # 追従バチフェーズ
        onlooker_bee(o_bee_matrix)
        o_bee_matrix = change(o_bee_matrix, o_bee)

        # 偵察バチフェーズ
        X = scout_bee(X)

        # 最良食料源の記録
        pos = best_sol(pos)

        if valueF(pos) < esp:
            print(generation)
            break

        print(f"{generation}世代:偏回帰係数={pos}")
        