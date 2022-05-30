# -*- coding: utf-8 -*-

from unittest import result
import random
import GeneticAlgorithm as ga
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt

def create_genom(length):
    """
    :param length: 遺伝子情報の長さ
    :return: 生成した個体集団genomClass
    """
    genome_list = []
    for i in range(length):
        genome_list.append(random.randint(0, 1))
    return ga.genom(genome_list, 0)


def evaluation(ga, number):
    """評価関数
    :param ga: 評価を行うgenomClass
    :return: 評価処理をしたgenomClass
    """
    valueA = 0
    valueB = 0
    genom_list = ga.getGenom()
    for i,state in enumerate(genom_list):
        if state == 0:
            valueA += np.sqrt(i+1)
        else:
            valueB += np.sqrt(i+1)
    value = (-1)*np.abs(valueA - valueB)

    return Decimal(value)


def select(ga, elite):
    """エリート選択
    評価が高い順番にソート
    :param ga: 選択を行うgenomClassの配列
    :return: 選択処理をした一定のエリート、genomClass
    """
    # 現行世代個体集団の評価を高い順番にソート
    sort_result = sorted(ga, reverse=True, key=lambda u: u.evaluation)
    # 一定の上位を抽出
    result = [sort_result.pop(0) for i in range(elite)]
    return result


def crossover(ga_one, ga_second):
    """二点交叉
    :param ga: 交叉させるgenomClassの配列
    :return: 二つの子孫genomClassを格納したリスト
    """
    # 子孫を格納するリストを生成
    genom_list = []
    # 入れ替える二点の点を設定
    cross_one = random.randint(0, GENOM_LENGTH)
    cross_second = random.randint(cross_one, GENOM_LENGTH)
    # 遺伝子を取り出す
    one = ga_one.getGenom()
    second = ga_second.getGenom()
    # 交叉
    progeny_one = one[:cross_one] + second[cross_one:cross_second] + one[cross_second:]
    progeny_second = second[:cross_one] + one[cross_one:cross_second] + second[cross_second:]
    # genomClassインスタンスを生成して子孫をリストに格納
    genom_list.append(ga.genom(progeny_one, 0))
    genom_list.append(ga.genom(progeny_second, 0))
    return genom_list


def next_generation_gene_create(ga, ga_elite, ga_progeny):
    """世代交代
    :param ga: 現行世代個体集団
    :param ga_elite: 現行世代エリート集団
    :param ga_progeny: 現行世代子孫集団
    :return: 次世代個体集団
    """
    # 現行世代個体集団の評価を低い順番にソート
    next_generation_geno = sorted(ga, reverse=False, key=lambda u: u.evaluation)
    # 追加するエリート集団と子孫集団の合計ぶんを取り除く
    for i in range(0, len(ga_elite) + len(ga_progeny)):
        next_generation_geno.pop(0)
    # エリート集団と子孫集団を次世代集団を次世代へ追加
    next_generation_geno.extend(ga_elite)
    next_generation_geno.extend(ga_progeny)
    return next_generation_geno


def mutation(ga, induvidual_mutation, genom_mutation):
    """突然変異
    :param ga: genomClass
    :return: 突然変異処理をしたgenomClass
    """
    ga_list = []
    for i in ga:
        # 個体に対して一定の確率で突然変異
        if induvidual_mutation > (random.randint(0, 100) / Decimal(100)):
            genom_list = []
            for i_ in i.getGenom():
                # 個体の遺伝子情報一つ一つに対して突然変異
                if genom_mutation > (random.randint(0, 100) / Decimal(100)):
                    genom_list.append(random.randint(0, 1))
                else:
                    genom_list.append(i_)
            i.setGenom(genom_list)
            ga_list.append(i)
        else:
            ga_list.append(i)
    return ga_list


if __name__ == '__main__':
    # 遺伝子情報の長さ
    GENOM_LENGTH = 64
    # 遺伝子集団の大きさ
    MAX_GENOM_LIST = 1000
    # エリート選択数
    SELECT_GENOM = 1
    # 個体突然変異確率
    INDIVIDUAL_MUTATION = 0.01
    # 遺伝子突然変異確率
    GENOM_MUTATION = 0.01
    # 繰り返す世代数
    MAX_GENERATION = 1000

    #  結果プロット
    y_max = []
    y_min = []
    y_mean = []
    y_median = []
    max_value = 0

    # 一番最初の現行世代個体集団を生成
    current_generation_individual_group = []
    for i in range(MAX_GENOM_LIST):
        current_generation_individual_group.append(create_genom(GENOM_LENGTH))

    for count_ in range(1, MAX_GENERATION + 1):
        # 現行世代個体集団の遺伝子を評価し、genomClassに代入
        for i in range(MAX_GENOM_LIST):
            evaluation_result = evaluation(current_generation_individual_group[i], GENOM_LENGTH)
            current_generation_individual_group[i].setEvaluation(evaluation_result)
        # エリート個体を選択
        elite_genes = select(current_generation_individual_group,SELECT_GENOM)
        # エリート遺伝子を交叉させ、リストに格納
        progeny_gene = []
        next_generation_individual_group = next_generation_gene_create(current_generation_individual_group,
                                                                       elite_genes, progeny_gene)
        # 次世代個体集団全ての個体に突然変異
        next_generation_individual_group = mutation(next_generation_individual_group,INDIVIDUAL_MUTATION,GENOM_MUTATION)

        # 1世代の進化的計算終了

        # 各個体適用度を配列化
        fits = [i.getEvaluation() for i in current_generation_individual_group]

        # 進化結果を評価
        min_    = min(fits)
        max_    = max(fits)
        avg_    = sum(fits) / Decimal(len(fits))
        median_ = np.median(fits)
        # 最大値が更新されたときのみ表示
        if max_value < max_ or count_ == 1:
            print ("-----第{}世代の結果-----".format(count_))
            print ("  Min   :{}".format(min_))
            print ("  Max   :{}".format(max_))
            print ("  Avg   :{}".format(avg_))
            print ("  Median:{}".format(avg_))
            max_value = max_

        y_max.append(max_)
        y_min.append(min_)
        y_mean.append(avg_)
        y_median.append(median_)
        # 現行世代と次世代を入れ替え
        current_generation_individual_group = next_generation_individual_group

    # 最終結果出力
    print ("最も優れた個体は{}".format(elite_genes[0].getGenom()))

    x = np.arange(0,1000)
    plt.plot(x, y_mean, label="mean", color="green")
    plt.plot(x, y_median,label="median", color="blue")
    plt.plot(x, y_min, label="min", color="red")
    plt.plot(x, y_max, label="max", color="cyan")
    plt.xlabel("$Generation$", fontsize=15)
    plt.ylabel("$Objective Function$", fontsize=15)
    plt.legend()
    plt.savefig("analysis_gene.jpg")
    plt.show()
