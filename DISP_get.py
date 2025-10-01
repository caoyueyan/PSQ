import numpy as np
import random
'''这个函数的作用就是根据需要的碱基分配顺序长度和总测试次数，随机生成相应的碱基分配顺序，生成的结果的类型是numpy数组'''
def DISP_get(length, Ntest):
    disp1 = np.zeros((1, length))  # 先生成一个空数组
    for i in range(Ntest):
         nucleotides = ["A", "T", "C", "G"]
         disp = [random.sample(nucleotides, 1)]
         for j in range(1, length):
             lastone = disp[-1]
             thedif = [x for x in nucleotides if x not in lastone]
             c = random.sample(thedif, 1)
             disp.append(c)
             disp = [i for k in disp for i in k]  # 将列表中的其他列表去掉，只包含一个大列表
         disp1 = np.vstack([disp1, disp])  # 将新生成的碱基分配顺序列表添加到碱基分配顺序汇总数组中。
    disp1 = np.delete(disp1, 0, 0)
    return(disp1)
    print(disp1)



