import time
start = time.perf_counter()
import pandas as pd
import numpy as np
import random
from DISP_get import DISP_get
import itertools
#random.seed(66)
length = 50
Ntest = 1000
disp1= DISP_get(length, Ntest) #生成随机碱基分配顺序

##生成相应的测序信号#
file = "D://cyy/data/scripts/datasets/sequence_list16.csv" #导入待测序列csv文件
read = open(file) #打开待测序列csv文件
seq = pd.read_csv(read) #输入待生成碱基分配顺序的序列，转换为dataframe格式
seq1 = np.array(seq).T  #将待测序列转为数组，每一行为一个MH
N_disp = disp1.shape[0]  #碱基分配顺序条数
N_MH = seq1.shape[0]  #微单倍型个数,0,1,2
N_allele_max = seq1.shape[1]  #等位基因个数,0,1,2,3,4,5,6
disp_list1 = np.zeros([1, length])  #生成一个空数组,用来存放信号"""不完整"""的碱基分配顺序
disp_list2 = np.zeros([1,length])  #生成一个空数组，在phase步骤中用来存放"""完整"""的碱基分配顺序
disp_list = pd.DataFrame()  #生成一个数据框，用来存放最终可以生成完整信号，并且可以phase的碱基分配顺序，并显示出其相关系数。
signal_sum1 = np.zeros([N_disp, N_MH, N_allele_max, length]) # 创建一个空数组,用来存放所有碱基分配顺序的测序信号
sum_signal = np.zeros([N_disp, N_MH, N_allele_max, 1]) # 创建一个空数组,用来存放所有碱基分配顺序的测序信号之和
for l in range(N_disp):
    dispensation = disp1[l] #导入碱基分配顺序
    signal_sum = np.zeros([N_MH,N_allele_max,length])  # 创建一个空数组,用来存放测序信号
    for p in range(N_MH): #在每个MH中遍历，生成测序信号
        for q in range(N_allele_max):
            signal = [0]*length  #焦磷酸测序信号的长度,并在此基础上进行相应信号的计数
            seq2 = seq1[p][q] #获取每一条待测序列
            len_allele = len(seq1[p][0]) #获取每一条待测序列的碱基个数
            j = 0 #待测序列上碱基个数计数
            for i in range(1, len(dispensation)+1): #对碱基分配顺序的碱基个数进行计数（range为左闭右开）
                count = 0  #碱基分配顺序上测序信号进行计数（0-若干）
                try:
                    while seq2[j] == dispensation[i - 1]:  # 当待测序列当前碱基和分配顺序当前碱基一致时，进行信号的叠加
                        count += 1  # 焦磷酸测序信号+1
                        signal[i - 1] = count  # 信号加1
                        j = j + 1  # 待测序列往后移一位
                        if j >= len(seq2):  # 当待测序列的碱基计数超过待测序列最大长度时
                            break  # 跳出while循环
                    if j >= len(seq2):  # 当待测序列的碱基计数超过待测序列最大长度时
                        break  # 跳出for循环
                except TypeError:
                    break
            signal_sum[p, q] = signal
            sum_signal[l,p,q] = sum(signal)
    signal_sum1[l] = signal_sum

##生成信号不完整的碱基分配顺序相应信号
for b in range(N_disp):
    dispensation1 = disp1[b]
    #print("b 是： %s" %b)
    for p in range(N_MH):
        #print("p 是： %s" % p)
        len_allele = len(seq1[p][0])
        for q in range(N_allele_max):
            #print("q 是： %s" % q)
            if 0 < sum_signal[b,p,q] < len_allele:
                #print("b,p,q is:")
                #print(b,p,q)
                #print("sum_signal 是：%s" %(sum_signal[b,p,q]))
                disp_list1 = np.vstack([disp_list1, dispensation1])
                b = b-1
                break
        else:
            continue
        break
disp_list1 = np.delete(disp_list1, 0, 0)  #disp_list1是信号不完整的碱基分配顺序

##生成信号完整的碱基分配顺序
disp1_rows = disp1.view([('', disp1.dtype)]*disp1.shape[1])
disp_list1_rows = disp_list1.view([('', disp_list1.dtype)]*disp_list1.shape[1])
disp_list1 = np.setdiff1d(disp1_rows, disp_list1_rows).view(disp1.dtype).reshape(-1, disp1.shape[1])
print("信号完整的碱基分配顺序是：%s" %disp_list1)
print(len(disp_list1))

##生成信号完整的碱基分配顺序的相应信号
print("phase检验~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
N_disp1 = disp_list1.shape[0]  #碱基分配顺序条数
signal_sum3 = np.zeros([N_disp1, N_MH, N_allele_max, length]) # 创建一个空数组,用来存放所有碱基分配顺序的测序信号
sum_signal1 = np.zeros([N_disp1, N_MH, N_allele_max, 1]) # 创建一个空数组,用来存放所有碱基分配顺序的测序信号之和
for l in range(N_disp1):
    dispensation2 = disp_list1[l] #导入碱基分配顺序
    signal_sum2 = np.zeros([N_MH,N_allele_max,length])  # 创建一个空数组,用来存放测序信号
    for p in range(N_MH): #在每个MH中遍历，生成测序信号
        for q in range(N_allele_max):
            signal2 = [0]*length  #焦磷酸测序信号的长度,并在此基础上进行相应信号的计数
            seq3 = seq1[p][q] #获取每一条待测序列
            len_allele2 = len(seq1[p][0]) #获取每一条待测序列的碱基个数
            j = 0 #待测序列上碱基个数计数
            for i in range(1, len(dispensation2)+1): #对碱基分配顺序的碱基个数进行计数（range为左闭右开）
                count = 0  #碱基分配顺序上测序信号进行计数（0-若干）
                try:
                    while seq3[j] == dispensation2[i - 1]:  # 当待测序列当前碱基和分配顺序当前碱基一致时，进行信号的叠加
                        count += 1  # 焦磷酸测序信号+1
                        signal2[i - 1] = count  # 信号加1
                        j = j + 1  # 待测序列往后移一位
                        if j >= len(seq3):  # 当待测序列的碱基计数超过待测序列最大长度时
                            break  # 跳出while循环
                    if j >= len(seq3):  # 当待测序列的碱基计数超过待测序列最大长度时
                        break  # 跳出for循环
                except TypeError:
                    break
            signal_sum2[p, q] = signal2
            sum_signal1[l,p,q] = sum(signal2)
    signal_sum3[l] = signal_sum2
#print("最终最终数组格式的信号汇总是: %s" % (signal_sum3))
#print(signal_sum3.shape)

##进行phase检验
n = []  #用来存放信号加和的大列表
dim0 = signal_sum3.shape[0]  #第一维度，即碱基分配顺序个数，2
dim1 = signal_sum3.shape[1]  #第二维度，即MH个数，3
for i in range(dim0):  #在每个碱基分配顺序下
    x = 0  ##用来对能够phase的MH进行计数，x = dim1时，认为可以phase
    for j in range(dim1):  #在每个MH下
        d = []  # 用来存放信号加和的中列表
        e = signal_sum3[i][j]
        #print(e)
        if np.all(e[-1] == 0):
            print("第%s个碱基分配顺序" % i)
            print("第%s个MH" % j)
            h = []  # 用来存放含全零行的信号加和小列表
            a1_rows = e.view([('', e.dtype)] * e.shape[1])
            a2_rows = e[-1].view([('', e[-1].dtype)] * e[-1].shape[0])
            c = np.setdiff1d(a1_rows, a2_rows).view(e.dtype).reshape(-1, e.shape[1])
            f = list(itertools.combinations(c, 2))
            #print("f是：%s" %f)
            #print("等位基因组合个数为： %s" %len(f))
            for l in f:
                add = l[0] + l[1]
                d.append(add)
            #print(d)
            #print(len(d))
            w1 = []  #对“一个MH下各个等位基因信号两两加和”进行去重
            for y in d:
                #print(y)
                #print(type(y))
                y = y.tolist()
                #print(y)
                #print(type(y))
                if y not in w1:
                    w1.append(y)
            if len(w1) == len(d):
                print("可以实现phase!")
                x = x+1
            else:
                print("不可以实现phase!")



        else:
            print("第%s个碱基分配顺序" % i)
            print("第%s个MH" % j)
            g = list(itertools.combinations(e, 2))
            #print("g是：%s" % g)
            #print("等位基因组合个数为： %s" % len(g))
            for l in g:
                add = l[0] + l[1]
                d.append(add)
            #print(d)
            #print(len(d))
            w2 = []  #对“一个MH下各个等位基因信号两两加和”进行去重
            for z in d:
                #print(z)
                z = z.tolist()
                if z not in w2:
                    w2.append(z)
                    #print(w2)
            if len(w2) == len(d):
                print("可以实现phase!")
                x = x + 1  ##用来对能够phase的MH进行计数，x = dim1时，认为可以phase
            else:
                print("不可以实现phase!")

        #n.append(d)
        #print(n)
        #print(len(n))
    if x == dim1:
        disp_list2 = np.vstack([disp_list2, disp_list1[i]])  ##把可以phase的碱基分配顺序保存下来

#n.append(d)
#print("总信号如下：")
#print(n)
print("disp!!!!!!!!!!!!disp")
disp_list2 = np.delete(disp_list2, 0, 0)
print(disp_list2)
print("可用的碱基分配顺序有%s条！" %len(disp_list2))
end = time.perf_counter()
print('Running time: %s Seconds' %(end-start))