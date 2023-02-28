# 假设反应速率远大于滴加的速率，每一滴加入的酰氯都完全反应，酰氯和不同取代数的胺 反应性相同。
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation


def count_ingredients(ingre_list):
    cnt = 0
    for _ in ingre_list:
        if _ == 0:
            cnt += 1
    return cnt + 1


def simulation(T=10000, n0=1, n1=0.5):
    """
    假设速率常数很大，滴加进去的酰氯立即被反应，不同物种的取代反应等活性假设的简单模拟
    :param T: 迭代步长，数字越大曲线越平滑
    :param n0: N8的物质量，单位为mmol
    :param n1: 酰氯的物质量，单位为mmol
    :return:
    """
    N8 = [n0]
    N8C1 = [0]
    N8C2 = [0]
    N8C3 = [0]
    N8C4 = [0]
    N8C5 = [0]
    N8C6 = [0]
    Chloride = [n1]
    for num in range(T):
        t = num + 1
        last_ingredients = [N8[num], N8C1[num], N8C2[num], N8C3[num], N8C4[num], N8C5[num]]
        ingre_num = count_ingredients(last_ingredients)
        Chloride.append(Chloride[num] - 1 / T * n1)
        chloride_all = N8C1[num] + 2 * N8C2[num] + 3 * N8C3[num] + 4 * N8C4[num] + 5 * N8C5[num] + 6 * N8C6[num]
        if chloride_all[num] > 0:
            N8.append(N8[num] - 1 / T * N8[num] * n1 / ingre_num)
            N8C1.append(N8C1[num] + 1 / T * n1 * (N8[num] - N8C1[num]) / ingre_num)
            N8C2.append(N8C2[num] + 1 / T * n1 * (N8C1[num] - N8C2[num]) / ingre_num)
            N8C3.append(N8C3[num] + 1 / T * n1 * (N8C2[num] - N8C3[num]) / ingre_num)
            N8C4.append(N8C4[num] + 1 / T * n1 * (N8C3[num] - N8C4[num]) / ingre_num)
            N8C5.append(N8C5[num] + 1 / T * n1 * (N8C4[num] - N8C5[num]) / ingre_num)
            N8C6.append(N8C6[num] + 1 / T * n1 * N8C5[num] / ingre_num)
        else:
            N8.append(N8[num])
            N8C1.append(N8C1[num])
            N8C2.append(N8C2[num])
            N8C3.append(N8C3[num])
            N8C4.append(N8C4[num])
            N8C5.append(N8C5[num])
            N8C6.append(N8C6[num])
    print(Chloride)
    return N8, N8C1, N8C2, N8C3, N8C4, N8C5, N8C6


def kinetic_simulation(T=10000, n0=1, n1=0.5, k=0.0008, speed_ratio=(1, 1, 1, 1, 1, 1)):
    """
    不同物种的取代反应等活性假设，动力学模拟的迭代过程
    :param T: 迭代步长，数字越大曲线越平滑
    :param n0: N8的物质量，单位为mmol
    :param n1: 酰氯的物质量，单位为mmol
    :param k: 取代反应的速率常数
    :param speed_ratio: 速率常数的相对大小，默认等活性假设，六个元素的列表，相对参比为1
    :return:
    """
    N8 = [n0]
    N8C1 = [0]
    N8C2 = [0]
    N8C3 = [0]
    N8C4 = [0]
    N8C5 = [0]
    N8C6 = [0]
    k1, k2, k3, k4, k5, k6 = speed_ratio
    # 体系外还未滴加的酰氯
    Chloride = [n1 - 1 / T * n1]
    # 体系内剩余的，参加动力学方程的酰氯
    Chloride_sys = [1 / T * n1]
    for num in range(3*T):
        t = num + 1
        Chloride.append(Chloride[num] - 1 / T * n1)
        chloride_all = N8C1[num] + 2*N8C2[num] + 3*N8C3[num] + 4*N8C4[num] + 5*N8C5[num] + 6*N8C6[num]
        if chloride_all < n1:
            N8.append(N8[num] - Chloride_sys[num] * N8[num] * k * k1)
            N8C1.append(N8C1[num] + Chloride_sys[num] * (N8[num] - N8C1[num]) * k * k2)
            N8C2.append(N8C2[num] + Chloride_sys[num] * (N8C1[num] - N8C2[num]) * k * k3)
            N8C3.append(N8C3[num] + Chloride_sys[num] * (N8C2[num] - N8C3[num]) * k * k4)
            N8C4.append(N8C4[num] + Chloride_sys[num] * (N8C3[num] - N8C4[num]) * k * k5)
            N8C5.append(N8C5[num] + Chloride_sys[num] * (N8C4[num] - N8C5[num]) * k * k6)
            N8C6.append(N8C6[num] + Chloride_sys[num] * N8C5[num] * k)
        else:
            N8.append(N8[num])
            N8C1.append(N8C1[num])
            N8C2.append(N8C2[num])
            N8C3.append(N8C3[num])
            N8C4.append(N8C4[num])
            N8C5.append(N8C5[num])
            N8C6.append(N8C6[num])
        if Chloride[num] > 0:
            Chloride_sys.append(Chloride_sys[num] + 1 / T * n1 * (1 - N8[num] * k))
        elif Chloride_sys[num] > 0:
            Chloride_sys.append(Chloride_sys[num] * (1 - N8[num] * k))
    return N8, N8C1, N8C2, N8C3, N8C4, N8C5, N8C6


def plot_kinetic_simulation(n0=1, n1=6, k=0.0008, axis='default', speed_ratio=(1, 1, 1, 1, 1, 1)):
    """

    :param n0:
    :param n1:
    :param k:
    :param axis: 坐标轴的表示方法，default为线性，其他选项均为模拟对数
    :return:
    """
    a, b, c, d, e, f, g = kinetic_simulation(n0=n0, n1=n1, k=k, speed_ratio=speed_ratio)
    with open('kinetic_k_{}_n0_{}_n1_{}.csv'.format(k, n0, n1), 'w') as file:
        file.writelines('N8,N8R1,N8R2,N8R3,N8R4,N8R5,N8R6\n')
        for _ in range(len(a)):
            file.writelines('{:.5f},{:.5f},{:.5f},{:.5f},{:.5f},{:.5f},{:.5f}\n'.format(a[_], b[_], c[_], d[_], e[_], f[_], g[_]))
    plt.figure(figsize=(8, 4), dpi=300)
    plt.subplots_adjust(left=0.1, right=0.93, top=0.9, bottom=0.13)
    # plt.plot(a, label='N8')
    plt.plot(np.array(b)*100, label='N8C1')
    plt.plot(np.array(c)*100, label='N8C2')
    plt.plot(np.array(d)*100, label='N8C3')
    plt.plot(np.array(e)*100, label='N8C4')
    plt.plot(np.array(f)*100, label='N8C5')
    plt.plot(np.array(g)*100, label='N8C6')
    plt.legend()
    plt.ylabel('Theoretical yield(or conversion)(%)')
    if axis != 'default':
        plt.xscale('symlog')
        plt.xlim([100, 10000])
        plt.xticks(ticks=[100, 200, 300, 500, 800, 1000, 1200, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000],
                   labels=[0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10])
        plt.xlabel('Molar Ratio (n1/n0)')
        plt.title('Kinetic simulation at different molar ratio')
        plt.grid()
        plt.savefig('together_k_{}.png'.format(k, n0, n1))
    else:
        plt.xlim([0, 10000])
        plt.xlabel('Progress of chloride addition(10000 pieces)')
        plt.title('Kinetic simulation at molar ratio {}:{}(k={})'.format(n0, n1, k))
        plt.savefig('kinetic_k_{}_n0_{}_n1_{}.png'.format(k, n0, n1))
    plt.cla()


def conversion_ratio(n0=1, n1=0.5, species=1):
    """
    Calculate the conversion rate of different species in the reaction
    :param n0: the amount of N8 in mole unit
    :param n1: the amount of chloride in mole unit
    :param species: the species index wanted, given in indexes, e.g. 1 means N8R1, 2 means N8R2.
    :return:
    """
    res = kinetic_simulation(n0=n0, n1=n1)
    sum = np.array(res[6])
    for _ in range(5):
        sum += np.array(res[_+1])
    ans = np.divide(np.array(res[species]), sum, out=np.zeros_like(np.array(res[species])), where=sum != 0) * 100
    return ans


def plot_simulation(n0=1, n1=6):
    a, b, c, d, e, f, g = kinetic_simulation(n0=n0, n1=n1)
    plt.subplots_adjust(left=0.1, right=0.93, top=0.9, bottom=0.1)
    # plt.plot(a, label='N8')
    plt.plot(np.array(b)*100, label='N8C1')
    plt.plot(np.array(c)*100, label='N8C2')
    plt.plot(np.array(d)*100, label='N8C3')
    plt.plot(np.array(e)*100, label='N8C4')
    plt.plot(np.array(f)*100, label='N8C5')
    plt.plot(np.array(g)*100, label='N8C6')
    plt.legend()
    plt.xlabel('Progress of chloride addition(10000 pieces)')
    plt.ylabel('Theoretical yield(or conversion)(%)')
    plt.xlim([0, 10000])
    plt.title('Recursive simulation at molar ratio {}:{}'.format(n0, n1))
    # plt.show()
    plt.savefig('Recursive_simulation_n0_{}_n1_{}.png'.format(n0, n1))
    plt.cla()


# for _ in [0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 20]:
#     plot_kinetic_simulation(n0=1, n1=_)

# 绘制不同速率常数下各个投料比的理论转化率动画

def animate_speed(step=np.logspace(-4, -2, 100)):
    fig, ax = plt.subplots(figsize=(12, 9))
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
    a, b, c, d, e, f, g = kinetic_simulation(n0=1, n1=10, k=step[0])
    lnb, = plt.plot(np.array(b) * 100, label='N8C1', animated=True)
    lnc, = plt.plot(np.array(c) * 100, label='N8C2', animated=True)
    lnd, = plt.plot(np.array(d) * 100, label='N8C3', animated=True)
    lne, = plt.plot(np.array(e) * 100, label='N8C4', animated=True)
    lnf, = plt.plot(np.array(f) * 100, label='N8C5', animated=True)
    lng, = plt.plot(np.array(g) * 100, label='N8C6', animated=True)
    ax.set_xlabel('Molar Ratio (n1/n0)', fontsize=20)
    ax.set_ylabel('Theoretical yield(or conversion)(%)', fontsize=20)
    ax.set_title('Kinetic simulation at K = {:.6f}'.format(step[0]), fontsize=25)
    ax.set_xscale('symlog')
    ax.set_xlim([100, 10000])
    ax.set_ylim([0, 100])
    plt.xticks(ticks=[100, 200, 300, 500, 800, 1000, 1200, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000],
                labels=[0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]
                )
    plt.grid()
    fig.legend(loc=(0, 50))

    def update(data):
        a, b, c, d, e, f, g = data[0]
        k_ = data[-1]
        lnb.set_ydata(np.array(b) * 100)
        lnc.set_ydata(np.array(c) * 100)
        lnd.set_ydata(np.array(d) * 100)
        lne.set_ydata(np.array(e) * 100)
        lnf.set_ydata(np.array(f) * 100)
        lng.set_ydata(np.array(g) * 100)
        ax.set_title('Kinetic simulation at K = {:.6f}'.format(k_))
        return lnb,

    def data_gen():
        for _ in step[1:]:
            yield kinetic_simulation(n0=1, n1=10, k=_), _

    ani = animation.FuncAnimation(fig,
                                  update,
                                  frames=data_gen,
                                  interval=15,
                                  blit=True)
    ani.save('test_animation.gif', writer='imagemagick')


# animate_speed()

plot_kinetic_simulation(n0=1, n1=10, axis='new', speed_ratio=(1, 0.9, 0.85, 0.8, 0.75, 0.7))

# 给定速率常数下绘制不同摩尔比下的N8R1转化率随滴加进度变化
plt.subplots_adjust(left=0.1, right=0.93, top=0.9, bottom=0.1)
# plt.plot(a, label='N8')
plt.plot(np.array(kinetic_simulation(n0=1, n1=0.2)[1])*100, label='1:0.2')
plt.plot(np.array(kinetic_simulation(n0=1, n1=0.5)[1])*100, label='1:0.5')
plt.plot(np.array(kinetic_simulation(n0=1, n1=0.8)[1])*100, label='1:0.8')
plt.plot(np.array(kinetic_simulation(n0=1, n1=1)[1])*100, label='1:1')
plt.plot(np.array(kinetic_simulation(n0=1, n1=2)[1])*100, label='1:2')
plt.plot(np.array(kinetic_simulation(n0=1, n1=4)[1])*100, label='1:4')
plt.plot(np.array(kinetic_simulation(n0=1, n1=6)[1])*100, label='1:6')
plt.legend()
plt.xlabel('Progress of chloride addition(10000 pieces)')
plt.ylabel('Theoretical yield(or conversion)(%)')
plt.xlim([0, 10000])
plt.title('Conversion of N8R1 at different molar ratio.(k={})'.format(0.0008))
# plt.show()
plt.savefig('MolarRatioComparison.png')
plt.cla()

# 给定速率常数下绘制不同摩尔比下的N8R1含量对于所有转化产物的比例
plt.subplots_adjust(left=0.1, right=0.93, top=0.9, bottom=0.1)
plt.plot(np.array(conversion_ratio(n0=1, n1=0.2)), label='1:0.2')
plt.plot(np.array(conversion_ratio(n0=1, n1=0.5)), label='1:0.5')
plt.plot(np.array(conversion_ratio(n0=1, n1=0.8)), label='1:0.8')
plt.plot(np.array(conversion_ratio(n0=1, n1=1)), label='1:1')
plt.plot(np.array(conversion_ratio(n0=1, n1=2)), label='1:2')
plt.plot(np.array(conversion_ratio(n0=1, n1=4)), label='1:4')
plt.plot(np.array(conversion_ratio(n0=1, n1=6)), label='1:6')
plt.legend()
plt.xlabel('Progress of chloride addition(10000 pieces)')
plt.ylabel('Conversion composition(%)')
plt.xlim([0, 10000])
plt.ylim([0, 100])
plt.title('Composition of N8R1 in all conversions at different molar ratio.')
# plt.show()
plt.savefig('ConversionComposition.png')
plt.cla()
"""
res1 = simulation(n0=1, n1=1)
res2 = simulation(n0=1, n1=0.8)
res3 = simulation(n0=1, n1=0.6)
res4 = simulation(n0=1, n1=0.4)
plt.plot(res1[1], label='1:1')
plt.plot(res2[1], label='1:0.8')
plt.plot(res3[1], label='1:0.6')
plt.plot(res4[1], label='1:0.4')
plt.legend()
plt.show()
"""