import matplotlib.pyplot as plt
import math

def expectedValue(vec):
    val = 0
    for i in vec:
        val = val + i
    val = val / len(vec)
    return val

f = open("D:\\cpp_projects\\rnt\\Debug\\stOutputLc.dat", "r")
#f = open("D:\\cpp_projects\\rnt\\Debug\\stOutputAdd.dat", "r")
#f = open("D:\\cpp_projects\\rnt\\Debug\\stOutput5p.dat", "r")
#f = open("D:\\cpp_projects\\rnt\\Debug\\stOutputLfsr.dat", "r")
#f = open("D:\\cpp_projects\\rnt\\Debug\\stOutputNfsr.dat", "r")
#f = open("D:\\cpp_projects\\rnt\\Debug\\stOutputMt.dat", "r")
#f = open("D:\\cpp_projects\\rnt\\Debug\\stOutputRc4.dat", "r")
#f = open("D:\\cpp_projects\\rnt\\Debug\\stOutputRsa.dat", "r")
#f = open("D:\\cpp_projects\\rnt\\Debug\\stOutputBbs.dat", "r")

currElem = []
mVec = []
dVec = []
nVec = []
acc = 0
allElems = list(map(float, (f.read()).split(',')))
for a in allElems:
    currElem.append(float(a))
    if (len(currElem) % 50 == 0):
        m = expectedValue(currElem)
        dev = 0
        for i in currElem:
            dev = dev + ((i - m) ** 2)
        dev = dev / (len(currElem) - 1)
        d = math.sqrt(dev)
        acc = acc + 50
        mVec.append(round(m, 3))
        dVec.append(round(d, 3))
        nVec.append(acc)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,10))
for i in range(len(nVec)):
    ax1.scatter(nVec[i], mVec[i], c = 'b', s = 1)
    ax2.scatter(nVec[i], dVec[i], c = 'b', s = 1)
ax1.set_title("Математическое ожидание")
ax2.set_title("Среднеквадратическое отклонение")
plt.show()
