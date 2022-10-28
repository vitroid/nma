#!/usr/bin/env python
# coding: utf-8

# 積極的に小さい関数を作り、テストをしやすくする。

# In[1]:


import numpy as np
import math

# TIP4Pの座標系

angle,holen,colen = 104.52,0.9572,0.15
wh = 1.0
wo = 16.0
wm = wh + wh + wo
rangle = angle*math.pi/360
ohz=holen*math.cos(rangle)
hyl=holen*math.sin(rangle)
hzl=wo*ohz/wm
ol=-ohz+hzl
cl=colen+ol

# 分子内の原子位置
intra = np.array([[0,hyl,hzl],
                  [0,-hyl,hzl],
                  [0,0,ol],
                  [0,0,cl]]).T

# 分子内座標の、軸を交換したもの
intra2 = np.zeros_like(intra)
intra2[0] = intra[1]
intra2[1] = intra[2]
intra2[2] = intra[0]


# In[2]:


# 分子の座標を読みこむ。
file = open("test.nx3a")
# file = open("1234.0.nx3a")
while True:
    line = file.readline()
    if len(line) == 0:
        break
    if "@NX3A" in line:
        line = file.readline()
        nmol = int(line)
        pos = []
        for i in range(nmol):
            line = file.readline()
            pos.append([float(x) for x in line.split()])
        pos = np.array(pos)
    if "@BOX3" in line:
        line = file.readline()
        cell = np.diag(np.array([float(x) for x in line.split()]))
        celli = np.linalg.inv(cell)


# In[3]:


#回転行列と、その角度微分をSympyで生成する     
from sympy import *

# おまじない
init_printing()

a,b,c=symbols('a b c')

B = Matrix([[ cos(c),-sin(c),0],
            [ sin(c), cos(c),0],
            [0,0,1]])
C = Matrix([[1,0,0],
            [0, cos(a),-sin(a)],
            [0, sin(a), cos(a)]])
D = Matrix([[ cos(b),-sin(b),0],
            [ sin(b), cos(b),0],
            [0,0,1]])
AI = D*C*B


# 解析的に回転行列を角度で微分し、それをpythonの関数とする。
A1 = lambdify([a,b,c],AI)
A1a = lambdify([a,b,c],diff(AI,a))
A1b = lambdify([a,b,c],diff(AI,b))
A1c = lambdify([a,b,c],diff(AI,c))
A1aa = lambdify([a,b,c],diff(AI,a,a))
A1bb = lambdify([a,b,c],diff(AI,b,b))
A1cc = lambdify([a,b,c],diff(AI,c,c))
A1ab = lambdify([a,b,c],diff(AI,a,b))
A1bc = lambdify([a,b,c],diff(AI,b,c))
A1ca = lambdify([a,b,c],diff(AI,c,a))


# In[4]:


# 回転行列とその微分の値を算出する。
# 分子は動かないので、回転行列もその微分も定数である。



#for i in range(nmol):
#    for j in range(i):
#        print("distance {0}-{1}".format(i,j))
#        print(np.linalg.norm(pos[i][0:3]-pos[j][0:3]))

def dSite(pos, intra, intra2):

    def conv_euler(ea,eb,ec):
        sina=math.sin(ea)
        sinb=math.sin(eb)
        sinc=math.sin(ec)
        cosa=math.cos(ea)
        cosb=math.cos(eb)
        cosc=math.cos(ec)

        ap = math.acos(sina*sinc)
        cp = math.atan2(sina*cosc, cosa)
        bp = math.atan2(-sinb*sinc*cosa+cosb*cosc, -sinb*cosc-sinc*cosa*cosb)
        return ap,bp,cp
    
    nmol = pos.shape[0]
    # 回転行列に分子内座標をかけたもの。原子の重心からの空間座標
    D0 = np.zeros([nmol,3,4])
    # 回転行列の、角度による一階微分に、分子内座標をかけたもの。
    # D1 == dR_i / da_i
    D1 = np.zeros([nmol,3,3,4]) # molecule, euler angle, vectors, sites
    # 回転行列の、角度による二階微分に、分子内座標をかけたもの。
    # D2 == d^2 R_i / da db
    D2 = np.zeros([nmol,3,3,3,4]) # molecule, euler angle1, angle2, vectors, sites

    # Euler角が特異点に近付く場合には、分子内座標をとりなおし、それにあわせて角度も再計算する。
    lc = np.zeros(nmol, dtype=int)

    for i in range(nmol):
        ea,eb,ec = pos[i][3:6]
        # print("criteria", abs(math.sin(ea)))
        if abs(math.sin(ea)) > 0.0:
            D0[i] = A1(ea,eb,ec) @ intra
            D1[i,0] = A1a(ea,eb,ec) @ intra
            D1[i,1] = A1b(ea,eb,ec) @ intra
            D1[i,2] = A1c(ea,eb,ec) @ intra
            D2[i,0,0] = A1aa(ea,eb,ec) @ intra
            D2[i,1,1] = A1bb(ea,eb,ec) @ intra
            D2[i,2,2] = A1cc(ea,eb,ec) @ intra
            D2[i,0,1] = A1ab(ea,eb,ec) @ intra
            D2[i,1,2] = A1bc(ea,eb,ec) @ intra
            D2[i,2,0] = A1ca(ea,eb,ec) @ intra
            D2[i,1,0] = D2[i,0,1]
            D2[i,2,1] = D2[i,1,2]
            D2[i,0,2] = D2[i,2,0]
        else:
            ea,eb,ec = conv_euler(ea,eb,ec)
            pos[i][3:6] = ea,eb,ec
            lc[i] = 1 # rotation order of second kind
            D0[i] = A1(ea,eb,ec) @ intra2
            D1[i,0] = A1a(ea,eb,ec) @ intra2
            D1[i,1] = A1b(ea,eb,ec) @ intra2
            D1[i,2] = A1c(ea,eb,ec) @ intra2
            D2[i,0,0] = A1aa(ea,eb,ec) @ intra2
            D2[i,1,1] = A1bb(ea,eb,ec) @ intra2
            D2[i,2,2] = A1cc(ea,eb,ec) @ intra2
            D2[i,0,1] = A1ab(ea,eb,ec) @ intra2
            D2[i,1,2] = A1bc(ea,eb,ec) @ intra2
            D2[i,2,0] = A1ca(ea,eb,ec) @ intra2
            D2[i,1,0] = D2[i,0,1]
            D2[i,2,1] = D2[i,1,2]
            D2[i,0,2] = D2[i,2,0]
    return lc, D0, D1, D2


# In[5]:


lc, D0, D1, D2 = dSite(pos, intra, intra2)

# Tanaka programでの配列要素(1,2)は、Pythonでは[0,1]になる。
# 以後、縦横を精査。


# In[6]:


# 相互作用関数の距離での微分はsympyで自動生成する。

a1,a2,r = symbols('a1 a2 r')
F = 1/r
phic=lambdify([r], F)
phic_r=lambdify([r], diff(F,r))
phic_rr=lambdify([r], diff(F,r,r))
F = a1*r**(-12) - a2*r**(-6)
phiLJ=lambdify([r,a1,a2], F)
phiLJ_r=lambdify([r,a1,a2], diff(F,r))
phiLJ_rr=lambdify([r,a1,a2], diff(F,r,r))


# ## 切断関数
# 
# 相互作用は無限遠まで計算するのではなく、どこかで打ち切る。打ち切る距離のすこし手前から、スムーズに0に収束させるために、切断関数(Truncation function)をかける。
# 
# Truncation functionは、Site-site間ではなく重心間の距離の関数なので、切断関数を重心間距離やEuler角で微分するのは容易。ただし、相互作用関数にこれをかけた上で、微分すると、いろんな微分の組み合わせが発生して、かなり面倒な式になる。
# 

# In[7]:


rout = 8.5
trans = 2.0
rin  = rout - trans

x,y = symbols('x y')
# F0 is a truncation function ranged in [0,1]
F0 = -(x-1)**3*(6*x**2+3*x+1)
# 関数区間を変更
F = F0.subs(x,(y-rin)/trans)
# 解析的な微分を、数値関数に変換する

# Truncation function, raw
Trunc = lambdify([y], F)

# Truncation function, first derivative
Trunc_r = lambdify([y], diff(F,y))

# Truncation function, second derivative
Trunc_rr = lambdify([y], diff(F,y,y))

# Trunc(6.5), Trunc(8.5), Trunc_r(6.5), Trunc_r(8.5), Trunc_rr(6.5), Trunc_rr(8.5)


# In[20]:


# 相互作用と分子形状に関するブロック


# 当面、田中プログラムとの照合のために、エネルギーの単位は、
# 水素2原子が1 Aにいる時のCoulomb力を1とする。

# SI unit 2019
Na=6.02214076e23
ee=1.60217662e-19
E0=8.8541878128e-12
UJ = 4.184
qe0 = ee**2/(4*math.pi*E0)*Na/UJ*1e7
qeT = 332.17752e0
sw = (qeT*UJ*1e6/18.0)**0.5/(6.0*math.pi)
# Tanaka's Units
TanakaScale = qeT / qe0
TokJmol = ee**2/(4*math.pi*E0*1e-7)*Na*TanakaScale
print("TanakaScale", TanakaScale)

# 質点の質量
mass = np.array([wh,wh,wo,0.0])

Inertia = np.array([[mass @ np.sum(intra[1:3]**2, axis=0),
                     mass @ np.sum(intra[0:3:2]**2, axis=0),
                     mass @ np.sum(intra[0:2]**2, axis=0)],
                    [mass @ np.sum(intra2[1:3]**2, axis=0),
                     mass @ np.sum(intra2[0:3:2]**2, axis=0),
                     mass @ np.sum(intra2[0:2]**2, axis=0)]])

charge = [0.52,0.52,0,-1.04]
LJpairs = {(2, 2):(6.0e5*UJ / TokJmol, 6.10e2*UJ / TokJmol)}
CCpairs = dict()
for i,ci in enumerate(charge):
    for j,cj in enumerate(charge):
        cc = ci*cj
        if cc:
            CCpairs[i,j] = cc
            


# ## Hessian
# 
# Hessianは、全相互作用$E_p=\sum_{i,j}^N \phi$を座標変数(剛体の場合は重心位置とEuler角の6変数のベクトル)$\mathbf{q}$で二階微分したもの。
# 実際の相互作用は、重心間ではなく、サイト間(原子間)で定義されるので、微分はこみいったものになる。
# 
# 全相互作用が、分子間相互作用の総和で書ける場合は、ヘシアンも分子間に分解できる。
# 
# $$E_p=\sum_{i,j}^N \phi(i,j)$$
# $$H=\sum_{i,j}^N{\partial^2 \phi\over \partial\mathbf{q}\partial\mathbf{q'}}=\sum_{i,j}^N h(i,j)$$
# $$h(i,j) = \left(\begin{array}{cc}h_{tt}(i,j)&h_{tr}(i,j)\\ h_{rt}(i,j)&h_{rr}(i,j)\end{array}\right)$$
# 
# ### 並進-並進ヘシアン
# $h_{tt}(i,j)$は分子$i$と分子$j$の相互作用$\phi(i,j)$を、分子$i$または分子$j$の重心座標変数$\mathbf{r}=\{x_i,y_i,z_i,x_j,y_j,z_j\}$の2つで微分したものを要素とする、$6\times 6$行列である。
# 
# (下の関数`Hessian_tt`がこれにあたる。)
# 
# 
# $\phi(i,j)$は、さらに分子内のサイト対での相互作用の総和として書かれる。
# 
# $$\phi(i,j)=\sum_k\sum_l \phi_{kl}(r_{kl})$$
# 
# ただし$r_{kl}$はサイト-サイト距離である。
# 
# $$\mathbf{r}_{kl}=\mathbf{r}_{ij}+\mathbf{s}_{ik}-\mathbf{s}_{jl}$$
# 
# $\mathbf{r}_{ij}$は分子$i,j$の重心間ベクトル、$\mathbf{s}_{ik}$は分子$i$のサイト$k$の、重心からの相対位置ベクトルである。
# 
# さらに、分子内のサイト位置$\mathbf{s}_{ik}$は、分子内座標系でのサイト位置$\mathbf{w}_{k}$を回転行列$\mathbf{R}(a_i,b_i,c_i)$で回転したものである。
# 
# $$\mathbf{s}_{ik}=\mathbf{R}_i\cdot \mathbf{w}_k$$
# 
# 
# 

# #### 式展開
# 
# ##### 一階微分
# $\phi_{kl}$を$x_i$で微分する。サイト間($\phi_{kl}$の$r_{kl}$での微分)$\cdot$($r_{kl}$の$x_i$での微分)、という形で計算する。
# 
# ${d\phi(r)\over dr}$を`phi_r`のように書く。
# 
# $${\partial \phi\over \partial x_i}={\partial r\over\partial  x_i}\cdot{\partial \phi\over \partial r}$$
# $$={\partial \left|\mathbf{r}\right|\over \partial x_i}\cdot{\partial \phi\over \partial r}$$
# $$={r_x\over r}\cdot\phi_r$$

# In[21]:


### 説明のためのブロックです。実際の計算ではここを実行する必要はありません。

from sympy import *
# おまじない
init_printing()

# f(r)をベクトルrで微分する。
f_r, r_x, r_y, r_z = symbols('f_r, r_x, r_y, r_z')
r = Matrix(3,1,[r_x,r_y,r_z])
r_2 = r_x**2 + r_y**2 + r_z**2
D = diff(sqrt(r_2), r) * f_r
Di = D.subs((r_x**2+r_y**2+r_z**2)**0.5,'r')
Di


# In[22]:


# 距離に関する関数を、相対ベクトルで微分する。

def di(r, f_r):
    """
    粒子i,jの距離r = |r_i - r_j|の関数f(r)を、iの位置ベクトルr_iで一階微分する。
    f_r,は距離での一階微分
    答は3-ベクトル
    """
    return f_r*(r/np.linalg.norm(r))


# ##### 二階微分
# 
# $${\partial^2 \phi\over \partial x_i\partial x_j}={\partial\over\partial x_j}\left({\partial r\over\partial  x_i}\cdot{\partial \phi\over \partial r}\right)$$
# $$={\partial^2 r\over\partial  x_i\partial  x_j}\cdot{\partial \phi\over \partial r}+{\partial r\over\partial  x_i}\cdot\left({\partial\over\partial x_j}{\partial \phi\over \partial r}\right)$$
# $$={\partial^2 r\over\partial  x_i\partial  x_j}\cdot\phi_r-{r_x\over r}\cdot{r_x\over r}\phi_{rr}$$

# In[23]:


### 説明のためのブロックです。実際の計算ではここを実行する必要はありません。

# f(r)をベクトルrで2回微分する。
f_rr = symbols('f_rr')
Dii = Matrix(3,3,diff(sqrt(r_2), r.T, r) * f_r)          # first term
Dii += diff(sqrt(r_2), r) * diff(sqrt(r_2), r).T * f_rr  # second term
Dii = Dii.subs((r_x**2+r_y**2+r_z**2)**0.5,'r')          # simplify
Dii


# #### 実装
# 
# `phis`はサイトサイト間相互作用と、それの距離による多階微分を含む辞書であり、下の`derivatives()`関数で計算される、
# 
# また、切断関数$T$についてはここには詳しく書かないが、$\phi\cdot T$のHessianの計算には$T$の高階微分が含まれるので、かなり面倒。導出は`DerivTrunc.ipynb`を参照。

# In[24]:


# 距離に関する関数を、相対ベクトルで二階微分する。

# このあたりの関数も、できればsympyに生成させたい。
def dxixi(r, f_r, f_rr):
    # r is a vector
    rL2 = r@r
    rL = rL2**0.5
    return f_r*(1 - r**2 / rL2) / rL + f_rr*r**2 / rL2

def dxiyi(r, f_r, f_rr):
    # r is a vector
    rL2 = r@r
    rL = rL2**0.5
    yz = -f_r*r[1]*r[2]/(rL2*rL) + f_rr*r[1]*r[2]/rL2
    zx = -f_r*r[2]*r[0]/(rL2*rL) + f_rr*r[2]*r[0]/rL2
    xy = -f_r*r[0]*r[1]/(rL2*rL) + f_rr*r[0]*r[1]/rL2
    return yz,zx,xy

def dii(r, f_r, f_rr):
    """
    粒子i,jの距離r = |r_i - r_j|の関数f(r)を、iの位置ベクトルr_iで二階微分する。
    f_r, f_rrは距離での一階、二階微分
    答は3x3行列。
    """
    H = np.zeros((3,3))
    vxx,vyy,vzz = dxixi(r, f_r, f_rr)
    H[0,0] += vxx
    H[1,1] += vyy
    H[2,2] += vzz

    vyz,vzx,vxy = dxiyi(r, f_r, f_rr)
    H[0,1] += vxy
    H[1,2] += vyz
    H[0,2] += vzx
    H[1,0] += vxy
    H[2,1] += vyz
    H[2,0] += vzx
    return H


# In[25]:


# 剛体対ヘシアン


def Hessian_tt(phis, trunc):
    """
    Trans-Trans Hessian for a pair of molecules
    """
    H = np.zeros([3,3])
    for (si,sj),(r, (phi, phi_r, phi_rr)) in phis.items():
        h = dii(r, phi_r, phi_rr)
        #print("sdxx",si,sj,h[0,0]/0.52**2)
        #print("ri",1/np.linalg.norm(r))
        if trunc:
            #print("ramp",ramp, Tvv)
            v, vL, ramp, ramp_r, ramp_rr, Tv, Tvv = trunc 
            rL = np.linalg.norm(r)
            h = ramp*h + Tvv*phi
            h += ramp_r*phi_r/(vL*rL) * (np.outer(r,v)+np.outer(v,r))
            #print("OUTER",(np.outer(r,v)+np.outer(v,r))/(vL*rL))
        H += h
    Htt = np.zeros([2,3,2,3])
    Htt[0, :, 0, :] = H
    Htt[1, :, 1, :] = H
    Htt[0, :, 1, :] = -H
    Htt[1, :, 0, :] = -H
    # print(im,jm,Htt)
    return Htt


# ### 回転-回転ヘシアン
# 
# $h_{rr}(i,j)$は分子$i$と分子$j$の相互作用$\phi(i,j)$を、分子$i$または分子$j$の角度変数$\mathbf{\theta}=\{a_i,b_i,c_i,a_j,b_j,c_j\}$の2つで微分したものを要素とする、$6\times 6$行列である。
# 
# (下の関数`Hessian_rr`がこれにあたる。)
# 

# #### 式展開
# 
# ##### 一階微分
# 
# サイト間相互作用$\phi$を分子$i$のオイラー角$a_i$で微分する。
# 
# $${\partial\phi\over\partial a_i}={\partial\phi\over\partial r}\cdot{\partial r\over\partial a_i}$$
# 
# $r$はサイト間距離。$\mathbf{r}$をサイト間相対位置ベクトルとすると、
# 
# $${\partial r\over\partial a_i}={\partial r\over\partial\mathbf{r}}\cdot{\partial \mathbf{r}\over\partial a_i}={\mathbf{r}\over r}\cdot{\partial \mathbf{r}\over\partial a_i}$$
# 
# $${\partial \mathbf{r}\over\partial a_i}={\partial \mathbf{R}_i\over \partial a_i}\cdot \mathbf{w}_k$$
# 
# ${\partial \mathbf{R}_i\over \partial a_i}\cdot \mathbf{w}_k$を$\mathbf{s}_a$と書く。最終的に、
# 
# $${\partial\phi\over\partial a_i}=\phi_r\cdot {\mathbf{r}\over r}\cdot\mathbf{s}_a$$
# 
# 

# ##### 二階微分
# 
# サイト間相互作用$\phi$を分子$i$のオイラー角$a_i$と分子$j$のオイラー角$b_j$で微分する。
# 
# $${\partial^2\phi\over\partial a_i\partial b_j}={\partial\phi\over\partial r}\cdot{\partial\over\partial b_j}{\partial r\over\partial a_i}+{\partial\over\partial b_j}{\partial\phi\over\partial r}\cdot{\partial r\over\partial a_i}$$
# $$=\phi_r\cdot{\partial^2 r\over\partial a_i\partial b_j}+{\partial r\over\partial b_j}{\partial r\over\partial a_i}\phi_{rr}$$
# 
# ${\partial \mathbf{R}_i\over \partial b_j}\cdot \mathbf{w}_l$を$\mathbf{t}_b$と書くと、
# 
# $${\partial^2\phi\over\partial a_i\partial b_j}=\phi_r\cdot{\partial^2 r\over\partial a_i\partial b_j}+\left({\mathbf{r}\over r}\cdot \mathbf{s}_a\right)\left({\mathbf{r}\over r}\cdot \mathbf{t}_b\right)\phi_{rr}$$
# 

# #### 実装
# 
# 分子$i$のベクトル${\partial \mathbf{R}_i\over \partial a_i}\cdot \mathbf{w}_k=\mathbf{s}_a$は下のコードでは`D1i[p,q,r]`である。
# ただし、
# 
# * `p`は微分する$a_i$の軸(0〜2)
# * `q`はベクトルの軸要素
# * `r`はサイト番号(水なら0〜3)
# 
# 分子$i$のベクトル${\partial^2 \mathbf{R}_i\over \partial a_i\partial b_i}\cdot \mathbf{w}_k=\mathbf{s}_a$は下のコードでは`D2i[p1,p2,q,r]`である。
# 
# * `p`は微分する$a_i$,$b_i$の軸(0〜2)
# * `q`はベクトルの軸要素
# * `r`はサイト番号(水なら0〜3)
# 

# In[26]:


def Hessian_rr(phis, D1i, D2i, D1j, D2j, ramp):
    """
    Rot-Rot Hessian for a pair of molecules
    """
    Hrr = np.zeros([2,3,2,3])
    for (i,j),(r, (phi, phi_r, phi_rr)) in phis.items():
        rL = np.linalg.norm(r)

        # 式の導出はDerivative.ipynbを参照
        # **3 同じ分子の異なる角度変数の組みあわせ
        for d1 in range(3):
            for d2 in range(d1,3):
                h = (phi_r/rL**3*(rL**2 * (r @ D2i[d1,d2,:,i] + D1i[d1,:,i] @ D1i[d2,:,i]) -
                                                      (r @ D1i[d1,:,i])*(r @ D1i[d2,:,i])) +
                                         phi_rr/rL**2*(r @ D1i[d1,:,i])*(r @ D1i[d2,:,i]))
                Hrr[0,d1,0,d2] += h*ramp
                h = (phi_r/rL**3*(rL**2 * (-r @ D2j[d1,d2,:,j] + D1j[d1,:,j] @ D1j[d2,:,j]) -
                                                      (r @ D1j[d1,:,j])*(r @ D1j[d2,:,j])) +
                                         phi_rr/rL**2*(r @ D1j[d1,:,j])*(r @ D1j[d2,:,j]))
                Hrr[1,d1,1,d2] += h*ramp
        # **4 異なる分子の角度変数の組みあわせ
        for id in range(3):
            for jd in range(3):
                h = -(phi_r/rL**3*(rL**2 * (D1i[id,:,i] @ D1j[jd,:,j]) -
                                                       (r @ D1i[id,:,i])*(r @ D1j[jd,:,j])) +
                                          phi_rr/rL**2*(r @ D1i[id,:,i])*(r @ D1j[jd,:,j]))
                Hrr[0,id,1,jd] += h*ramp
        # 下半分はあとでコピーすればいい。
    # Hrrを一旦3nmol x 3nmolの平たい上半行列に戻し、
    Hrr = Hrr.reshape(6,6)
    # 対称化し、
    Hrr = Hrr + Hrr.T - np.diag(Hrr.diagonal())
    # また行列の配列に戻す。
    Hrr = Hrr.reshape(2,3,2,3)
    return Hrr


# ### 並進-回転ヘシアン
# 
# $h_{tr}(i,j)$は分子$i$と分子$j$の相互作用$\phi(i,j)$を、分子$i$または分子$j$の角度変数$\mathbf{\theta}=\{a_i,b_i,c_i\}$と重心座標の2つで微分したものを要素とする、$6\times 6$行列である。
# 
# (下の関数`Hessian_tr`がこれにあたる。)
# 
# #### 式展開
# 
# ##### 二階微分
# 
# $${\partial^2 \phi\over \partial x_i\partial a_j}={\partial^2 r\over\partial x_i\partial a_j}\cdot{\partial \phi\over \partial r}+{\partial r\over\partial  x_i}\cdot{\partial\over\partial a_j}{\partial \phi\over \partial r}$$
# $$={\partial^2 r\over\partial x_i\partial a_j}\cdot\phi_r+{\partial r\over\partial  x_i}\cdot{\partial r\over\partial a_j} \phi_{rr}$$
# $$={\partial^2 r\over\partial x_i\partial a_j}\cdot\phi_r+{r_x\over r}\cdot\left({\mathbf{r}\over r}\cdot \mathbf{t}_a\right)\phi_{rr}$$
# 
# $${\partial^2 r\over\partial x_i\partial a_i}={\partial \over \partial x_i}{\partial r\over \partial a_j}$$
# $$={\partial \over \partial x_i}{\mathbf{r}\over r}\mathbf{t}_a$$
# $$={\partial \over \partial x_i}\left({1\over r}(r_xt_{ax}+r_yt_{ay}+r_zt_{az})\right)$$
# $$={\partial r\over \partial x_i} {-1\over r^2} \mathbf{r}\cdot\mathbf{t}_a + {1\over r}\cdot t_{ax}$$
# $$={-r_x\over r^3} \mathbf{r}\cdot\mathbf{t}_a + {t_{ax}\over r} $$

# #### 実装
# 

# In[27]:


def Hessian_tr(phis, D1i, D1j, ramp, Tv):
    """
    Trans-Rot Hessian for a pair of molecules
    """
    Htr = np.zeros([2,3,2,3])
    for (i,j),(r, (phi, phi_r, phi_rr)) in phis.items():
        rL = np.linalg.norm(r)

        # **5 同じ分子
        for it in range(3):
            for ir in range(3):
                h = (phi_r/rL**3 * (rL**2*D1i[ir,it,i] - r[it]*(r @ D1i[ir,:,i])) +
                                          phi_rr/rL**2 * r[it]*(r @ D1i[ir,:,i]))
                if ramp != 1.0:
                    h = ramp*h + Tv[it]*phi_r/rL*(r @ D1i[ir,:,i])
                Htr[0,it,0,ir] += h
        for jt in range(3):
            for jr in range(3):
                h = (phi_r/rL**3 * (rL**2*D1j[jr,jt,j] - r[jt]*(r @ D1j[jr,:,j])) +
                                          phi_rr/rL**2 * r[jt]*(r @ D1j[jr,:,j]))
                if ramp != 1.0:
                    h = ramp*h + Tv[jt]*phi_r/rL*(r @ D1j[jr,:,j])
                Htr[1,jt,1,jr] += h
        # **6 異なる分子
        for it in range(3):
            for jr in range(3):
                h = (phi_r/rL**3 * (rL**2*D1j[jr,it,j] - r[it]*(r @ D1j[jr,:,j])) +
                                          phi_rr/rL**2 * r[it]*(r @ D1j[jr,:,j]))
                if ramp != 1.0:
                    h = ramp*h + Tv[it]*phi_r/rL*(r @ D1j[jr,:,j])
                Htr[0,it,1,jr] -= h
        for jt in range(3):
            for ir in range(3):
                h = (phi_r/rL**3 * (rL**2*D1i[ir,jt,i] + r[jt]*(-r @ D1i[ir,:,i])) +
                                          phi_rr/rL**2 * -r[jt]*(-r @ D1i[ir,:,i]))
                if ramp != 1.0:
                    h = ramp*h + Tv[jt]*phi_r/rL*(r @ D1i[ir,:,i])
                Htr[1,jt,0,ir] -= h
    return Htr


def E_pair(phis, ramp):
    ep = 0.0
    for (si,sj),(r, (phi, phi_r, phi_rr)) in phis.items():
        ep += ramp*phi
    return ep


# In[28]:


# 相互作用をサイト間距離で微分したもの。

def derivatives(CCpairs, LJpairs, v, D0i, D0j):
    phis = dict()

    for (i,j),cc in CCpairs.items():
        r = D0i[:,i] - D0j[:,j] + v
        if (i,j) not in phis:
            phis[i,j] = [r, np.zeros(3)]
        rL = np.linalg.norm(r)
        phi    = cc*phic(rL)
        phi_r  = cc*phic_r(rL)
        phi_rr = cc*phic_rr(rL)
        phis[i,j][1] += np.array([phi, phi_r, phi_rr])

    for (i,j),(a1,a2) in LJpairs.items():
        r = D0i[:,i] - D0j[:,j] + v
        if (i,j) not in phis:
            phis[i,j] = [r, np.zeros(3)]
        rL = np.linalg.norm(r)
        phi    = phiLJ(rL,a1,a2)
        phi_r  = phiLJ_r(rL,a1,a2)
        phi_rr = phiLJ_rr(rL,a1,a2)
        phis[i,j][1] += np.array([phi, phi_r, phi_rr])

    return phis


# In[29]:


## Hessianの計算。もっと小さい関数に分けたいが、かなり手強い。

def Hessian(pos, cell, D0, D1, D2, CCpairs, LJpairs, Trunc, Trunc_r, Trunc_rr, rin, rout):
    # 3x3行列の配列。
    Htt = np.zeros([nmol, 3, nmol, 3])
    Hrr = np.zeros([nmol, 3, nmol, 3])
    Htr = np.zeros([nmol, 3, nmol, 3])
    ep = 0
    
    celli = np.linalg.inv(cell)
    for im in range(nmol):
        for jm in range(im+1,nmol):
            # 重心間相対ベクトル
            v = pos[im][:3] - pos[jm][:3]
            # 周期境界条件
            v -= cell @ np.floor( celli @ v + 0.5 )
            # 重心間距離
            vL = np.linalg.norm(v)
            # print("vL",vL)

            # 切断関数を準備。
            trunc = False
            if vL > rout:
                # 遠すぎる対はパス
                continue
            if vL > rin:
                # カットオフ関数がかかる場合
                ramp = Trunc(vL)
                ramp_r = Trunc_r(vL)
                ramp_rr = Trunc_rr(vL)
                # Tのvによる一階微分(ベクトル)
                Tv = di(v, ramp_r)
                # Truncのvによる二階微分(テンソル)
                Tvv = dii(v, ramp_r, ramp_rr)
                trunc = v, vL, ramp, ramp_r, ramp_rr, Tv, Tvv
            else:
                # 近距離
                ramp = 1.0
                Tv = None

            # サイト間の相互作用関数の、相対位置ベクトルによる多階微分をあらかじめ配列phisに入れておく。
            phis = derivatives(CCpairs, LJpairs, v, D0[im], D0[jm])

            # Check EP again
            ep += E_pair(phis, ramp)

            # 並進・並進
            Htt0 = Hessian_tt(phis, trunc)
            Htt[im, :, im, :] += Htt0[0, : , 0, :]
            Htt[jm, :, jm, :] += Htt0[1, : , 1, :]
            Htt[im, :, jm, :] += Htt0[0, : , 1, :]
            Htt[jm, :, im, :] += Htt0[1, : , 0, :]

            # 次は回転・回転
            Hrr0 = Hessian_rr(phis, D1[im], D2[im], D1[jm], D2[jm], ramp)
            Hrr[im, :, im, :] += Hrr0[0, : , 0, :]
            Hrr[jm, :, jm, :] += Hrr0[1, : , 1, :]
            Hrr[im, :, jm, :] += Hrr0[0, : , 1, :]
            Hrr[jm, :, im, :] += Hrr0[1, : , 0, :]

            # 最後に並進・回転
            Htr0 = Hessian_tr(phis, D1[im], D1[jm], ramp, Tv)
            Htr[im, :, im, :] += Htr0[0, : , 0, :]
            Htr[jm, :, jm, :] += Htr0[1, : , 1, :]
            Htr[im, :, jm, :] += Htr0[0, : , 1, :]
            Htr[jm, :, im, :] += Htr0[1, : , 0, :]
            
    return Htt, Hrr, Htr, ep


# In[30]:


get_ipython().run_cell_magic('time', '', '\nHtt, Hrr, Htr, ep = Hessian(pos, cell, D0, D1, D2, CCpairs, LJpairs, Trunc, Trunc_r, Trunc_rr, rin, rout)\n\nprint("energy (kJ/mol)=",ep*TokJmol/nmol)\n\n# 系が小さい場合には、圧倒的にここが遅い。')


# In[32]:


get_ipython().run_cell_magic('time', '', '# np.set_printoptions(4)\n\n\n\n# ***1 謎のテンソルの計算\neival = np.zeros([nmol,3])\neivec = np.zeros([nmol,3,3])\n\nfor i in range(nmol):\n    # 分子内座標での慣性テンソル\n    Ii = Inertia[lc[i]] / wm\n# wmで割ることで、並進には質量を掛ける必要がなくなる。\n# verified rix, riy,riz\n\n    ea,eb,ec = pos[i][3:6]\n    sina = math.sin(ea)\n    cosa = math.cos(ea)\n    sinb = math.sin(eb)\n    cosb = math.cos(eb)\n    sinc = math.sin(ec)\n    cosc = math.cos(ec)\n    ev = np.array([[Ii[0]*cosc**2 + Ii[1]*sinc**2, (Ii[0]-Ii[1])*sina*sinc*cosc, 0.0],\n                   [(Ii[0]-Ii[1])*sina*sinc*cosc,  (Ii[0]*sinc**2+Ii[1]*cosc**2)*sina**2+Ii[2]*cosa**2, Ii[2]*cosa],\n                   [0.0,  Ii[2]*cosa, Ii[2]]])\n    #print("***6")\n    #print(ev)\n    # ev verified\n    ival, ivec = np.linalg.eig(ev)\n    idx = np.argsort(-ival)\n    ival = ival[idx]\n    ivec = ivec[:,idx]\n    eival[i] = 1/ival**0.5\n    eivec[i] = ivec\n    # eig verified\n    #print("eival")\n    #print(eival[i])\n    #print("eivec")\n    #print(eivec[i])\n\n# rot-rot\nfor i in range(nmol):\n    for j in range(nmol):\n        # elements of rot-rot\n        t = Hrr[j,:,i,:]\n        #print("***4",i,j)\n        #print(t)\n        t = eivec[j].T @ t @ eivec[i]\n        sd = t * np.outer(eival[j],eival[i])\n        Hrr[j,:,i,:] = sd\n        #print("***5",i,j)\n        #print(sd)\n        #print()\n\n# rot-trans\nfor i in range(nmol):\n    for j in range(nmol):\n        # elements of trans-rot\n        t = Htr[i,:,j,:].T\n        #print("***7")\n        #print(t)\n        t = Htr[i,:,j,:]\n        t = t @ eivec[j]\n        sd = t * eival[j]\n        Htr[i,:,j,:] = sd\n        #print("***8")\n        #print(sd)\n        #print()\n\nHessian = np.zeros([6*nmol, 6*nmol])\nHessian[0*nmol:3*nmol,0*nmol:3*nmol] = Htt.reshape(3*nmol,3*nmol)\nHessian[3*nmol:6*nmol,3*nmol:6*nmol] = Hrr.reshape(3*nmol,3*nmol)\nHessian[0*nmol:3*nmol,3*nmol:6*nmol] = Htr.reshape(3*nmol,3*nmol)\nHessian[3*nmol:6*nmol,0*nmol:3*nmol] = Htr.reshape(3*nmol,3*nmol).T\n\nnp.set_printoptions(4)\n\n#print(Hessian/0.52**2)')


# In[33]:


get_ipython().run_cell_magic('time', '', 'ival, ivec = np.linalg.eigh(Hessian)\n\nfor iv in sorted(ival,reverse=True):\n    if iv > 0:\n        iv = iv**0.5 * sw\n    else:\n        iv = -(-iv)**0.5 * sw\n    print(iv)')

95.18011714359442
44.2507979270288
27.18980399678063
9.887480827730782e-07
5.79563716924399e-07
-6.564938741762731e-07
-3.2567213126119685
-7.482131380620787
-28.6285571437099
-37.195860241136856
-50.656183591544234
-72.48723707456553
CPU times: user 630 µs, sys: 235 µs, total: 865 µs
Wall time: 658 µs
# ## 完了
# * 第二回転規準
#    * 全部第二規準にすると、全部第二規準とした田中プログラムと同じ答えになることを確認。
#    * 一方で、田中プログラム自体、規準を変更すると振動数が変わるのがとても気になる。
#    * 回転行列を二種類準備する代わりに、分子内座標を2種類準備することにした。
#      * xyz→yzxに変えるようにした。
#      * 角度変換規則(euler_conv)は田中のものと結果的には同じであった。たぶん発想は同じ。
#      * その場合でも、振動数が微妙に変化する。変化のしかたも同様なので、同じことをやっているにすぎないと思われる。
#      * ということは、分子内座標で、分子をどういう向きに置くかによって、結果が変わるということ。
#      * Hessianの計算に関しては、あまり疑う余地がないので、mass weightのあたりしか考えられないのだが…。
#      * Quenchしてあれば問題がおきないということで合意。Quenchしない場合に変になる理由はわからない。
# * LJ
# * 多分子
# * 周期境界条件
# * Truncation
# * 特殊な単位系の廃止
# * 小さい系(CS2)であれば、現実的な速度で正確な振動を得ることができた。2019-09-26
# * 整理し、読みやすくする。
# 
# ## 未完
# 
# * Hessianの式の自動生成
# * 多成分
# * 直線分子
# * 単原子分子
# 
# ## 別計画
# 
# * 数値微分によるHessian計算はどれぐらい現実味があるのか。
#    * このプログラムがうまく動くことが確認できれば、比較対象に使える。 
# 

# In[ ]:


intra


# In[ ]:


intra[1:3] * intra[1:3]


# In[ ]:


np.sum(intra[1:3]**2, axis=0) @ mass


# In[ ]:


np.outer(np.array([1,2,3]),np.array([1,2,3]))


# In[ ]:


np.diag(np.array([1,1,1]))


# In[121]:


np.cross(np.array([1,2,3]),np.array([1,2,3]))


# In[ ]:




