# 解読

## 第二回転規準

lc2の行列の構造から、回転順序を推定する。

lc2のsp??変数は回転行列の要素を示しているので、これを抽出して並べてみる。

```
? coscp*cosbp-cosap*sinbp*sincp -(sincp*cosbp+cosap*sinbp*coscp)
? coscp*sinbp+cosap*cosbp*sincp -sincp*sinbp+cosap*cosbp*coscp
? sinap*sincp sinap*coscp

```

?の部分は、水の場合には必要がないので、計算すらしていない。

これを睨みながら、これにあう回転行列の構成を試行錯誤する。

```
ap, cp, bp=symbols('ap cp bp')
E = Matrix([[1,0,0],
            [0, cos(cp),-sin(cp)],
            [0, sin(cp),cos(cp)]])
F = Matrix([[-sin(ap),0,cos(ap)],
            [0, 1, 0],
            [cos(ap), 0,sin(ap)]])
G = Matrix([[-sin(bp),cos(bp),0],
            [cos(bp),sin(bp),0],
            [0,0,1]])

G*F*E
```

つまり、最初にx軸、次にy軸、最後にz軸で回転する。しかも、角度の基準が変。FとGは角度0でも単位行列にならないのでおかしい。でも回転順序が問題にはならないので、このまま進める。

lc1のほうは、下に書いてある通りで、Goldsteinの定義のまま。

ただ、2つの行列を見比べれば、角度を読み替える方法はわかる。

1. 要素(3,1)を比較し、`cos(ap) = sin(a)*sin(c)`から`ap`を求める。
2. 要素(3,3)を比較し、`cos(a) = sin(ap)*cos(cp)`から`cp`を求める。
3. 要素(1,1)を比較し、`sin(ap)*sin(bp) = -sin(b)*sin(c)*cos(a) + cos(b)*cos(c)`から`bp`を求める。

という手順で`(ap,bp,cp)`を求めている。この手続きはこのまま採用するしかないだろう。

## 相互作用の一階微分

相互作用の記号演算を考える。

回転行列$\mathbf{R}_i$は角度$\mathbf{\theta}_i=(a_i,b_i,c_i)$でできている。

原子$j$の分子内座標を$\mathbf{w}_{ij}$、分子の重心位置を$\mathbf{v}_i=(x_i,y_i,z_i)^t$とすると、原子$j$の空間位置$\mathbf{r}_{ij}$は$\mathbf{r}_{ij}=\mathbf{v}_i+\mathbf{R}_i\cdot \mathbf{w}_{ij}$と書ける。

さらに相互作用は$\mathbf{r}_{ij}$の関数として書かれる。こうして書かれた相互作用は、$\mathbf{v}_i$や$\mathbf{R}_i$の各要素で微分することができる。

ポテンシャルエネルギーの関数を$\phi$としよう。$i$分子と、$k$分子の相互作用はこんな感じで書ける。

$$U_{ik}=\sum_j\sum_l \phi_{jl}(\mathbf{r}_{ij}-\mathbf{r}_{kl})$$

### 並進
分子の並進座標$\mathbf{r}_i$による一階微分は力。原子間距離を$\mathbf{r}_{ij}-\mathbf{r}_{kl}=\mathbf{r}$、$|\mathbf{r}|=r$と略記すると、

$${\partial \phi_{jl}\over \partial \mathbf{r}_i}={\partial r \over \partial \mathbf{r}_i}\cdot{\partial \phi_{jl}\over \partial r}$$
$${\partial r \over \partial \mathbf{r}_i}={\partial \left|\mathbf{r}\right|^2 \over \partial \mathbf{r}_i}\cdot{\partial \sqrt{r^2} \over \partial (r^2)}$$
$${\partial \left|\mathbf{r}\right|^2 \over \partial \mathbf{r}_i}=2\mathbf{r}$$
$${\partial \sqrt{r^2} \over \partial (r^2)}={1\over 2r}$$

つまり、
$${\partial r \over \partial \mathbf{r}_i}={\mathbf{r} \over r}$$

これらをあわせると、
$${\partial \phi_{jl}\over \partial \mathbf{r}_i}={\mathbf{r} \over r}\cdot{\partial \phi_{jl}(r)\over \partial r}$$

### 回転

次に角度について。オイラー角三つをまとめて$\mathbf{\theta}_i$と書く。分子の回転角$\mathbf{\theta}_i$による一階微分はトルク。

$${\partial \phi_{jl}\over \partial \mathbf{\theta}_i}={\partial r \over \partial \mathbf{\theta}_i}\cdot{\partial \phi_{jl}\over \partial r}$$
$${\partial r \over \partial \mathbf{\theta}_i}={\partial \left|\mathbf{r}\right|^2 \over \partial \mathbf{\theta}_i}\cdot{\partial \sqrt{r^2} \over \partial (r^2)}$$
$${\partial r^2 \over \partial \mathbf{\theta}_i}={\partial \over \partial \mathbf{\theta}_i}\left(\mathbf{v}_i+\mathbf{R}_i\cdot \mathbf{w}_j-\mathbf{v}_k-\mathbf{R}_k\cdot \mathbf{w}_l\right)^2$$
$$=2\mathbf{v}_i{\partial \mathbf{R}_i\over \partial \mathbf{\theta}_i}\cdot \mathbf{w}_j+2{\partial \mathbf{R}_i\over \partial \mathbf{\theta}_i}\cdot \mathbf{w}_j\cdot \mathbf{R}_i\cdot \mathbf{w}_j-2{\partial \mathbf{R}_i\over \partial \mathbf{\theta}_i}\cdot \mathbf{w}_j\cdot(\mathbf{v}_k+\mathbf{R}_k\cdot \mathbf{w}_l)$$
$$=2\left(\mathbf{v}_i+\mathbf{R}_i\cdot \mathbf{w}_j-\mathbf{v}_k-\mathbf{R}_k\cdot \mathbf{w}_l\right)\cdot{\partial \mathbf{R}_i\over \partial \mathbf{\theta}_i}\cdot \mathbf{w}_j$$
$$=2\mathbf{r}\cdot{\partial \mathbf{R}_i\over \partial \mathbf{\theta}_i}\cdot \mathbf{w}_j$$

つまり、
$${\partial r \over \partial \mathbf{\theta}_i}={ \mathbf{r} \over  r}\cdot{\partial \mathbf{R}_i\over \partial \mathbf{\theta}_i}\cdot \mathbf{w}_j$$

これらをあわせると、
$${\partial \phi_{jl}\over \partial \mathbf{\theta}_i}={ \mathbf{r} \over  r}\cdot{\partial \mathbf{R}_i\over \partial \mathbf{\theta}_i}\cdot \mathbf{w}_j{\partial \phi_{jl}(r)\over \partial r}$$

回転行列のオイラー角による微分はRotation.ipynbですでに求めた。

## 相互作用の二階微分

Hessianの計算は面倒だが、ここまでに使った式を流用できるはず。

### 並進・並進

#### 同一変数の場合
ベクトルでの微分で表記するとややこしいので、成分表記にする。

$${\partial^2 \phi_{jl}(r)\over \partial x_i^2}=\
{\partial\over\partial x_i}\left({r_x \over r}\cdot{\partial \phi_{jl}(r)\over \partial r}\right)$$

ただし$r_x$は$\mathbf{r}$の$x$成分。

$${\partial r_x\over\partial x_i}=1$$
$${\partial r\over\partial x_i}={x_i\over r}$$
$${\partial\over\partial x_i}\left({\partial \phi_{jl}(r)\over \partial r}\right)=\
{\partial r\over\partial x_i}\left({\partial^2 \phi_{jl}(r)\over \partial r^2}\right)$$
$$={x_i\over r}\left({\partial^2 \phi_{jl}(r)\over \partial r^2}\right)$$
これらを組みあわせると、
$${\partial^2 \phi_{jl}(r)\over \partial x_i^2}={1\over r^2}\left(r+{r_x x_i\over r}\right){\partial \phi_{jl}(r)\over \partial r}+{r_x x_i\over r^2}{\partial^2 \phi_{jl}(r)\over \partial r^2}$$

にわかには信じられない。

`Sympy`でこの展開をできないか。

#### 異なる変数の場合


# 感想・展望など

集中すれば一週間ぐらいで移植できそうだ。これができれば、適用できる分子や混合物の範囲がひろがる。




## 簡潔にするために

### drvtv

回転行列は、軸まわりの回転行列3つを掛けて作れる。回転行列の生成と、それの角度微分の生成は、いずれも記号処理を使って自動生成できる。

MathematicaはFortranコードを吐ける。

SymPyはFortranコードを吐けるか? [もちろん!](https://docs.sympy.org/latest/modules/utilities/codegen.html)

こんな感じでSympyで行列を生成させてみよう。(コードは`Rotation.ipynb`に記載)

```python
from sympy import *

a, b, c=symbols('a b c')

B = Matrix([[ cos(c),sin(c),0],
            [-sin(c),cos(c),0],
            [0,0,1]])
C = Matrix([[1,0,0],
            [0, cos(a),sin(a)],
            [0,-sin(a),cos(a)]])
D = Matrix([[ cos(b),sin(b),0],
            [-sin(b),cos(b),0],
            [0,0,1]])

# 
A=B*C*D
A.T            
```

```
diff(A.T, φ, θ)
```

### 行列と配列に関する注意

Fortranでは、配列の添え字は数学表記(行、列)と同じでよかった、はず。
matmulによるかけ算もその順番でやってくれるのかどうかを確認しておく。

test_matmul.f90で見る限りそれで問題ない。

## Pythonで書けるか?

固有値計算が最大のボトルネックになると予想されるが、そこはNumpyを使えばFortanなみ(GPUなどが使えるようになったらそれ以上)に速いので、ほかの手間を考えれば、Fortranで書く必然性はないと言えよう。(いざとなったらそこだけFortranにすればいい。)

