# 感想・展望など

## 解読

### 第二回転規準

lc2の行列の構造から、回転順序を推定する。

lc2のsp??変数は回転行列の要素を示しているので、これを抽出して並べてみる。

```
? coscp*cosbp-cosap*sinbp*sincp -(sincp*cosbp+cosap*sinbp*coscp)
? coscp*sinbp+cosap*cosbp*sincp -sincp*sinbp+cosap*cosbp*coscp
? sinap*sincp sinap*coscp

```

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

つまり、最初にx軸、次にy軸、最後にz軸で回転する。しかも、角度の基準が変。

lc1のほうは、下に書いてある通りで、Goldsteinの定義のまま。

ただ、2つの行列を見比べれば、角度を読み替える方法はわかる。

1. 要素(3,1)を比較し、`cos(ap) = sin(a)*sin(c)`から`ap`を求める。
2. 要素(3,3)を比較し、`cos(a) = sin(ap)*cos(cp)`から`cp`を求める。
3. 要素(1,1)を比較し、`sin(ap)*sin(bp) = -sin(b)*sin(c)*cos(a) + cos(b)*cos(c)`から`bp`を求める。

という手順で`(ap,bp,cp)`を求めている。この手続きはこのまま採用するしかないだろう。

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

