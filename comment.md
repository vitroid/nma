# 感想・展望など

## 簡潔にするために

### drvtv

回転行列は、軸まわりの回転行列3つを掛けて作れる。回転行列の生成と、それの角度微分の生成は、いずれも記号処理を使って自動生成できる。

MathematicaはFortranコードを吐ける。

SymPyはFortranコードを吐けるか? [もちろん!](https://docs.sympy.org/latest/modules/utilities/codegen.html)

こんな感じでSympyで行列を生成させてみよう。(コードは`Rotation.ipynb`に記載)

```python
from sympy import *

θ, ψ, φ=symbols('θ ψ φ')

B = Matrix([[ cos(ψ),sin(ψ),0],
            [-sin(ψ),cos(ψ),0],
            [0,0,1]])
C = Matrix([[1,0,0],
            [0, cos(θ),sin(θ)],
            [0,-sin(θ),cos(θ)]])
D = Matrix([[ cos(φ),sin(φ),0],
            [-sin(φ),cos(φ),0],
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

