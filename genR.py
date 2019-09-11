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

A=B*C*D
R=A.T

from sympy.utilities.codegen import codegen
result = codegen([('R', R),
                  ('Ra', diff(R,a)),
                  ('Rb', diff(R,b)),
                  ('Rc', diff(R,c)),
                  ('Raa', diff(R,a,a)),
                  ('Rab', diff(R,a,b)),
                  ('Rac', diff(R,a,c)),
                  ('Rbb', diff(R,b,b)),
                  ('Rbc', diff(R,b,c)),
                  ('Rcc', diff(R,c,c))],
                 language='f95',
                 project='NMA',
                 to_files=True)
