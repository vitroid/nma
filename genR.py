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
result = codegen([('RotMat', R),
                  ('Rot_a', diff(R,a)),
                  ('Rot_b', diff(R,b)),
                  ('Rot_c', diff(R,c)),
                  ('Rot_aa', diff(R,a,a)),
                  ('Rot_ab', diff(R,a,b)),
                  ('Rot_ac', diff(R,a,c)),
                  ('Rot_bb', diff(R,b,b)),
                  ('Rot_bc', diff(R,b,c)),
                  ('Rot_cc', diff(R,c,c))],
                 prefix='Rotation',
                 language='f95',
                 project='NMA',
                 to_files=True)
