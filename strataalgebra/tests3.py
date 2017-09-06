from StrataAlgebra import *
from sage.all import QQ, var

#s = StrataAlgebra(QQ,1,(1,2,3))

#a = s.psi(2) * s.irr()

s = StrataAlgebra(QQ,1,(1,1))
print s.psi(1)*s.psi(1)


