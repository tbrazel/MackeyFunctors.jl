using AbstractAlgebra
M = free_module(ZZ,1)
id_hom = hom(M,M,ZZ[1;])
twoZ, _ = image(ZZ(2)*id_hom)
threeZ, _ = image(ZZ(3)*id_hom)

Zmodtwo , _ = quo(M,twoZ)
Zmodthree , _ = quo(M,threeZ)
D,_,_ = direct_sum(Zmodtwo,Zmodthree)
invariant_factors(D)