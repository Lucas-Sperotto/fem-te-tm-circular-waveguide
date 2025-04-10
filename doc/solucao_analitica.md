# Solução Analítica dos Modos TE e TM em Guias de Onda Circulares

A solução analítica dos modos em um guia de onda circular parte da equação de Helmholtz, resultante da separação das equações de Maxwell no vácuo. Em coordenadas cilíndricas, os campos elétrico e magnético assumem formas que envolvem funções de Bessel.

---

## Modos TM (Transversais Magnéticos)

Para os modos TM, a componente longitudinal do campo elétrico $E_z$ é diferente de zero, e o campo magnético não possui componente na direção axial ($H_z = 0$).

A solução geral para $E_z$ é:

$$
E_z(r, \phi) = J_n(k_c r) \cos(n\phi) \quad 	\text{ou} \quad J_n(k_c r) \sin(n\phi)
$$

A condição de contorno de campo elétrico nulo na parede metálica exige:

$$
J_n(k_c a) = 0
$$

---

## Modos TE (Transversais Elétricos)

Nos modos TE, o campo elétrico longitudinal é nulo ($E_z = 0$) e o campo magnético $H_z$ possui a forma:

$$
H_z(r, \phi) = J_n(k_c r) \cos(n\phi) \quad 	\text{ou} \quad J_n(k_c r) \sin(n\phi)
$$

A condição de contorno, obtida da continuidade da componente tangencial do campo magnético, leva a:

$$
J_n'(k_c a) = 0
$$

---

## Frequência de Corte

O número de onda transversal $k_c$ é relacionado à frequência de corte $f_c$ por:

$$
f_c = \frac{c}{2\pi} k_c = \frac{c}{2\pi} \frac{x_{nm}}{a}
$$

onde:
- $x_{nm}$ é a m-ésima raiz da função de Bessel $J_n$ (TM) ou da derivada $J_n'$ (TE),
- $a$ é o raio do guia de onda,
- $c$ é a velocidade da luz no vácuo.

Essas expressões fornecem os valores teóricos que podem ser comparados com os obtidos numericamente por elementos finitos.

