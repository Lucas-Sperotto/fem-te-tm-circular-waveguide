# Fundamentos Teóricos dos Modos TE e TM em Guias de Onda Circulares

Guias de onda circulares são estruturas metálicas cilíndricas utilizadas para conduzir ondas eletromagnéticas. A solução das equações de Maxwell nesse tipo de guia leva a dois tipos principais de modos propagantes:

- **Modos Transversais Elétricos (TE)**: o campo elétrico não possui componente ao longo do eixo do guia ($E_z = 0$).
- **Modos Transversais Magnéticos (TM)**: o campo magnético não possui componente axial ($H_z = 0$).

A separação de variáveis e a simetria cilíndrica permitem resolver as equações de Maxwell usando funções de Bessel.

As frequências de corte e os campos associados dependem das **raízes das funções de Bessel**:
- Modos **TE**: raízes de $J'_n(k_c r) = 0$.
- Modos **TM**: raízes de $J_n(k_c r) = 0$.

Essas raízes determinam os valores próprios do número de onda transversal $k_c$ os quais se relacionam diretamente com as frequências de corte por:

$$
f_c = \frac{c}{2\pi} k_c = \frac{c}{2\pi} \frac{x_{n,m}}{a}
$$

onde:
- $c$ é a velocidade da luz no vácuo,
- $a$ é o raio do guia de onda,
- \$x_{n,m}$ é a m-ésima raiz da função de Bessel $J_n$ ou sua derivada, conforme o modo.

