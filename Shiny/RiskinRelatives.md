---
title: "test"
author: "Bart Baselmans"
date: "21/08/2020"
output: html_document
---
## Liability Threshold Model

When many factors contribute to the risk of disease it can be helpful to consider a model of disease where there is a latent distribution of liability to disease underlying the dichotomous Case/Control status.

Since this latent liability comprises many genetic and other risks it is reasonable to assume that the liability distribution is approximately normally distributed (since many things added together will make a bell-shaped distribution in a population sample) and that those affected by disease have the combination of genetic and other factors that it places them in the top end of the liability distribution. 

Hence, those whose  liability to disease is above the threshold are affected, and so this model is sometimes called the liability threshold model. The threshold can be calculated using the lifetime risk ($K$) parameter as $K$ represents the proportion of individuals in the general population have the disorder of interest. The threshold (or boundary) that determines the area of $K$ (indicated in blue) in a normal distribution is given by:
$$
T=\Phi^{-1}(1-K)
$$

Conversely, the life-time risk ($K$) is derived from the threshold parameter
$$
K = 1-\Phi(T)
$$

where $\Phi^{-1}(x)$ is the inverse of the cumulative standard normal distribution function. 


## Risk in relatives and heritability of liability

The lifetime risk of a disease in relatives ($K_{R}$) of those affected by the disease is expected to be greater, or equal to, the lifetime risk of the disease in a population sample. It is logical to assume that the threshold in liability associated with disease has the same value in the relatives of those affected. As a result, the liability distribution in first degree members must be shifted (in the direction of increased liability) compared to the general population to an extent consistent with the observed higher risk of $K_{R}$.

We can then ask: what proportion of the variation in liability must be attributable to genetic factors (i.e., what is the heritability?) in order to generate this observed increased risk in relatives, given the known coefficient of relationship between the relatives, $a_R$, i.e., 0.5 for 1st degree relatives.

The difference in mean liabilities between population ($M_{pop}$ = 0) and relatives ($M_{R}$), can be shown to be equivalent to the difference in thresholds when the liability distribution of relatives is scaled back to a N(0,1) distribution. i.e., $M_{R}$ - $M_{pop}$ = $a_Rih^2$ = $T_R - T$, where   $T_R=\Phi^{-1}(1-K_R)$.

However, Falconer's derivation assumed that the genetic (and hence liability) variance amongst the relatives was not changed by ascertainment on disease in the probands. [Reich et al.,(1972)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1469-1809.1972.tb00767.x) showed that while it is OK to assume that the distribution of disease liability in relatives of those affected is approximately normal, it should be recognised that the variance in liability is slightly reduced in the relatives as a result of ascertaining the proband as been affected. They showed that the variance in liability is reduced by a factor of $a_{R}^2h^4i(i-T)$ 
so that the scaling of the $T_R$ to a standard N(0,1) distribution has a denominator of the standard deviation of liability in the relatives given that the probands have disease:

$$
T_{R} = \frac{T -a_Rih^2}{\sqrt{1-a_{R}^2h^4i(i-T)}}
$$

Since heritability is in the numerator and denominator, making $h^2$ the subject of the equation requires solving of a quadratic equation to give:

$$
h^2 = \frac{T - T_R \sqrt{1 -(1-\frac{T}{i})(T^2-T^2_R)}}{a_R(i+(i-T)T^2_R)}
$$

Although, this equation looks complicated, in fact everything on the right hand side again can be calculated from two observations $K$ and $K_R$, and $a_R$ (and holds as long as $a_R$ < 1), and  if ascertainment is based on the disease status of one proband. A special case is when ascertainment is based on both parents being affected, hence the $K$ and $K_{2PAR}$ is estimated, with risk ratio $RR_{2PAR} = \frac{K_{2PAR}}{K}$

In this case it can be shown in [Wray & Gottesman (2012)](https://www.frontiersin.org/articles/10.3389/fgene.2012.00118/full)

$$
T_{2PAR} = \frac{T-ih^2}{\sqrt{1-0.5h^4i(i-T)}}
$$

and solving this quadratic equation gives:

$$
h^2 = 2T-\sqrt{2} T_{2PAR} \frac{\sqrt{2 -(T^2 - T_{2PAR})(1-\frac{T}{i})}}{2i+(i-T)T_{2PAR}}
$$

All calculations assume that only genetic factors contribute to shared risk in relatives.


