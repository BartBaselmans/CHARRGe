---
title: "test"
author: "Bart Baselmans"
date: "21/08/2020"
output: html_document
---
## Liability Threshold Model

When many factors contribute to the risk of disease it can be helpful to consider a model of disease where there is a latent distribution of liability to disease underlying the dichotomous case/control status.

Since this latent liability comprises many genetic and other risks it is reasonable to assume that the liability distribution is approximately normally distributed (since many things added together will make a bell-shaped distribution in a population sample) and that those affected by disease have the combination of genetic and other factors that it places them in the top end of the liability distribution. 

Hence, those whose liability to disease is above the threshold are affected, and so this model is sometimes called the liability threshold model. The threshold can be calculated using the lifetime risk ($K$) parameter as $K$ represents the proportion of individuals in the general population have the disorder of interest. The threshold (or boundary) that determines the area of $K$ (indicated in blue) in a normal distribution is given by:
$$
T=\Phi^{-1}(1-K)
$$

Conversely, the life-time risk ($K$) is derived from the threshold parameter
$$
K = 1-\Phi(T)
$$

where $\Phi^{-1}(x)$ is the inverse of the cumulative standard normal distribution function. 

The lifetime risk of a disease in relatives ($K_{R}$) of those affected by the disease is expected to be greater, or equal to, the lifetime risk of the disease in a population sample. It is logical to assume that the threshold in liability associated with disease has the same value in the relatives of those affected. As a result, the liability distribution in first degree members must be shifted (in the direction of increased liability) compared to the general population to an extent consistent with the observed higher risk of $K_{R}$.\

We can then ask: what proportion of the variation in liability must be attributable to genetic factors (i.e., what is the heritability?) in order to generate this observed increased risk in relatives, given the known coefficient of relationship between the relatives, $a_R$, i.e., 0.5 for 1st degree relatives.

With both $T$ and $T_R$ known as well $a_R$, it is possible to calculate the heritability using:

$$
h^2 = \frac{T - T_R \sqrt{1 -(1-\frac{T}{i})(T^2-T^2_R)}}{a_R(i+(i-T)T^2_R)}
$$

Although, this equation looks complicated, in fact everything on the right hand side can be calculated from two observations $K$ and $K_R$, and $a_R$ (and holds as long as $a_R$ < 1), and  if ascertainment is based on the disease status of one proband. 

Note, that $i$ is the mean phenotypic liability. which can be calculated using $i$ = $z/K$, where $z$ is the height of the standard normal curve at threshold $T$.

