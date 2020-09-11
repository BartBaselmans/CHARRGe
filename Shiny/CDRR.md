Calculating the cross-disorder risk in relatives can be done in a rather similar way as described in the other sections following [Falconer, 1965](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1469-1809.1965.tb00500.x) and [Wray and Gottesman, (2012)](https://www.frontiersin.org/articles/10.3389/fgene.2012.00118/full) using bivariate extensions. 

Here, the threshold $T_{Rx,y}$ bisects the normal distribution for the proportion $K_{Rx,y}$ which is the lifetime risk of disease $y$ in relatives of probands with disease $x$, and which relates to the threshold for disease *x* ($T_{x}$) defined as:

$$
T_{Rx,y} = \frac{T_x - a_{R}r_{g}h_{x}h_{y}i_{y}}{\sqrt{1-a_{R}^2 r^2_{g}h^2_{x}h^2_{y}i_y(i_{y}-T_{y})}},
$$

where $r_{g}$ is the genetic correlation between the traits, and $r_{g}h_{x}h_{y}$ is the co-heritability between the traits (sometimes denoted $h_{x,y}$), and where and $i_x$ and $i_y$ , are the mean phenotypic liabilities of cases for the two diseases. 

Note that the genetic correlation can be estimated using:

$$
r_g = \frac{T_y-T_{Ry,x}\sqrt{1-(1-T_x/i_x)(T_{y}^2-T_{Ry,x}^2)}}{a_R (i_x+(i_x-T_x )T_{Ry,x}^2)\sqrt{h^2_x h^2_y}}
$$

Subsequently, the co-heritability can be estimated using:

$$
h_{x,y} = \sqrt h^2_{x} \times r_g \times \sqrt h^2_{y}
$$

If the heritabilities of two diseases, their lifetime risks and the genetic correlation between them are known, then this equation can be used to estimate $K_{Rx,y}$ as:

$$
K_{Rx,y} = \Phi^{-1}(T_{Rx,y})
$$

where $\Phi^{-1}(x)$ is the inverse of the cumulative standard normal distribution function. And the risk ratio for relatives is defined as: 

$$
RR_{Rx,y} = \frac{K_{Rx,y}}{K_x}
$$

These calculations assume that only genetic factors contribute to the increased risk of disease y in relatives of those with disease x.