# Shock Tubes using Finite Volume Methods

One-dimensional finite volume solver to simulate a shock tube using the Roe scheme.
A first-order (in both space and time) one-dimensional finite volume scheme to solve shock-tube
problems using Roe's approximate Riemann solver with the Entropy fix.
(b) Run the standard Sod's shock tube case with [L; uL; pL] = [1:0; 0; 1:0] and [R; uR; pR] = [0:125; 0; 0:1]
using N = 80 volumes. Plot the nal density, velocity, and pressure distributions.
(c) Repeat using N = 160, 320, and 640 and comment on the resolution of the expansion fan, contact
discontinuity, and shock at each level of mesh renement.

Figure 5 shows the density, velocity, and pressure distributions at N=2000.‎ The resolution of the expansion fan, contact discontinuity, and shock increases as the level of mesh refinement increases. Contact discontinuities (surfaces that separate zones of different density and temperature) becomes sharper as N increases. The expansion fan becomes sharper too. The shock become sharp and its distance decreases as N increases.


The results of solution at nx=ny=33 are shown in figure 1 (a) - (d).

**Figure 1: Velocity & pressure contours and line plot of velocity components at nx=ny =33.**

**Figure 1 (a):** ![image](https://user-images.githubusercontent.com/89004966/152657677-c2e91f6f-b573-4920-a4d9-470f995c4126.png)

**Figure 1 (b):** ![image](https://user-images.githubusercontent.com/89004966/152657745-0452d52c-776e-4ae5-a5bc-4fbf2504a084.png)

**Figure 1 (c):** ![image](https://user-images.githubusercontent.com/89004966/152657788-7275f569-4821-43c4-84ce-4d528baba985.png)

**Figure 1 (d):** ![image](https://user-images.githubusercontent.com/89004966/152657807-4dff1bc7-1d2e-48a2-8bf6-d1bb5ff41ced.png)





