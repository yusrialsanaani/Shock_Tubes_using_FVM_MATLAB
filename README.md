# Shock Tubes using Finite Volume Methods

One-dimensional finite volume solver to simulate a shock tube using the Roe scheme.
A first-order (in both space and time) one-dimensional finite volume scheme to solve shock-tube
problems using Roe's approximate Riemann solver with the Entropy fix.
(b) Run the standard Sod's shock tube case with [L; uL; pL] = [1:0; 0; 1:0] and [R; uR; pR] = [0:125; 0; 0:1]
using N = 80 volumes. Plot the nal density, velocity, and pressure distributions.
(c) Repeat using N = 160, 320, and 640 and comment on the resolution of the expansion fan, contact
discontinuity, and shock at each level of mesh renement.

The conservative form of 1D Euler equations [1]:
∂W/∂t+∂F/∂x=0	(1)
W=(█(ρ@ρu@ρE)),F=(█(uρ@p+ρu^2@uρE)),p=ρRT or  p=ρ(γ-1)(E-1/2 u^2 ),H=E+p/ρ
ρE=p/((γ-1))+ρ/2 u^2  or E=p/((γ-1)ρ)+1/2 u^2
where ρ is the density, p is the pressure, u is the velocity (horizontal), E is the internal energy, H is the static enthalpy, and γ is the ratio of specific heats. Equation (1) can be rewritten as follows:
∂W/∂t+A ∂W/∂x=0,A=∂F/∂x
(2)
A is the convective flux Jacobian matrix which is given as:
A=(■(0&1&0@1/2 u^2 (-1+γ)-u^2&u(3-γ)&γ-1@u(1/2 u^2 (-1+γ)-H)&H-u^2 (-1+γ)&uγ))	(3)
A can be expressed as A=RΛL where Λ is the diagonal matrix of real eigenvalues of A, L and R are the left and right matrices of eigenvectors of A; L=R^(-1).
L=(■(1-(u^2 (γ-1))/(2α^2 )&u(γ-1)/α^2 &(1-γ)/α^2 @1/2 u^2 (γ-1)-uα&α-u(γ-1)&γ-1@uα+1/2 u^2 (γ-1)&-α-u(γ-1)&γ-1)),
R=(■(1&β&β@u&(u+α)β&(u-α)β@u^2/2&(H+uα)β&(H-uα)β)),α^2=(γ p)/ρ,β=1/(2α^2 )	(4)
where a is the speed of sound. 
With 〖W(x,0)=W〗_0 (x), equation (1) can be written as:
(∂W_i)/∂t+(F_(i+1⁄2)-F_(i-1⁄2))/Δx=0	(5)
W_i is the cell-centered values stored at nods and F_(i±1⁄2)  are the fluxes at left and right cell interface.
F_(i+1⁄2)=1/2 (〖F^R〗_(i+1⁄2)+〖F^L〗_(i+1⁄2) )-1/2 RΛL(〖W^R〗_(i+1⁄2)+〖W^L〗_(i+1⁄2) )	(6)
Λ is the diagonal matrix of absolute eigenvalues of A and L and R. The Roe solver is an approximate Riemann solver and it uses below formulation to find the numerical fluxes at the interface:
F_(i+1⁄2)=1/2(〖F^R〗_(i+1⁄2)+〖F^L〗_(i+1⁄2))-1/2  R┴-  Λ┴-  L┴- (〖W^R〗_(i+1⁄2)+〖W^L〗_(i+1⁄2))	(7)
The bars denote the Roe’s averages between the left and right interfaces. To get A┴-, u┴-, H┴-, and a┴-  should be calculated first using Roe average based on expansion of Taylor series for  F about W^L,W^R  points.
F(W)=F(W^L )+A(W^L )((W-W^L )),F(W)=F(W^R )+A(W^R )((W-W^R ))
F(W^L)-F(W^R)=A(W^L-W^R)+(A(W^R)-A(W^L))W
Roe’s A┴- , A┴- (W^L-W^R )   goes to A(W) as W┴-  goes to W.
Therefore: F(W^L )-F(W^R )=A┴- (W^L-W^R )	(8)
The Roe averaging formulas to compute approximate values for constructing R┴- and L┴- are given below
u┴-=(u_R √(ρ_R )+u_L √(ρ_L ))/(√(ρ_R )+√(ρ_L )),   H┴-=(H_R √(ρ_R )+H_L √(ρ_L ))/(√(ρ_R )+√(ρ_L )),   α┴-=√((γ-1)[H┴--1/2 u┴-^2 ] )
The eigenvalues of the Jacobian matrix are
λ_1=u┴-,λ_2=u┴-+α,λ_3=u┴--α
(9)
After obtaining A_(i+1⁄2), a first-order upwind Roe’s scheme is used. The conservative form of flux for Roe’s solver can be written as:
F(〖W^t〗_L,〖W^t〗_R)=1/2(F(〖W^t〗_L)+F(〖W^t〗_R)-1/2  R┴-  Λ┴-  L┴- (〖W^t〗_L+〖W^t〗_R)	(10)
The first-order upwind Roe’s scheme can be written as:
〖W^(t+1)〗_i=〖W^t〗_i-Δt/Δx (F(〖W^t〗_L,〖W^t〗_R )-F(〖W^t〗_L,〖W^t〗_R ))	(11)
Figure 5 shows the density, velocity, and pressure distributions at N=2000.‎ The resolution of the expansion fan, contact discontinuity, and shock increases as the level of mesh refinement increases. Contact discontinuities (surfaces that separate zones of different density and temperature) becomes sharper as N increases. The expansion fan becomes sharper too. The shock become sharp and its distance decreases as N increases.


The results of solution at nx=ny=33 are shown in figure 1 (a) - (d).

**Figure 1: Velocity & pressure contours and line plot of velocity components at nx=ny =33.**

**Figure 1 (a):** ![image](https://user-images.githubusercontent.com/89004966/152657677-c2e91f6f-b573-4920-a4d9-470f995c4126.png)

**Figure 1 (b):** ![image](https://user-images.githubusercontent.com/89004966/152657745-0452d52c-776e-4ae5-a5bc-4fbf2504a084.png)

**Figure 1 (c):** ![image](https://user-images.githubusercontent.com/89004966/152657788-7275f569-4821-43c4-84ce-4d528baba985.png)

**Figure 1 (d):** ![image](https://user-images.githubusercontent.com/89004966/152657807-4dff1bc7-1d2e-48a2-8bf6-d1bb5ff41ced.png)





