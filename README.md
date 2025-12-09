# Bio-project-slime-mold

# Bio-project
This project implements a simplified network model for memory encoding in tube diameter hierarchy of P. Polycephalum. A chemical softening agent flows through network edges and modifies their radius depending on local concentration. 

## Overview

Each edge represents a cylindrical tube with radius (a_i). The concentration of softening agent locally reduces Young’s modulus -> radius adaptation while keeping total volume constant.

Simulation steps:

1. inject softening agent at a root edge
2. redistribute chemical to neighboring edges (advective approximation)
3. update Young’s modulus as a function of concentration
4. update radii under a global volume constraint

---

## The model

### Concentration

$c_i=\frac{n_i}{a_i^2}$

### Young modulus

$E_i = E_0 - \frac{dE  c_i}{c_i + c_0}$

### Mechanical equilibrium

$(a_i-a_0)E_i=\alpha$

### Volume conservation

$\sum_i a_i^2 = V$

### Minimization and update

Minimize
$\sum_i\left((a_i-a_0)-\frac{\alpha}{E_i}\right)^2$
with volume constraint, which gives
$a_i=\frac{a_0+\alpha/E_i}{1+h}$
where (h) is the Lagrange multiplier ensuring $(\sum_i a_i^2 = V)$.

### Advection of softening agent

Redistribution is proportional to cross-sectional area:
$\Delta n_{i\to j}\propto a_j^2$

## Parameters
$\alpha$: mechanical sensitivity , E: base Young modulus ,  dE: softening strength , c0: reference concentration , Nk: injected softening agent per step 

---
## Interpretation


 Edges with higher concentration soften more and increase their radius; larger edges transport more chemical agent, reinforcing their size. 
