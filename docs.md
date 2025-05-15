# Documentation: Standard Map Analysis Solution

This document explains the MATLAB implementation and mathematical concepts
behind the analysis of the Chirikov Standard Map, a chaotic dynamical system.

## Problem Overview

The Chirikov Standard Map (or simply Standard Map) is an area-preserving map on
the square [0, 2π] × [0, 2π] defined by:

$$F(I, \theta) = \begin{pmatrix} I + K \cdot \sin(\theta) \mod 2\pi \\ \theta + I + K \cdot \sin(\theta) \mod 2\pi \end{pmatrix}$$

With recursion equations:

$$I_{n+1} = I_n + K \cdot \sin(\theta_n) \mod 2\pi$$

$$\theta_{n+1} = \theta_n + I_{n+1} \mod 2\pi$$

Where $K \geq 0$ is a bifurcation parameter controlling the degree of chaos in
the system.

## Solution Approach

The solution uses a pure functional programming paradigm in MATLAB, with the
following components:

1. Phase portrait generation for different K values
2. Calculation of Lyapunov exponents
3. Analytical verification of Lyapunov exponent properties

## Part (a): Phase Portraits

### Implementation Details

1. **Core Functions**:
   - `standard_map_step`: Computes one iteration of the standard map
   - `generate_trajectory`: Creates a trajectory from an initial state
   - `generate_phase_portrait`: Generates multiple trajectories for
     visualization
   - `plot_phase_portrait`: Visualizes the phase portrait

2. **Parameters**:
   - 50 different trajectories per phase portrait
   - Each trajectory has a length of 1000 iterations
   - Three different K values chosen from specified ranges

3. **Visualization**:
   - Each phase portrait shows trajectories in the $\theta$-$I$ plane
   - Axes are labeled with appropriate tick marks (0, π/2, π, 3π/2, 2π)

### Physical Interpretation

As K increases:

- For small K (≤ 0.6): Mostly regular behavior with visible KAM
  (Kolmogorov-Arnold-Moser) tori
- For medium K (0.9-1.1): Partial breakdown of regular structures, some chaotic
  regions appear
- For large K (1.4-2.0): Significant chaotic behavior with few remaining islands
  of stability

This illustrates the transition from regular to chaotic dynamics as the
nonlinearity parameter K increases.

## Part (b): Lyapunov Exponents

### Implementation Details

1. **Core Functions**:
   - `jacobian`: Calculates the Jacobian matrix (system matrix) for the standard
     map
   - `compute_lyapunov`: Calculates Lyapunov exponents for a single K value
   - `compute_lyapunov_spectrum`: Computes Lyapunov exponents across a range of
     K values

2. **Calculation Method**:
   - QR decomposition algorithm to avoid numerical overflow
   - Initial transient iterations discarded for better convergence
   - Both Lyapunov exponents calculated and plotted against K

3. **Parameters**:
   - K ranging from 0 to 4 with 100 discrete values
   - 1000 iterations for exponent calculation after 200 transient iterations

### Physical Interpretation

Lyapunov exponents quantify the rate of separation of infinitesimally close
trajectories:

- Positive exponent: Exponential divergence of nearby trajectories (chaos)
- Zero exponent: Linear divergence or no divergence
- Negative exponent: Convergence of nearby trajectories

In the plot:

- For small K: Both exponents are near zero (regular dynamics)
- As K increases: One exponent becomes increasingly positive while the other
  becomes increasingly negative
- The sum of exponents remains zero (area-preserving property)

## Part (c): Analytical Verification

### Questions and Answers

1. **Value of the determinant of the system matrix DF(I,θ)**:

   The system matrix is:
   $$DF(I, \theta) = \begin{pmatrix} 1 & K\cos(\theta) \\ 1 & 1+K\cos(\theta) \end{pmatrix}$$

   The determinant is:
   $$\det(DF(I, \theta)) = 1 \cdot (1+K\cos(\theta)) - K\cos(\theta) \cdot 1 = 1$$

   The determinant equals 1, independent of I and θ.

2. **Property of orthogonal matrices Qn**:

   Orthogonal matrices have the property that $Q^T Q = I$ (the identity matrix).
   This implies that $\det(Q_n) = \pm 1$. Since QR decomposition typically
   produces matrices with positive determinants, $\det(Q_n) = 1$.

3. **Consequence for det(An)**:

   From the recurrence relation: $$A_{n+1} = DF(I_n, \theta_n) \cdot Q_n$$

   Using the multiplicative property of determinants:
   $$\det(A_{n+1}) = \det(DF(I_n, \theta_n)) \cdot \det(Q_n) = 1 \cdot 1 = 1$$

4. **Consequence for det(Rn)**:

   From the QR decomposition: $$A_{n+1} = Q_{n+1} \cdot R_{n+1}$$

   Therefore: $$\det(A_{n+1}) = \det(Q_{n+1}) \cdot \det(R_{n+1})$$

   Since $\det(A_{n+1}) = 1$ and $\det(Q_{n+1}) = 1$, we have:
   $$\det(R_{n+1}) = 1$$

5. **Relationship to diagonal elements of Rn**:

   For an upper triangular matrix R, the determinant is the product of the
   diagonal elements: $$\det(R_n) = (R_n)_{11} \cdot (R_n)_{22} = 1$$

   Therefore, $(R_n)_{11} \cdot (R_n)_{22} = 1$

6. **Value of ln((Rn)11) + ln((Rn)22)**:

   Taking logarithms of both sides:
   $$\ln((R_n)_{11} \cdot (R_n)_{22}) = \ln(1) = 0$$

   By the properties of logarithms: $$\ln((R_n)_{11}) + \ln((R_n)_{22}) = 0$$

   This holds for all n.

7. **Consequence for Lyapunov exponents λ1 and λ2**:

   The Lyapunov exponents are defined as:
   $$\lambda_i = \lim_{N \to \infty} \frac{1}{N} \sum_{n=0}^{N-1} \ln((R_n)_{ii})$$

   Since $\ln((R_n)_{11}) + \ln((R_n)_{22}) = 0$ for all n, summing over n:
   $$\sum_{n=0}^{N-1} \ln((R_n)_{11}) + \sum_{n=0}^{N-1} \ln((R_n)_{22}) = 0$$

   Dividing by N:
   $$\frac{1}{N}\sum_{n=0}^{N-1} \ln((R_n)_{11}) + \frac{1}{N}\sum_{n=0}^{N-1} \ln((R_n)_{22}) = 0$$

   As N approaches infinity: $$\lambda_1 + \lambda_2 = 0$$

   This shows that the sum of the two Lyapunov exponents must be zero.

8. **Verification in the diagram**:

   In our plot of Lyapunov exponents vs. K, we can observe that the sum of λ1
   and λ2 is consistently zero across all K values. This confirms the
   theoretical result derived above and validates our numerical implementation.

   This property is a direct consequence of the area-preserving nature of the
   Standard Map, which maintains the phase space volume under time evolution.

## Conclusion

The Chirikov Standard Map demonstrates the transition from regular to chaotic
dynamics as the nonlinearity parameter K increases. The analysis of phase
portraits and Lyapunov exponents provides quantitative measures of this
transition. The area-preserving property of the map is reflected in the
relationship between the Lyapunov exponents, where their sum remains zero
regardless of the level of chaos in the system.

This property is a manifestation of Liouville's theorem in Hamiltonian dynamics,
which states that the phase space volume is conserved under time evolution.
