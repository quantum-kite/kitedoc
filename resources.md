

<img src=https://user-images.githubusercontent.com/39924384/41094707-9e4ead6e-6a25-11e8-9e16-070a3236c8da.png width="100">

 
##1. <a href="#kite">KITE: bird's-eye view </a>
##2. <a href="#tb">Tight-binding models</a> 
##3. <a href="#se">Spectral expansions</a>
##4. <a href="#app">Applications</a>

* * * 
* * *

## <a name="kite">KITE: bird's-eye view </a> 

KITE evaluates generic electronic response functions and spectral properties of large-scale molecular and solid-state systems by means of extremely accurate spectral expansions of products of single-particle Green's functions [1]. KITE's code uses as input lattice models (tight-binding matrices) of arbitrary complexity that can be imported from standard formats or defined directly via its versatile and user-friendly Pybinding interface.

At the heart of the KITE software is an exact spectral expansion of broadened lattice Green's functions discovered independently by A. Ferreira (KITE's team) in collaboration with E. Mucciolo (U Central Florida) [2] and by A. Braun and P. Schmitteckert (Karlsruhe Institute of Technology) [3]. A large-RAM "single-shot" recursive algorithm developed by A. Ferreira and M. D. Costa (National University of Singapore) enables the evaluation of zero-temperature response functions in systems with multi billions of orbitals $ N \sim 10^{10}$.

In lattice models with small coordination number [$ Z = O(1)$], evaluations of response functions at fixed Fermi energy in large-memory nodes take only a few hours even when billions of spectral coefficients (Chebyshev moments) are retained. This gives access to accuracy and energy resolutions several orders of magnitude beyond previous approaches [4]. To assess generic response functions at finite temperature/frequency, KITE implements the Green's function spectral approach to carry out a direct evaluation of the Kubo-Bastin formula as proposed by L. Covaci and T. G. Rappoport (KITE's team) in collaboration with J. H. García (ICN2) [5].

The pre-release of KITE contains the following functionalities:

- average density of states (DOS) and local DOS;

- generic multi-orbital local (on-site) and bond disorder;
- generic linear response functions for generic orbital observables;
- linear and non-linear optical (AC) conductivity;

To optimize multi-threading and speed up spectral expansions,  KITE provides the option to thread pre-defined partitions in real space (i.e., lattice domains) by means of a domain decomposition algorithm developed by J. Lopes (KITE's team).

For more details about the current pre-release (including a to-do list) refer to the documentation section.

**References**

[1] *Kernel polynomial method*, A. Weiße, G. Wellein, A. Alvermann and H. Fehske. https://link.aps.org/doi/10.1103/RevModPhys.78.275

[2] Critical delocalization of chiral zero energy modes in graphene, A. Ferreira and E. Mucciolo, [Phys. Rev. Lett. 115, 106601 (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.106601).

[3] Numerical evaluation of Green's functions based on the Chebyshev expansion, A. Braun and P. Schmitteckert, [Phys. Rev. B 90, 165112 (2014)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.165112).

[4] Efficient multiscale lattice simulations of strained and disordered graphene, N. Leconte, A. Ferreira, and J. Jung. [Semiconductors and Semimetals 95, 35 (2016)](https://www.sciencedirect.com/science/article/pii/S0080878416300047).

[5] Real-Space Calculation of the Conductivity Tensor for Disordered Topological Matter, J. H. García, L. Covaci, and T. G. Rappoport, [Phys. Rev. Lett. 114, 116602 (2015)](https://doi.org/10.1103/PhysRevLett.114.116602).

* * *

## <a name="tb">Tight-binding models</a>

The tight-binding (TB) method (and the closely related linear combination of atomic orbitals (LCAOs) method in quantum chemistry) is a computationally fast and robust approach to handle large-scale molecular and condensed matter systems. In the TB approximation, electrons are assumed to be strongly bound to the nuclei. The one-particle wavefunctions $ \{\psi_{\alpha}(\mathbf{x})\}$ are approximated by linear combinations of Slater-Koster-type states (i.e., LCAOs) for isolated atoms, i.e.,

$$ \psi_{\alpha}(\mathbf{x})=\frac{1}{\sqrt{N}}\sum_{i=1}^{N}a_{\alpha}^{(i)}\phi_{\textrm{SK}}(\mathbf{x}-\mathbf{x}_{i}),\,\,  (1) $$

where $ i=\{1...N\}$ runs over all sites and orbitals. The one-particle states $ |\psi_{\alpha}\rangle$ are eigenvectors of the parametrized Hamiltonian matrix (the TB Hamiltonian), $\hat H = \sum_{i,j} t_{i,j} |i\rangle\langle j|$. The TB matrix elements—encoding on-site energies ($i=j$) and hopping integrals between different atomic orbitals ($i \neq j$)—can be estimated by means of other methods (e.g., Slater-Koster approach) or by matching the spectrum to that obtained by first-principles calculations in a suitable reference system [1].

Parameterized TB models provide an accurate description of molecular orbitals in molecules and Bloch wavefunctions in many solids. The complexity of TB models grows only linearly with the number of atomic orbitals, providing a basis for large-scale calculations of a plethora of equilibrium and non-equilibrium physical properties, including optical absorption spectra, simulations of amorphous solids, and wave-packet propagation. Disorder, interfaces, and defects can be conveniently added to a TB model by modifying on-site energies and hopping integrals, and adding auxiliary sites. Such a multi-scale approach has proven very successful in describing impurity scattering [2], moiré patterns [3], complex interactions induced by ad-atoms [4], optical conductivity of disordered 2D materials with up to tens of millions of atoms [5], and geometrical properties, vibrational frequencies and interactions of large molecular systems [6].

**References**

[1] *Simplified LCAO Method for the Periodic Potential Problem*, J. C. Slater and G. F. Koster, [Phys. Rev. 94, 1498 (1954)](https://link.aps.org/doi/10.1103/PhysRev.94.1498);

[2] *Elementary prediction of linear combination of atomic orbitals matrix elements*, S. Froyen and W.A. Harrison, [Phys. Rev. B 20, 2420 (1979)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.20.2420);

[3] *Tight-binding modelling of materials*, C. M. Goringe, D. R. Bowler, and E. Hernández, [Rep. Prog. Phys. 60, 1447 (1997)](http://iopscience.iop.org/article/10.1088/0034-4885/60/12/001/pdf); *The Slater–Koster tight-binding method: a computationally efficient and accurate approach*, D. A. Papaconstantopoulos and M. J. Mehl, J[ournal of Physics: Condensed Matter 15, R413 (2003)](http://iopscience.iop.org/volume/0953-8984/15)

[2] *Resonant scattering by realistic impurities in graphene*, T. O. Wehling et al. [Phys. Rev. Lett. 105, 056802 (2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.105.056802); *Unified description of the dc conductivity of monolayer and bilayer graphene at finite densities based on resonant scatterers*, A. Ferreira et al., [Phys. Rev. B 83, 165402 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.165402).

[3] *Ab initio theory of moiré superlattice bands in layered two-dimensional materials*, J. Jung, A. Raoux, Z. Qiao, and A. H. MacDonald, [Phys. Rev. B 89, 205414 (2014)](https://link.aps.org/doi/10.1103/PhysRevB.89.205414).

[4] *Impact of complex adatom-induced interactions on quantum spin Hall phases*. F. J. dos Santos et al., pre-print: a[rXiv:1712.07827 (2017)](https://arxiv.org/abs/1712.07827).

[5] *Numerical calculation of the Casimir-Polder interaction between a graphene sheet with vacancies and an atom*. T. Cysne et al., [Phys. Rev. B 94, 235405 (2016)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.235405).

[6] *A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1−86)*, S. Grimme , C. Bannwarth, and P. Shushkov, [J. Chem. Theory Comput., 13 , 1989 (2017)](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00118)

* * *

## <a name="se">Spectral expansions </a>

The electronic properties of complex molecular and condensed systems are encoded in the eigenvalues and eigenfunctions of Hamiltonian matrices $\hat H$ with large dimension $D$. Direct evaluation of spectral properties and correlation functions generally requires memory of the order of $ D^2$ and the number of floating point operations scale as $ D^3$. Such a large resource consumption clearly restricts the type and size of systems that can be handled by exact diagonalization techniques. Spectral methods offer a powerful and increasingly popular alternative. In the spectral approach, the target functions of interest are decomposed into a spectral series

$ f(E)\propto\sum_{n=0}^{\infty} f_{n}\,P_{n}(E),\,\,  (1)$

where $ P_{n}(E)$ are orthogonal polynomials. The elegance and power of spectral decomposition lie in the fact that the expansion moments $ f_n$ can be obtained by means of  a highly-efficient and stable recursive scheme. Chebyshev polynomials of the first kind

$ T_0(x) = 1,\,\,\, T_1(x)=x,\,\,\, T_2(x) = 2x^2 -1,\,\,\, ...\,,\,\,\, T_{n+1}(x)=2 x T_{n}(x)-T_{n-1}(x),\,\,  (2)$

and $ x\in [-1:1]\equiv \mathcal L$ are a widely used basis functions to approximate generic (non-periodic) functions defined on finite intervals given their unique convergence properties and relation to the Fourier transform [1]. It is easy to verify that the Chebyshev polynomials of first kind satisfy the orthogonality relations

$ \int_{\mathcal{L}}dx\,\omega(x)\,T_{n}(x)\,T_{m}(x)=\frac{1+\delta_{n,0}}{2}\delta_{n,m},\,\,  \textrm{with}\,\,  \omega(x)=1/(\pi\sqrt{1-x^{2}}),\,\,(3)$

and thus form a complete set on $ \mathcal L$. These relations allow to define the numerically convenient spectral decomposition: $ f(x)=\omega(x)\sum_{n=0}^{\infty}\frac{2\mu_{n}}{1+\delta_{n,0}}\,T_{n}(x)$, where $ \mu_n$ are the so-called Chebyshev moments given by the overlap $ \mu_{n}=\int_{\mathcal{L}}dx\,f(x)\,T_{n}(x)$. Efficient numerical implementations evaluate Chebyshev moments "on-the-fly" exploiting the  highly stable recursive procedure [Eq. (2)]. The recursive procedure is stopped when the number of calculated moments allow to retrieve the target function with the desired accuracy.

The extension of these concepts to operators (matrices) allows TB calculations to be carried out for extremely large systems bypassing direct diagonalization. For example, the Chebyshev expansion of the familiar "spectral operator" $ \delta(E-\hat H)$ is given by [2]

$ \delta(E-\hat{H})=\frac{1}{\pi\sqrt{1-E^{2}}}\sum_{n=0}^{\infty}\frac{2}{1+\delta_{n,0}}T_{n}(E)\,\mathcal T_{n}(\hat{H}),\,\,  (4)$

where $ ||\hat H||\le 1$ has been re-scaled to guarantee that its spectrum falls into the  spectral interval $ E\in[-1:1]$. In the above, the operators $ \mathcal{T}_{n}(\hat{H})$ are defined by the matrix version of the standard Chebyshev recursion relations

$ \mathcal{T}_{0}=1, \,\,\,\,\,\, \mathcal{T}_{1}(\hat{H})=\hat{H}, \,\,\,\,\,\,\mathcal{T}_{n+1}(\hat{H})=2\hat{H}\cdot\mathcal{T}_{n}(\hat{H})-\mathcal{T}_{n-1}(\hat{H}).\,\, (5)$

The spectral decomposition (4) allows straightforward determination of several important quantities, for example, the DOS:

$ \rho(E)\equiv\frac{1}{D}\,\textrm{Tr}\,\delta(E-\hat{H})\simeq\frac{1}{\pi\sqrt{1-E^{2}}}\sum_{n=0}^{M-1}\mu_{n}\,T_{n}(E).\,\,\,  (6)$

The Chebyshev moments $ \mu_n = \textrm{Tr}\, T_n(\hat H)/[D(1+\delta_{n,0})/2]$ are evaluated recursively in two steps. First, a series of matrix-vector multiplications is carried out to construct the Chebyshev matrix polynomials using Eq. (5). The complexity of this step for sparse Hamiltonian matrices is $ Z \times D$ (per Chebyshev iteration). Secondly, at the end of a recursive step $ n \rightarrow n+1$, the overlap (trace) is evaluated using an efficient stochastic approach; see below. The knowledge of the Chebyshev moments then allows to reconstruct the DOS over a grid of energies.

Crucially, Chebyshev expansions have uniform resolution due to errors being distributed uniformly on $ \mathcal L$ [1]. In principle, the spectral resolution is only limited by the number of moments retained in the expansion.  As a rule of thumb, the resolution is inversely proportional to the number of moments used, $ \delta E \propto \Delta /M$, where $ M - 1$ is the highest polynomial order and $ \Delta$ is the original spectrum bandwidth (prior to re-scaling). In some situations, a higher number of moments  may be required to converge to a good accuracy (e.g., near a singularity in the DOS).

Truncated spectral representations often present spurious features known as Gibbs oscillations (Fourier expansions) and Runge phenomenon (polynomial expansions).  A well-known example is the "ringing" artifact  in the Fourier expansion of a square wave signal, which persists irrespective of the number of coefficients in the series. Gibbs oscillations can be cured using specialized filtering techniques. An increasingly popular approach in quantum chemistry and computational condensed matter physics is the kernel polynomial method (KPM) [2]. As the name suggests, the KPM makes uses of convolutions with a kernel to attenuate the Gibbs oscillations, e.g., for the DOS, $ \mu_{n}\rightarrow\mu_{n}\times g_{n}$. The kernel $ g_{n}$ must satisfy a number of general conditions to guarantee that the approximate function $ f_M(x)$ converges to the target function $ f(x)$ uniformly as $ M\rightarrow \infty$. A good example is the Lorentz kernel, $ g_n^L = \sinh (\lambda(1-n/M))/ \sinh(\lambda)$, where $ \lambda$ is an adjustable parameter. It has the property that approximates nascent Dirac-delta functions $ \delta_\eta(x)$ by a Lorentzian with resolution $ \eta = \lambda / M$, and thus has been employed to damp Gibbs oscillations in spectral decomposition of  broadened Green's functions [2].

A powerful alternative is given by the Chebyshev polynomial Green's function (CPGF) method [3], which is based on the exact spectral decomposition of the resolvent operator:

$ \hat{\mathcal{G}}(E+i\eta) = \sum_{n=0}^{\infty}g_{n}(E+\imath\eta)\,\mathcal{T}_{n}(\hat{H}),\,\,\,\,\textrm{with        }{g}_{n}(z)\equiv\frac{2i^{-1}}{1+\delta_{n,0}}\frac{\left(z-i\sqrt{1-z^{2}}\right)^{n}}{\sqrt{1-z^{2}}}.\,\,\,(7).$

Differently from KPM, the spectral coefficients depend on the energy. Also, being an exact asymptotic expansion, the convergence of the $ M$th-order approximation to $ \mathcal G(E)$ to a given accuracy is guaranteed provided $ M$ is large enough. Figure 1 shows the convergence of the DOS at the band center $ E=0$ of a giant honeycomb lattice with 3.6 billion sites and dilute random defects. The CPGF method is seen to converge faster than the KPM. More importantly, the CPGF expansion converges asymptotically at all energies (indeed, away from the band center, the Lorentz kernel leads to small errors, which can potentially suppress convergence near a "critical point").
<figure>
  <img src="https://i0.wp.com/quantum-kite.com/wp-content/uploads/2018/02/fig_resources_01.jpg?resize=1024%2C548&ssl=1" alt="convergence"/>
  <figcaption>
Fig 1. Convergence of $ M$-order approximation to the DOS at the band center of a giant honeycomb lattice with $ 60000 \times 60000$ sites and vacancy defect concentration of 0.4% at selected values of energy resolution. As a guide to the eye, we plot the DOS normalized to its converged value (to 0.1% accuracy).  As comparison, the DOS obtained from a KPM expansion with a Lorentz kernel is shown for $ \eta = 1$ meV.
</figcaption>
</figure>

The number of required moments $ M$ drastically increases as electronic states are probed with finer energy resolutions. This poses difficulties when evaluating complex quantities, such as linear response functions given by products of two Green functions (e.g., longitudinal conductivity [4]), especially at finite temperature/frequency (where off-Fermi surface processes are relevant) [5] or when the system approaches a quantum phase transition [3]. To overcome this difficulty, KITE integrates a number of state-of-the-art algorithms:

large-scale "single-shot" algorithm for direct evaluation of zero-temperature DC response functions at  fixed Fermi energy. The algorithm bypasses the expensive recursive calculation of moments and thus can treat giant systems ($ N \propto 10^{10}$) with fine resolution ($ M$ up to hundreds thousands). It can be applied to other quantities, such as  DOS and local DOS; a detailed description is given in Ref. [3];
a Chebyshev-moment approach based on double Chebyshev expansion of the Kubo-Bastin formula developed in [5] gives access to finite temperature response functions of large systems up to ten millions of  orbitals (with $ M$ up to tens thousands [3]); a detailed description of the algorithm is given in Ref. [5] and its large-scale implementation in Ref. [3].
To speed up the evaluation of trace operation $ \textrm{Tr}\{\,\mathcal T_n (\hat H) ... T_m (\hat H)\,\}$, KITE implements the stochastic trace evaluation technique (STE) [2]

$ \rho_{\textrm{STE}}(E)=\sum_{r=1}^{R}\langle r|\,\delta(E-\hat{H})\,|r\rangle,\,\,  (8)$
with random vectors $ |r\rangle=\sum_{i=1}^{N}\chi_{r,i}|i\rangle$. The random variables $ \chi_{r,i}$ are real- or complex-valued [2] (depending on the symmetries of the TB Hamiltonian) and fulfill "white noise" statistics: $ \langle\langle \chi_{r,i} \rangle\rangle=0$, $ \langle\langle \chi_{r,i} \chi_{r^\prime,i^\prime} \rangle\rangle=0$ and $ \langle\langle \chi_{r,i}^* \chi_{r^\prime,i^\prime} \rangle\rangle=\delta_{r,r^\prime}\delta_{i,i^\prime}$.

The STE is extremely accurate for sparse matrices of large dimension (only a few random vectors are needed to converge to many decimal places), which allows substantial savings in computational time. For example, the evaluation of Chebyshev moments of the DOS function requires a total number of operations scaling as

$ P_{\textrm{DOS}} = Z \times N \times M \times R.\, (9)$

The required number of random vectors depends on sparsity of the Chebyshev polynomial matrices $ \mathcal T_n (\hat H)$. For typical TB problems with $ Z\propto O(1)$, in the "thermodynamic limit" ($ N \gg 1$), a single random vector is enough to achieve accuracy of 1% or better [3]. In fact, for sparse matrices, the STE relative error has the favorable scaling $ 1/\sqrt{R N}$. On the other hand, the number of moments required to converge the expansion (say to the same accuracy of the STE) depends strongly on the desired resolution, $ \eta$. As a rule of thumb, $ M$ should not be smaller than a few times the linear dimension of the system $ N^{1/D}$, where $ D$ is the number of spatial dimensions, which then leads to

$ P_{\textrm{DOS}} \propto N^{1+1/D}, \,\,\,\, \textrm{for  }N \gg 1,\,\,\, (10)$

allowing a dramatic reduction in computational time w.r.t. direct diagonalization techniques, especially in $ D\ge 2$.

 **References**

[1] *Chebyshev and Fourier spectral methods*, John P. Boyd, [2nd Ed. Dover][5], New York (2001).

[2] *Kernel polynomial method*, Alexander Weiße, Gerhard Wellein, Andreas Alvermann, and Holger Fehske. [Rev. Mod. Phys. 78, 275 (2006)](https://link.aps.org/doi/10.1103/RevModPhys.78.275)

[3] *Critical delocalization of chiral zero energy modes in graphene*, A. Ferreira and E. Mucciolo, [Phys. Rev. Lett. 115, 106601 (2015)](https://link.aps.org/doi/10.1103/PhysRevLett.115.106601).

[4] *Unified description of the dc conductivity of monolayer and bilayer graphene at finite densities based on resonant scatterers*, A. Ferreira et al., [Phys. Rev. B **83**, 165402 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.83.165402).

[5] *Real-Space Calculation of the Conductivity Tensor for Disordered Topological Matter*, J. H. García, L. Covaci, and T. G. Rappoport,[ Phys. Rev. Lett. **114**, 116602 (2015)][7]. 

* * *

##<a name="app">Applications </a>

#### Optical Conductivity

The real part of the diagonal optical conductivity at zero temperature and finite frequency is given by  [1]

$ \Re\,\sigma(\omega)=\frac{\pi}{\omega\,\Omega}\int_{\mu-\hbar\omega}^{\mu}dE\,\,\textrm{Tr}\,\langle\hat{J}_{x}\,\hat{A}(E)\,\hat{J}_{x}\,\hat{A}(E+\hbar\omega)\rangle_{\textrm{c}},\,\,\, (1)$

where $ \hat{J}_{x}=(ite/\hbar)\sum_{\langle i,j\rangle}(x_{i}-x_{j})(\hat{a}_{i}^{\dagger}\hat{b}_{j}-\textrm{H.c.})$ is the $ x$-component of the current density operator, and

$ \hat{A}(E)=-\frac{1}{\pi}\,\Im\,\frac{1}{E-\hat{H}+i\eta},\,\,\, (2)$

is the spectral operator The symbol $ \langle...\rangle_{\textrm{c}}$ denotes configurational average, $ \Omega$ is the volume of the lattice, $ \mu$ is the chemical potential, and $ \eta$ is a small broadening parameter required for numerical convergence in a finite lattice. Interestingly, the broadening $ \eta=\hbar/\tau_{i}$ mimics the effect of uncorrelated inelastic scattering processes with lifetime  $ \tau$ (e.g., due to phonons) and hence can be viewed as an energy uncertainty due to coupling of electrons to a bath [2].

To obtain an accurate spectral decomposition of the optical conductivity, we use the CPGF operator identity given in Eq. (7) of Sec. 2 to write

$ \hat{A}=-\frac{1}{\pi W}\sum_{n=0}^{\infty}\Im[a_{n}(z)]\,\mathcal{T}_{n}(\hat{h}),\,\,\, (3)$

where $ \hat{h}=\hat H/W$ is the rescaled Hamiltonian (Sec. 2). To simplify notation, here we have assumed a symmetric spectrum $ E = [-W/2,W/2]$. The action of $ \hat{A}$ on a given state vector is computed by standard Chebyshev recursion (Sec. 2). In a numerical implementation, the sum in Eq. (3) is truncated when convergence to a given desired accuracy is achieved. The $ M$th-order approximation to the optical conductivity is therefore given by

$ \Re\:\sigma^{M}(\omega)=\frac{\pi}{\,\omega\,\Omega}\sum_{n,m=0}^{M-1}\sigma_{nm}\,A_{nm}(\mu,\omega),\,\,\, (4) $

where

$ \Re\:\sigma^{(M)}(\omega)=\frac{\pi}{\,\omega\,\Omega}\sum_{n,m=0}^{M-1}\sigma_{nm}\,A_{nm}(\mu,\omega),\,\,\, (5)$
 and

$ A_{nm}(\mu,\omega)=\frac{1}{\pi^{2}W^{2}}\int_{\mu-\hbar\omega}^{\mu}dE\,\alpha_{n}(E)\alpha_{m}(\epsilon+\hbar\omega),\,\,\, (6)$

In the above, $ \alpha_{n}(E)$ is a shorthand for $ \Im[a_{n}\left((E+i\eta)/W\right)]$. The original problem is essentially reduced to the evaluation of the moments $ \sigma_{nm}$, which contains the relevant dynamical information. Once the expansion moments have been determined, the optical conductivity can be quickly retrieved using Eq. (1).

Fig. 2 shows the optical conductivity of graphene with a dilute concentration of vacancies [4]
<figure>
  <img src="https://i2.wp.com/quantum-kite.com/wp-content/uploads/2018/06/optical.png?w=544&ssl=1>
  <figcaption>
Fig. 2 The optical conductivity of graphene with a dilute vacancy concentration of 0.4% at selected values of the chemical potential μ with η ≈ 8 meV.
</figcaption>
</figure>



####Transverse conductivity

The DC conductivity tensor for a given chemical potential and temperature can be written, accoding to Kubo-Bastin formula, as

$ \sigma_{\alpha\beta}(\mu,T)=\frac{\imath \hbar}{\Omega}\int_{-\infty}^{\infty}dE\,\,\textrm{Tr}\,\langle\hat{J}_{\alpha}\,\frac{dG^+}{dE}\hat{A}(E)\hat{J}_{\beta}-\,\hat{J}_{\alpha}\,\hat{A}(E)\hat{J}_{\alpha}\frac{dG^-}{dE}\rangle_{\textrm{c}},\,\,\, (7)$

where $ \hat{J}_{\alpha}$ is the $ \alpha$-component of the current density operator, and $ \hat{A}(E)$ is the spectral operator. The symbol $ \langle...\rangle_{\textrm{c}}$ denotes configurational average, $ \Omega$ is the volume of the lattice, $ \mu$ is the chemical potential $\ T$ is the temperature.

$G^\pm\frac{1}{E-\hat{H}\pm i\eta}$

and $ \eta$ is the small broadening parameter required for numerical convergence in a finite lattice.

The double Chebyshev expansion follows an approach that is equivalent to the one for the optical conductivity. 

The $ M$th-order approximation to the DC conductivity is therefore given by

$ \Re\:\sigma_{\alpha\beta}^{M}(\mu)=\frac{\imath}{ \hbar\,\Omega}\sum_{n,m=0}^{M-1}\sigma_{nm}\,A_{nm}(\mu),\,\,\, (8) $

where $\sigma_{nm}$ has precisely the same structure of the one for optical conductivity while the energy dependent part is now

$ A_{nm}(\mu)=\int_{-\infty}^{\infty}dE\,\Gamma_{m,n}(E)f(\mu,T),\,\,\, (9)$

As all the dependence on temperature is contained in eq. (9), the DC conductivity for different temperatures can be obtained with a single run of the double Chebyshev expansion.

Fig. 3 shows the longitudinal and transverse conductivities for a graphene lattice with $10^5$ unit cells [5].

<figure>
  <img src="https://i1.wp.com/quantum-kite.com/wp-content/uploads/2018/06/hall.png?w=487&ssl=1" alt="hall"/>
  <figcaption>
Fig 3. (a) Longitudinal (solid line) and transverse (dashed line) conductivities for graphene under external magnetic field for different temperatures (b) Density of states and (c) longitudinal conductivity away from the Dirac point where Shubnikov-de Haas oscillations can be observed for kB T /t = 0.002 (black) ,kB T /t = 0.004 (red), kB T /t = 0.008 (blue).
</figcaption>
</figure>


####  Large Systems

<figure>
  <img src="https://i2.wp.com/quantum-kite.com/wp-content/uploads/2018/06/zem.png?w=509&ssl=1" alt="zem"/>
  <figcaption>

DC conductivity for 0.4% vacancy concentration as a function of Fermi energy at selected values of η. The calculation required N = 6.4 × 10^7 Chebyshev moments. The inset shows a zoom of the peak at the Dirac point. Statistical fluctuations of the data are within≃ 1%.
</figcaption>
</figure>



**References**

[1] *Many-Particle Physics*, G. D. Mahan, 3rd Ed. Plenum Press, New York (2001).

[2] *Conductivity of the disordered linear chain*, D. J. Thouless and S. Kirkpatrick, [J. Phys. C **14**, 235 (1981)](http://iopscience.iop.org/article/10.1088/0022-3719/14/3/007/pdf).

[3] *Introduction to Mesoscopic Physics*, Y. Imry, 2nd ed., Oxford University Press, New York (2002).

[4] *Numerical calculation of the Casimir-Polder interaction between a graphene sheet with vacancies and an atom*, T. P. Cysne, T. G. Rappoport, Aires Ferreira, J. V. Lopes, N. M. R. Peres, [Phys. Rev. B **94**, 235405 (2016)]( https://doi.org/10.1103/PhysRevB.83.165402).

[5] *Real-Space Calculation of the Conductivity Tensor for Disordered Topological Matter*, J. H. García, L. Covaci, and T. G. Rappoport, [Phys. Rev. Lett. **114**, 116602 (2015)](https://doi.org/10.1103/PhysRevLett.114.116602). 

[6] *Critical delocalization of chiral zero energy modes in graphene*, A. Ferreira and E. Mucciolo, [Phys. Rev. Lett. 115, 106601 (2015)](https://link.aps.org/doi/10.1103/PhysRevLett.115.106601).
