---
layout: post
title: "T>T: How Electrons Behave Like Sheep Dogs: Bond length in the Hydrogen Molecular Ion"
date: 2021-11-04
excerpt: "We explore the effect of nuclear charge on the bond length in the hydrogen molecular ion."
tags: [science, mathematics, programming, hydrogen, molecular, ion, molecules, chemistry, physics, bond, length, H2+, sheep, dog]
comments: false
math: true
---

The inspiration for this post came from my good friend [Rob Ziolek](https://twitter.com/robziolek) who sent me an [interesting paper](https://pubs.rsc.org/en/content/articlelanding/1972/c3/c39720000568) by J. D. Dunitz and T. K. Ha titled:

*Non-empirical SCF Calculations on Hydrogen-like Molecules : the Effect of Nuclear Charge on Binding Energy and Bond Length.*

By considering the hydrogen molecule, H\\(_2\\) (2 protons, 2 electrons), they showed that the calculated bond length temporarily decreased as the charge of the nuclei increased. This is a counter intuitive result as like charges repel, so we might think larger magnitude like charges should exhibit greater repulsion.

Here we explore whether this behaviour is true for the hydrogen molecular ion (2 protons, 1 electron), H\\(_2^+\\).

# Hydrogen Molecular Ion

The hydrogen molecular ion, H\\(_2^+\\), also referred to as the dihydrogen cation represents the simplest molecular ion in nature, and can be formed from the ionization of the hydrogen molecule, H\\(_2\\). It holds a significant place in the history of physics as the first application of quantum mechanics to a molecule, predicting a bound state solution where classical mechanics could not.

The complete Schrodinger equation cannot be solved **exactly** but by invoking the [Born-Oppenheimer approximation](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Quantum_States_of_Atoms_and_Molecules_(Zielinksi_et_al)/10%3A_Theories_of_Electronic_Molecular_Structure/10.01%3A_The_Born-Oppenheimer_Approximation) the **electronic Schrodinger equation** can be solved exactly using confocal-elliptic coordinates. This allows for complete separation of the Hamiltonian resulting in two coupled differential equations solvable using a recurrence relation.

With the advancement in computer power, we can obtain **very accurate** energies and wave functions to these molecular ions **without assuming the Born Oppenheimer approximation** i.e. including the motion of the nuclei explicitly in the calculation. I have [previously published work which accomplishes this](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.94.042512) and in this post we will use these very high accuracy wavefunctions to explore the effect of nuclear charge on the bond length in hydrogen molecular-like ions with different nuclear charges.

# How to Calculate Bond Length?

The bond length of H\\(_2^+\\) is the distance between the two protons, but how do we calculate this distance? If you read quantum physics textbooks, they tell us the expected value of a physical quantity corresponds to the measured value from an experiment. 

The expectation value corresponds to the average value so can be thought of as an average of all the possible outcomes of a measurement as weighted by their likelihood. Something I found in my paper cited above is that the calculated expectation values of the proton-proton distance did not agree well with the experimental value, whilst the **most probable** value did. Have a look at the following table which lists the calculated values of \\(\langle \Psi \rvert r_3 \rvert \Psi \rangle\\) and the most probable value of \\(r_3\\) for H\\(_2^+\\).


|            | \\(\langle \Psi \rvert r_3 \rvert \Psi \rangle\\) / bohr     | Most probable \\(r_3\\) / bohr |
|            |:------------------------------------------------------------:|:------------------------------:|
|Calculated: | 2.063 913 867                                                | 1.989 685                      |
|Experiment\\(^{[1]}\\): | 1.987 99                                         | 1.987 99                       |

The calculated value of the most probable \\(r_3\\) agrees with experiment to three significant figures, whist the expectation value of \\(r_3\\) agrees to zero significant figures. For the purposes of our investigation we will calculate both of these quantities and see how they compare.

# Results

If you want to know more about how these high accuracy results were calculated see the Method section at the bottom of this post. If you just want to see the results, read on.

The following table shows the calculated expectation value (average value) and most probable value of \\(r_3\\) between the two protons with nuclear charges in the range \\(0.4 \le Z \le 1.3\\). Why 1.3? when the charge of the two 'protons' is 1.3, the system is unbound, shown by the sudden increase in bond length from approximately 2 bohr to > 140 bohr. Physically this means the electron no longer has the influence to hold the two protons in a bound state with itself.

| Nuclear charge / Z | \\(\langle \Psi \rvert r_3 \rvert \Psi \rangle\\) / bohr        | \\(4\pir^2\Psi\Psi\\) / bohr   | \\(\Psi\Psi\\) / bohr  |
|:------------------:|:---------------------------------------------------------------:|:------------------------------:|:----------------------:|
|        0.4         |                    2.296 558 484                                |    2.260 134                   |                        |
|        0.5         |                    2.158 214 212                                |    2.125 930                   |                        |
|        0.6         |                    2.078 038 239                                |    2.048 665                   |                        |
|        0.7         |                    2.035 824 320                                |    2.007 950                   |                        |
|        0.8         |                    2.021 769 795                                |    1.994 954                   |                        |
|        0.9         |                    **2.031 411 074**                            |    **2.004 474**               |                        |
|        1.0         |                    2.063 913 867                                |    2.035 328                   |                        |
|        1.1         |                    2.122 205 373                                |    2.092 229                   |                        |
|        1.2         |                    2.215 691 763                                |    2.185 469                   |                        |
|        1.21        |                    2.228 968 027                                |                                |                        |
|        1.22        |                    2.241 708 265                                |                                |                        |
|        1.23        |                    2.255 127 463                                |                                |                        |
|        1.24        |                    2.269 180 307                                |                                |                        |
|        1.241       |                    2.270 636 718                                |                                |                        |
|        1.242       |                    2.272 101 112                                |                                |                        |
|        1.243       |                    42.561 145 435                               |                                |                        |
|        1.3         |                                                                 |    138.154 519                 |                        |

The bond length can be seen to decrease as the nuclear charge increases between \\(0.4 < Z < 0.8\\) where it then begins to increase up to \\(Z=1.243\\) where the system becomes unbound.

By just considering the two protons the result seems counter-intuitive, but the presence of the electron explains the physics behind this behaviour. The electron has a charge of -1 so when the nuclear charges of the protons are < 1 the electron dominates the electrostatic interactions in the system. Even though the mass of the proton is ~ 1836 times greater than the mass of the electron, it is able to localize the protons to be closer together via its charge influence. As the charge of the protons increases to \\(Z=0.9\\) the charge advantage the electron had is now diminished where the dominance shifts to the protons which have a much higher mass and now comparable charge.

# Electron Density

How does the electron density respond to the increasing nuclear charge? We plot the radial electron density calculated as \\((4\pi r_1 r_2)^2 \Psi \Psi\\) and the electron density, \\(\Psi\Psi\\), around the two protons.


# Method

A brief overview of the method used to numerically solve the Schrodinger equation is now given, for a more indepth discussion consult this [review paper](https://www.researchgate.net/publication/323862960_The_Series_Solution_Method_in_Quantum_Chemistry_for_Three-Particle_Systems) we published. 

* Assume the following wavefunction form, consisting of a triple orthogonal set of Laguerre polynomials

	$$\Psi (u,v,w) =  e^{-\frac{1}{2}(u+v+w)} \times\sum\limits_{l,m,n=0}^{\infty}A(l,m,n) L_{l}(u)L_{m}(v)L_{n}(w),$$

	where \\(u,v,w\\) are perimetric coordinates, linear combinations of inter-particle coordinates.

	$$u = \alpha(r_2+r_{12}-r_1), \ v = \beta(r_1+r_{12}-r_2), \ w= \gamma(r_1+r_2-r_{12}),$$

	where \\(\alpha,\beta,\gamma\\) are non-linear variational parameters.

* Substitute the wavefunction form into the time independent, non-relativistic Schrodinger equation

	$$\left(-\frac{1}{2\mu_{13}}\nabla_1^2 - \frac{1}{2\mu_{23}}\nabla_2^2 -\frac{1}{m_3}\nabla_1\cdot\nabla_2 + \frac{Z_1Z_3}{r_1} + \frac{Z_2Z_3}{r_2} + \frac{Z_1Z_2}{r_{12}}\right)\Psi = E\Psi$$

* Substitute the Laguerre recursion relations

	$$
	\begin{align*}
		xL_n(x) &= -(n+1)L_{n+1}(x)+(2n+1)L_n(x)-nL_{n-1}(x) \\
		xL'_n(x) &= nL_n(x)-nL_{n-1}(x) \\
		xL''_n(x) &= (x-1)L'_n(x)-nL_n(x) 
	\end{align*}
	$$

* Which leads to a 57-term recursion relation

	$$\sum\limits_{a,b,c = -2}^{+2}C_{a,b,c}(l,m,n)A(l+a, m+b, n+c) = 0$$

* Collapsing the triple index, \\(l,m,n\\) down results in

	$$ \sum\limits_{k}^{}C_{ik}B_k = 0$$

* Which we solve as a generalized eigenvalue problem

	$$\sum\limits_{k}^{}\left(H_{ik} - ES_{ik}\right)B_k = 0,$$

	where \\(B_k\\) is the eigenvector (wavefunction).

Modest 2856-term wavefunctions were calculated for hydrogen molecular-like systems with nuclear charges in the range \\(0.4 \le Z \le 1.4\\) in steps of \\(Z = 0.1\\). The calculated energies were accurate to approximately the pico-hartree and the resultant wavefunctions were then used to calculate the average value of \\(r_3\\) and the most probable value of \\(r_3\\). The most probable value is the maximum position of the [**intracule density**](https://en.wikipedia.org/wiki/Intracule), the particle density between the two protons, calculated here using parallelized numerical integration in C++.

# References

[1] Taken from [NIST](http://webbook.nist.gov/chemistry/)