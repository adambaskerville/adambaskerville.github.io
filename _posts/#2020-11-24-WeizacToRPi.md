---
layout: post
title: "T>T: From WEIZAC to RPi: Quantum Mechanics on a Raspberry Pi Zero"
date: 2020-11-24
excerpt: "We investigate how Chaim. L. Pekeris accurately calculated the quantum mechanical properties of the helium atom using the WEIZAC computer and implement his method on a raspberry pi zero."
tags: [science, mathematics, programming, matplotlib, numpy, python, pekeris, quantum, helium]
comments: false
math: true
---

In 1949, Chaim, L, Pekeris a mathematician and geophysicist received a reply to a letter he sent to John Von Neumann, the foremost mathematician of the time, regarding a discussion on the properties of the helium atom:

"*The differential equations which you give as determining the problem of Helium II do not look vicious. How bad are they from the numerical point of view? They should certainly not be bad for high-speed computer. Quite apart from more ambitious devices, the ENIAC or the SSEC, 33 or the Harvard Mark II or III should be sufficient here ... It would certainly be good to discuss these and other similar problems personally when you are here.*"

Pekeris would later go on to develop a method to numerically solve the time independent, non-relativistic Schrodinger equation for the helium atom, where he used the Weizmann Automatic Computer, **WEIZAC**, the first computer built in Israel.

In this post we are going to write a simple Python program to replicate Pekeris' method and run it on the £9, 6cm wide **Raspberry Pi zero W** computer, a far cry from the room sized behemoth WEIZAC which in the 1950s cost $50,000 - equivalent to approximately half a million dollars in today's money.

# WEIZAC Computer

The WEIZAC (Weizmann Automatic Calculator) was the first computer to be constructed in Israel and one of the first large-scale, stored-program, electronic computers in the world. It was built at the Weizmann Institute during 1954-1955, based on the Institute for Advanced Study (IAS) architecture developed by John von Neumann. The project was initiated by Pekeris who originally wanted the computer to solve Laplace's tidal equations for Earth's oceans as well as to benefit the scientific community. Being 1953 it was difficult to source components for such a machine; some were imported, but others were clever adaptations, such as thin copper strips that came from a small local bicycle-part shop.

It performed its first calculation in October 1955 and remained in demand by Israeli scientists up until 1963 when it was replaced by a commercially-built computer, a [CDC 1604A](https://en.wikipedia.org/wiki/CDC_1604).

# Raspberry Pi Computer

Raspberry Pi's are small single-board computers developed in the United Kingdom by the [Raspberry Pi Foundation](https://www.raspberrypi.org/). These are typically used by computer and electronic hobbyists, due to their adoption of HDMI and USB devices. In this post we use the [Raspberry Pi Zero W](https://www.raspberrypi.com/products/raspberry-pi-zero-w/) which is one of their smallest and cheapest boards. They are designed to be compact and powerful enough for simple applications and have found use in a [variety of different applications](https://www.makeuseof.com/tag/different-uses-raspberry-pi/). I use one as a all-in-one astrophotography camera and viewer for my telescope where I wrote a GUI application called [AstroPitography](https://adambaskerville.github.io/tabs/astro/). 

# WEIZAC vs. Raspberry Pi Zero W

The following table compares two key specifications used in modern computation, the memory and instructions per second between both the WEIZAC and Raspberry Pi Zero W (RPi).

| Specification           | WEIZAC      |  RPi       |
|:-----------------------:|:-----------:|:----------:|
| Memory / Bytes          | ~ 8192      | 512 000 000|
| Instructions per second |   | 1 000 000 000

It is easy to forget how rapidly computation has advanced in the last 60 years, but the cheap, credit card sized Raspberry Pi Zero W **greatly** outperforms WEIZAC, one of the most powerful computers in the world back in the 1950s.

We now describe how Pekeris developed a technique to numerically solve the tiem independent, non-relativistic Schrodinger equation for two-electron atoms. There are two sections:

1. **Series Solution Method for Helium Atom: Brief** This gives a brief overview of the methodology if you are not feeling in a mathematically inclined mood.
2. **Series Solution Method for Helium Atom: Detailed** This gives a more detailed description of the method for those who want to crunch the numbers.
# Series Solution Method for Helium Atom: Brief

The series solution method developed by Pekeris is not trivial so this section provides a bullet point summary:

* The Schrodinger equation for > 2 particles is not exactly solvable.
* Helium contains three particles meaning standard mathematical techniques will not yield an exact solution.
* Pekeris re-derived the problem into perimetric coordinates, linear combinations of the sides of the triangle to circumvent certain mathematical difficulties owing to the triangle inequality.
* Pekeris assumed the wavefunction form from knowledge of the asymptotic behaviour of the hydrogen atom.
* He used the Laguerre recurrence relations to decompose all partial derivatives in the Schrodinger equation.
* This resulted in a 33-term recurrence relation (effectively a sum).
* This was used to find energy eigenvalues by calculating the vanishing of the matrix determinant.
# Series Solution Method for Helium Atom: Detailed

Returning to the focus of this post, Pekeris was interested in numerically solving the Schrodinger equation for the helium atom, a two-electron, three-particle system. The wavefunction of the ground state, \\(\Psi\\) only depends on three coordinates, \\(r_1, r_2, r_{12}\\). Using these coordinates, the Schrodinger equation becomes

\\[
\underbrace{\frac{\partial^2\Psi}{\partial r_1^2} + \frac{2}{r_1}\frac{\partial \Psi}{\partial r_1} + \frac{\partial^2\Psi}{\partial r_2^2} + \frac{2}{r_2}\frac{\partial\Psi}{\partial r_2} + 2\frac{\partial^2\Psi}{\partial r_{12}^2} + \frac{4}{r_{12}}\frac{\partial\Psi}{\partial r_{12}} + \frac{(r_1^2 - r_2^2 + r_{12}^2)}{r_1r_{12}}\frac{\partial^2\Psi}{\partial r_1\partial r_{12}} + \frac{(r_2^2 - r_1^2 + r_{12}^2)}{r_2r_{12}} \times \frac{\partial^2\Psi}{\partial r_2 \partial r_{12}}}_{\text{Kinetic energy}} + \underbrace{2\left(E + \frac{Z}{r_1} + \frac{Z}{r_2} - \frac{1}{r_{12}}\right)\Psi}_{\text{E \& Potential energy}} = 0
\\]

Using the inter-particle coordinates, \\(r_1, r_2, r_{12}\\) results in an undesired coupling between coordinate ranges, i.e. \\(r_1\\) and \\(r_2\\) may range from \\(0 \rightarrow \infty\\) but then \\(r_{12}\\) will range from \\(|r_1 - r_2| \rightarrow r_1+r_2\\) in order to correctly describe the shape of the triangle (triangle inequality). Numerically this is highly undesirable as we would need to write a complicated routine to handle this coupling. Pekeris realised that using **perimetric coordinates** diverts all these computational difficulties allowing for each of the perimetric coordinates to range from \\(0 \rightarrow \infty\\)

\\[
\begin{aligned}
    u &= \epsilon(r_2 + r_{12} - r_1) = \epsilon z_1 \\
    v &= \epsilon(r_1 + r_{12} - r_2) = \epsilon z_2 \\
    w &= 2\epsilon(r_1 + r_2 - r_{12}) = 2\epsilon z_3
\end{aligned}
\\]
Figure shows what the inter-particle and perimetric coordinates look like

The exactly solvable hydrogen atom shows that the expected asymptotic behaviour of the wavefunction, \\(\Psi\\), will decay exponentially

\\[
\Psi = e^{\frac{1}{2}(u+v+w)}F(u,v,w),
\\]

where \\(F(u,v,w)\\) is the part of the wavefunction we have not figured out yet. Pekeris expressed this unknown function as a series expansion of products of Laguerre polynomials

\\[
F = \sum\limits_{l,m,n=0}^{\infty}A(l,m,n)L_l(u)L_m(v)L_n(w),
\\]

where \\(L_n(x)\\) denotes the Laguerre polynomial

\\[
    L_n(x) = \sum\limits_{k=0}^{n}\binom{n}{k}\frac{(-x)^k}{k!}
\\]

All families of classical orthogonal polynomials were derived as solutions to linear differential equations, i.e. the Laguerre polynomials were discovered as solutions to a particular differential equation. They obey a set of mixed differential-recurrence relations

\\[
\begin{aligned}
    xL_n''(x) &= (x-1)L_n'(x) - nL_n(x) \\
    xL_n'(x)  &= nL_n(x) - nL_{n-1}(x) \\
    xL_n(x)   &= -(n+1)L_{n+1}(x) + (2n+1)L_n(x) - nL_{n-1}(x) 
\end{aligned}
\\]

where primes denote derivatives with respect to \\(x\\). Pekeris then substituted the expansion for \\(F(u,v,w)\\) in terms of yet-unknown coefficients, \\(A(l,m,n)\\) into the transformed Schrodinger equation expansion, resulting in a 33-term linear partial recurrence equation eith polynomial coefficients satisfied by \\(A(l,m,n)\\).

This recurrence yields a homogenous linear system of equations with \\(\infty^3\\) equations and \\(\infty^3\\) unknowns, but for eigenvalues, \\(\epsilon\\), the matrix determinant vanishes and there are solutions. The largest of these eigenvalues corresponds to the ground state energy of the atom.

No computer can handle infinite systems so Pekeris developed a means of truncating the problem, considering only those \\(l,m,n \ge 0\\) for which \\(l,m,n \le \omega\\) where \\(\omega\\) is a integer representing the **matrix order**.

Exploiting fermionic exchange in two-electron atoms allows for the system to be reduced once more, requiring either symmetry, \\(A(l,m,n) = A(m,l,n)\\) **para states** and \\(A(l,m,n) = -A(m,l,n)\\) **ortho states**.

In order to run all this on a computer, one needed a convenient way to order linearly all triplets of integers, \\(l,m,n\\). Pekeris devised a bijective map \\(k:\lbrace (l,m,n) \in \mathbb{N}_0^3 |l \le m \rbrace \rightarrow \mathbb{N}\\), the details of which are not important but just represent a systematic means of ordering three indices.

We now construct our Python program to implement the Pekeris series solution method.
# Implementation

Pekeris originally derived the 33-term recursion relation by hand, but here I pre-calculated it using Maple and saved it to a text file called `RR.txt` which our Python program will call.

## Loading 33-term Recursion Relation

There are a number of ways to read external data files into python, but here we use `NumPy`, specifically the [`genfromtxt`](https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html) command.

```python
from numpy import genfromtxt

RR = genfromtxt("RR.txt")
```


## Building Matrices

## Solving Generalized Eigenvalue Problem

## Complete Program


# Conclusions

# References

[1] Von Neumann to Pekeris, Feb. 23, 1950 (WIA 3-96-72).

[2] L. Carry and R. Leviathanl, *WEIZAC: An Israeli Pioneering Adventure in Electronic Computing (1945–1963)*, Springer, 2019, [10.1007/978-3-030-25734-7](10.1007/978-3-030-25734-7)

[3] C. Koutschan, and D. Zeilberger, *Math. Intelligencer*, **33**, 52--57, 2011.