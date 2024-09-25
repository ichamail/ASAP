# Subsonic Potential Aerodynamics

# A 3D Aerodynamic Potential-Flow Code

## Theoretical Model

![Potential Flow Theoretical Model](https://github.com/ichamail/Panel-Methods-3D/assets/107580530/a0a380f7-a822-4470-b100-ed94ea4238ed)

* A vector field $` \underline{V}: \mathbb{U} \to \mathbb{R}^n `$, where $` \mathbb{U} `$ is an open subset of $` \mathbb{R}^n `$, is said to be conservative if  there exists a $` \mathrm{C}^1 `$ (continuously differentiable) scalar field $` \phi `$ on $` \mathbb{U} `$ such that $` \underline{V} = \nabla \phi `$.

* According to Poincaré's Lemma, A continuously differentiable ($` \mathrm{C}^1 `$) vector field $` \underline{V} `$ defined on a simply connected subset $` \mathbb{U} `$ of $` \mathbb{R}^n `$  ($` \underline{V} \colon \mathbb{U} \subseteq \mathbb{R}^n \to \mathbb{R}^n `$), is conservative if and only if it is irrotational throughout its domain ($` \nabla \times \underline{V} = 0 `$, $` \forall \underline{x} \in \mathbb{U} `$).

* Circulation $` \Gamma = \oint_{C} \underline{V} \cdot \mathrm{d} \underline{l} = \iint_S \nabla \times \underline{V} \cdot \mathrm{d}\underline{S} `$.
* In a conservative vector field this integral evaluates to zero for every closed curve. $` \Gamma = \oint_{C} \underline{V} \cdot \mathrm{d} \underline{l} = \iint_S \nabla \times \underline{V} \cdot \mathrm{d}\underline{S} = \iint_S \nabla \times \nabla \phi \cdot \mathrm{d}\underline{S} = 0 `$ 



### Velocity field $` \underline{V} \colon \mathbb{U} \subseteq \mathbb{R}^3 \to \mathbb{R}^3 `$

$`
\begin{array}{l}
\bullet \text{ incompressible: } \nabla \cdot \underline{V} = 0  \\
\begin{drcases} 
\bullet \text{ irrotational: } \nabla \times \underline{V} = 0
\\
\bullet \text{ simply connected domain } \mathbb{U} 
\end{drcases} \implies \underline{V} = \nabla \phi \text{ (conservative)}
\end{array}
`$

### Laplace's Equation
$` \nabla \cdot \underline{V} = 0  \implies \nabla \cdot \nabla \phi = 0 \implies \nabla^2 \phi = 0`$

### Rotational Invariance of Laplace differential operator
A function defined on an inner product space is said to have rotational invariance if its value does not change when arbitrary rotations are applied to its argument. 
$` f(\underline{x}') = f(\mathbf{R}\underline{x}) = f(\underline{x}) `$ This also applies for an operator that acts on a function $` f : \mathbb{U} \subseteq \mathbb{R} \to \mathbb{U} `$.
In that case rotational invariance may also mean that the function commutes with rotations of elements in $` \mathbb{U} `$. An example is the Laplace differential operator 
$` \nabla^2 = \frac{\partial^2()}{\partial x^2} + \frac{\partial^2()}{\partial y^2} + \frac{\partial^2()}{\partial z^2} `$.


$` \nabla_\mathbf{x}^2 u(\mathbf{x}) = \nabla_\mathbf{y}^2 v(\mathbf{y}) `$, where $` \mathbf{y} = \mathbf{R} \mathbf{x} `$ 
and $`  u(\mathbf{x}) = v(\mathbf{y}) = v(\mathbf{R} \mathbf{x}) `$

### Fundamental Solution of the Laplace operator
Since  Laplace equation is invariant under rigid motions, it is natural to look for solutions to $` \nabla^2 \psi = 0 `$ which have rotational symmetry and the form 
$` \psi \colon \mathbb{R}^3 \to \mathbb{R} \colon \psi(\underline{r},\underline{r}_p) = \psi(\lVert \underline{r} - \underline{r}_p \rVert) `$

Assuming that $` \underline{r} \neq \underline{r}_p `$ , it is true that

```math
\begin{align*}
&\frac{\partial \psi}{\partial x} = \psi'(\lVert \underline{r} - \underline{r}_p \rVert) \frac{x-x_p}{\lVert \underline{r} - \underline{r}_p \rVert} \, , \qquad
\frac{\partial \psi}{\partial x} = \psi'(\lVert \underline{r} - \underline{r}_p \rVert) \frac{x-x_p}{\lVert \underline{r} - \underline{r}_p \rVert} \, , \qquad
\frac{\partial \psi}{\partial x} = \psi'(\lVert \underline{r} - \underline{r}_p \rVert) \frac{x-x_p}{\lVert \underline{r} - \underline{r}_p \rVert}
\\
&\frac{\partial^2 \psi}{\partial x^2} = \frac{(x-x_p)^2}{\lVert \underline{r} - \underline{r}_p \rVert^2} \psi''(\lVert \underline{r} - \underline{r}_p \rVert)
+
\frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \psi'(\lVert \underline{r} - \underline{r}_p \rVert)
-
\frac{(x-x_p)^2}{\lVert \underline{r} - \underline{r}_p \rVert^3} \psi'(\lVert \underline{r} - \underline{r}_p \rVert)
\\
&\frac{\partial^2 \psi}{\partial y^2} = \frac{(y-y_p)^2}{\lVert \underline{r} - \underline{r}_p \rVert^2} \psi''(\lVert \underline{r} - \underline{r}_p \rVert)
+
\frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \psi'(\lVert \underline{r} - \underline{r}_p \rVert)
-
\frac{(y-y_p)^2}{\lVert \underline{r} - \underline{r}_p \rVert^3} \psi'(\lVert \underline{r} - \underline{r}_p \rVert)
\\
&\frac{\partial^2 \psi}{\partial z^2} = \frac{(z-z_p)^2}{\lVert \underline{r} - \underline{r}_p \rVert^2} \psi''(\lVert \underline{r} - \underline{r}_p \rVert)
+
\frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \psi'(\lVert \underline{r} - \underline{r}_p \rVert)
-
\frac{(z-z_p)^2}{\lVert \underline{r} - \underline{r}_p \rVert^3} \psi'(\lVert \underline{r} - \underline{r}_p \rVert)
\end{align*}
```

```math
\nabla^2 \psi = \psi''(\lVert \underline{r} - \underline{r}_p \rVert) - \frac{2}{\lVert \underline{r} - \underline{r}_p \rVert} \psi'(\lVert \underline{r} - \underline{r}_p \rVert) = 0
\qquad \forall \underline{r} \in \mathbb{R}^3 - \{\underline{r}_p\}
```

```math
\begin{align*}
&\psi''(\lVert \underline{r} - \underline{r}_p \rVert) - \frac{2}{\lVert \underline{r} - \underline{r}_p \rVert} \psi'(\lVert \underline{r} - \underline{r}_p \rVert) = 0 \implies
\\
&\lVert \underline{r} - \underline{r}_p \rVert^2 \psi''(\lVert \underline{r} - \underline{r}_p \rVert)
- 2 \lVert \underline{r} - \underline{r}_p \rVert \psi'(\lVert \underline{r} - \underline{r}_p \rVert)
= 0
\implies
\left[ \lVert \underline{r} - \underline{r}_p \rVert^2 \psi'(\lVert \underline{r} - \underline{r}_p \rVert) \right]' = 0
\implies
\\
&\psi'(\lVert \underline{r} - \underline{r}_p \rVert) =
\frac{c}{\lVert \underline{r} - \underline{r}_p \rVert^2}
\implies
\psi(\lVert \underline{r} - \underline{r}_p \rVert) = - \frac{c}{\lVert \underline{r} - \underline{r}_p \rVert} + c' \, , \qquad c, c' \in \mathbb{R}
\end{align*}
```


When  $` c = \frac{1}{4 \pi} `$ and $` c' = 0 `$ are chosen, so that
$` \psi(\lVert \underline{r} - \underline{r}_p \rVert) = - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} `$, it can be shown that 

```math
\nabla^2 \psi(\underline{r}, \underline{r}_p) = \delta(\underline{r} - \underline{r}_p)
\qquad \text{and} \qquad \iiint_V f(\underline{r}) \delta(\underline{r} - \underline{r}_p) \mathrm{d}V = f(\underline{r}_p)
```

where $` \delta(\underline{r} - \underline{r}_p) `$ is the Dirac delta function defined in $` \mathrm{C}_c^\infty(V) `$ and, $` f \colon V \subseteq \mathbb{R}^3 \to \mathbb{R} `$ and $` f \in \mathrm{C}_c^\infty(V) `$.


$` \psi(\lVert \underline{r} - \underline{r}_p \rVert) = - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} `$ is called Fundamental Solution of the Laplace operator


### Gauss's Divergence Theorem

Let $` V \subset \mathbb{R}^3 `$ be a bounded domain, and his boundary $` \partial V `$.
Let $` \partial V `$ be a smooth hypersurface and $` \underline{n} `$ the outward unit normal vector to $` \partial V `$.
Supose $` \underline{F} \colon V \subseteq \mathbb{R}^3 \to \mathbb{R}^3 `$ and  $` F \in \mathrm{C}^1(V) \cap \mathrm{C}^0(\partial V) `$. It is true that

```math
\iiint_V \nabla \cdot \underline{F} \mathrm{d}V = \iint_{\partial V} \underline{F} \cdot \underline{n} \mathrm{d}S
```

### Green's 2nd Identity
Let $` V \subset \mathbb{R}^3 `$ be a bounded domain, and his boundary $` \partial V `$.
Let $` \partial V `$ be a smooth hypersurface and $` \underline{n} `$ the outward unit normal vector to $` \partial V `$.
Supose $` \underline{F} = \psi \nabla \phi - \phi \nabla \psi `$, 
where $` \psi \, , \phi \in \mathrm{C}^2(V) \cap \mathrm{C}^1(\partial V) `$. It is true that

```math
\begin{align*}
&\iiint_V \nabla \cdot \underline{F} \mathrm{d}V = \iint_{S} \underline{F} \cdot \underline{n} \mathrm{d}S
\implies
\iiint_V \nabla \cdot \left( \psi \nabla \phi - \phi \nabla \psi \right) \mathrm{d}V = 
\iint_{\partial V} \left( \psi \nabla \phi - \phi \nabla \psi \right) \cdot \underline{n} \mathrm{d}S
\implies
\\
&\iiint_V \left( \nabla \psi \cdot \nabla \phi + \psi \nabla^2 \phi  - \nabla \phi \cdot \nabla \psi - \phi \nabla^2 \psi \right) \mathrm{d}V = 
\iint_{\partial V} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
\implies
\end{align*}
```
```math
\begin{aligned}
\iiint_V \left(\psi \nabla^2 \phi - \phi \nabla^2 \psi \right) \mathrm{d}V &=
\iint_{\partial V} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
\\
&= \iint_{\partial V} \left( \psi  \frac{\partial \phi}{\partial n} - \phi \frac{\partial \psi}{\partial n} \right)\mathrm{d}S
\end{aligned}
```



### Integral Equation of velocity potential $\phi$
* Let $` V \subset \mathbb{R}^3 `$ be a bounded domain, and his boundary $` \partial V `$.
* Let $` \partial V =  S_\infty \cup S \cup S_w `$ be a smooth hypersurface and $` \underline{n} (= - \underline{e}_n ) `$ the outward unit normal vector to $` \partial V `$.
* Let $` \phi \in \mathrm{C}^2(V) \cap \mathrm{C}^1(\partial V) `$ and $` \psi(\lVert \underline{r} - \underline{r}_p \rVert) = - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} `$
* Let $` V_\epsilon = V - B[\underline{r}_p, \epsilon] `$. Then $` \partial V_\epsilon = \partial V \cup \partial B[\underline{r}_p, \epsilon] = S_\infty \cup S \cup S_w \cup S_\epsilon `$

Using Green's 2nd Identity we have
```math
\begin{align*}
&\iiint_{V_\epsilon} \left(\psi \nabla^2 \phi - \phi \nabla^2 \psi \right) \mathrm{d}V =
\iint_{\partial V_\epsilon} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
\xRightarrow{\nabla^2 \psi(\underline{r}, \underline{r}_p) = 0 \, , \underline{r} \neq \underline{r}_p}
\iiint_{V_\epsilon} \psi \nabla^2 \phi \mathrm{d}V =
\iint_{\partial V_\epsilon} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
\xRightarrow[ {\partial V_\epsilon = \partial V \cup \partial B[\underline{r}_p, \epsilon]} ]{ {V_\epsilon = V - B[\underline{r}_p, \epsilon]} }
\\
&
\iiint_{V} \psi \nabla^2 \phi \mathrm{d}V - \iiint_{B[\underline{r}_p, \epsilon]} \psi \nabla^2 \phi \mathrm{d}V =
\iint_{\partial V} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
+ \iint_{\partial B[\underline{r}_p, \epsilon]} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
\xRightarrow{\epsilon \to 0}
\end{align*}
```
```math
\iiint_{V} \psi \nabla^2 \phi \mathrm{d}V - \lim_{\epsilon \to 0} \iiint_{B[\underline{r}_p, \epsilon]} \psi \nabla^2 \phi \mathrm{d}V =
\iint_{\partial V} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
+
\lim_{\epsilon \to 0} \iint_{\partial B[\underline{r}_p, \epsilon]} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
```

* $` \lim\limits_{\epsilon \to 0} \iiint_{B[\underline{r}_p, \epsilon]} \psi \nabla^2 \phi \mathrm{d}V = 0 `$
* $` \lim\limits_{\epsilon \to 0} \iint_{\partial B[\underline{r}_p, \epsilon]} \psi (\underline{n} \cdot \nabla) \phi  \mathrm{d}S = 0 `$ 
* $` \lim\limits_{\epsilon \to 0} \iint_{\partial B[\underline{r}_p, \epsilon]} \phi (\underline{n} \cdot \nabla) \psi \mathrm{d}S = - \phi(\underline{r}_p) `$ 

#### Important notes
* Since $` \phi \in C^2(B[\underline{r}_p, \epsilon]) `$ and $` B[\underline{r}_p, \epsilon] `$ is a compact subset of $` \mathbb{R}^3 `$, then $` \nabla^2 \phi `$ is bounded in $` B[\underline{r}_p, \epsilon] `$ $` \left( \exists \, M \in \mathbb{R}^3: \lvert \nabla^2 \phi(\underline{r}) \rvert \leq M \quad \forall \underline{r} \in B[\underline{r}_p, \epsilon] \right) `$
* $` \psi(\lVert \underline{r} - \underline{r}_p \rVert) = - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} = const \, , \qquad \forall \underline{r} \in \partial B[\underline{r}_p, \epsilon] `$
* $` \underline{n} = \frac{\underline{r}_p - \underline{r}}{\lVert \underline{r}_p - \underline{r} \rVert} = - \frac{\underline{r} - \underline{r}_p}{\lVert \underline{r} - \underline{r}_p \rVert} = - \underline{e}_n \, , \qquad \forall \underline{r} \in \partial B[\underline{r}_p, \epsilon]  `$

#### Proof


```math
\begin{flalign*}
\bullet \quad \lim_{\epsilon \to 0} \iiint_{B[\underline{r}_p, \epsilon]} \psi \nabla^2 \phi  \mathrm{d}V
&= \nabla^2 \phi \lim_{\epsilon \to 0} \iiint_{B[\underline{r}_p, \epsilon]} \psi \mathrm{d}V
= \nabla^2 \phi \lim_{\epsilon \to 0} \iiint_{B[\underline{r}_p, \epsilon]} - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \mathrm{d}V 
= - \frac{1}{4 \pi} \nabla^2 \phi \lim_{\epsilon \to 0} \iiint_{B[\underline{r}_p, \epsilon]} \frac{1}{\rho} \mathrm{d}V
&& \\
&= - \frac{1}{4 \pi} \nabla^2 \phi \lim_{\epsilon \to 0} \int_0^\epsilon \iint_{\partial B[ 0, 1]} \frac{1}{\rho} \rho^2 \mathrm{d}S \, \mathrm{d}\rho
= - \frac{1}{4 \pi} \nabla^2 \phi \lim_{\epsilon \to 0} \int_0^\epsilon \rho \iint_{\partial B[0, 1]} \mathrm{d}S \, \mathrm{d}\rho 
= - \frac{1}{4 \pi} \nabla^2 \phi \lim_{\epsilon \to 0} \int_0^\epsilon 4\pi \rho \, \mathrm{d}\rho
&& \\
&= -  \nabla^2 \phi \lim_{\epsilon \to 0} \int_0^\epsilon \rho \, \mathrm{d}\rho
= - \nabla^2 \phi \lim_{\epsilon \to 0} \int_0^\epsilon \frac{1}{2} \frac{\mathrm{d} \rho^2}{\mathrm{d} \rho} \, \mathrm{d}\rho
= - \frac{1}{2} \nabla^2 \phi \lim_{\epsilon \to 0} \int_0^\epsilon  \, \mathrm{d}\rho^2
= - \frac{1}{2} \nabla^2 \phi \lim_{\epsilon \to 0} \epsilon^2 = 0
\end{flalign*}
```

```math
\begin{flalign*}
\bullet \quad \lim_{\epsilon \to 0} \iint_{\partial B[\underline{r}_p, \epsilon]} \psi (\underline{n} \cdot \nabla) \phi \, \mathrm{d}S
&= \lim_{\epsilon \to 0} \iint_{\partial B[\underline{r}_p, \epsilon]} - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} (\underline{n} \cdot \nabla) \phi \, \mathrm{d}S
= \lim_{\epsilon \to 0} - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert}  \iint_{\partial B[\underline{r}_p, \epsilon]}  (\underline{n} \cdot \nabla) \phi \, \mathrm{d}S
&& \\
&= \lim_{\epsilon \to 0} - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert}  \iiint_{B[\underline{r}_p, \epsilon]} \nabla \cdot \nabla \phi \, \mathrm{d}V
= \lim_{\epsilon \to 0} - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert}  \iiint_{B[\underline{r}_p, \epsilon]} \nabla^2 \phi \, \mathrm{d}V
= \lim_{\substack{\epsilon \to 0 \\ \underline{r} \to \underline{r}_p}} - \frac{1}{4 \pi} \frac{\nabla^2 \phi(\underline{r})}{\lVert \underline{r} - \underline{r}_p \rVert}  \iiint_{B[\underline{r}_p, \epsilon]}  \mathrm{d}V
&& \\
&= \lim_{\substack{\epsilon \to 0 \\ \underline{r} \to \underline{r}_p}} - \frac{1}{4 \pi} \frac{\nabla^2 \phi(\underline{r})}{\lVert \underline{r} - \underline{r}_p \rVert}  \frac{4}{3} \pi \lVert \underline{r} - \underline{r}_p \rVert^3
= - \frac{1}{3} \nabla^2 \phi(\underline{r}_p) \lim_{\substack{\epsilon \to 0 \\ \underline{r} \to \underline{r}_p}} \lVert \underline{r} - \underline{r}_p \rVert^2 = 0
\end{flalign*}
```

```math
\begin{flalign*}
\bullet \quad \lim_{\epsilon \to 0} \iint_{\partial B[\underline{r}_p, \epsilon]} \phi (\underline{n} \cdot \nabla) \psi \mathrm{d}S
&= \lim_{\epsilon \to 0} \iint_{\partial B[\underline{r}_p, \epsilon]} \phi \left(- \frac{\underline{r} - \underline{r}_p}{\lVert \underline{r} - \underline{r}_p \rVert} \cdot \nabla \right) \left(- \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \right) \mathrm{d}S
= \frac{1}{4 \pi} \lim_{\substack{\epsilon \to 0 \\ \underline{r} \to \underline{r}_p}} \iint_{\partial B[\underline{r}_p, \epsilon]} \phi  \left( - \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert^2} \right) \mathrm{d}S
&& \\
&= - \frac{1}{4 \pi} \lim_{\substack{\epsilon \to 0 \\ \underline{r} \to \underline{r}_p}} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert^2} \phi(\underline{r})  \iint_{\partial B[\underline{r}_p, \epsilon]} \mathrm{d}S
= - \frac{1}{4 \pi} \lim_{\substack{\epsilon \to 0 \\ \underline{r} \to \underline{r}_p}} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert}^2 \phi(\underline{r}) 4 \pi \lVert \underline{r} - \underline{r}_p \rVert^2
= - \phi(\underline{r}_p)
\end{flalign*}
```

By substituting the values of the limits into the integral equation, we have
```math
\iiint_{V} \psi \nabla^2 \phi \mathrm{d}V =
\iint_{\partial V} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
+
\phi(\underline{r}_p)
```

### Internal Potential $` \phi_i `$
* Let $` V_i \subset \mathbb{R}^3 `$ be a bounded domain, and his boundary $` \partial V_i `$.
* Let $` \partial V_i =  S \cup S_w `$ be a smooth hypersurface and $` \underline{n_i} (= \underline{e}_n ) `$ the outward unit normal vector to $` \partial V_i `$.
* Let $` \phi_i \in \mathrm{C}^2(V_i) \cap \mathrm{C}^1(\partial V_i) `$ and $` \psi(\lVert \underline{r} - \underline{r}_p \rVert) = - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} `$
* $` \underline{r}_p \notin V_i \quad ( \underline{r}_p \in V ) \implies \underline{r} \neq \underline{r}_p `$


Using Green's 2nd Identity we have
```math
\iiint_{V_i} \left(\psi \nabla^2 \phi_i - \phi_i \nabla^2 \psi \right) \mathrm{d}V =
\iint_{\partial V_i} \left[ \psi (\underline{n_i} \cdot \nabla) \phi_i - \phi_i (\underline{n_i} \cdot \nabla) \psi \right] \mathrm{d}S
\xRightarrow{\nabla^2 \psi(\underline{r}, \underline{r}_p) = 0 \, , \underline{r} \neq \underline{r}_p}
```

```math
\iiint_{V_i} \psi \nabla^2 \phi_i \mathrm{d}V =
\iint_{\partial V_i} \left[ \psi (\underline{n_i} \cdot \nabla) \phi_i - \phi_i (\underline{n_i} \cdot \nabla) \psi \right] \mathrm{d}S
```

* If $` phi `$ is the velocity potential of an conservative and incompressible velocity vector field $` V `$ then it it satisfies Laplace's equation $` \left( \nabla^2 \phi = 0 \right) `$
* If $` \phi_i `$ satisfies also Laplace's equation $` \left( \nabla^2 \phi_i = 0 \right) `$

Then combining these two integral equations we have

```math
\begin{align*}
&\iiint_{V} \psi \nabla^2 \phi \mathrm{d}V + \iiint_{V_i} \psi \nabla^2 \phi_i \mathrm{d}V =
\iint_{\partial V} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
+
\phi(\underline{r}_p)
+
\iint_{\partial V_i} \left[ \psi (\underline{n_i} \cdot \nabla) \phi_i - \phi_i (\underline{n_i} \cdot \nabla) \psi \right] \mathrm{d}S
\xRightarrow{\nabla^2 \phi = \nabla^2 \phi_i = 0}
\\
&\phi(\underline{r}_p)
+ \iint_{\partial V} \left[ \psi (\underline{n} \cdot \nabla) \phi - \phi (\underline{n} \cdot \nabla) \psi \right] \mathrm{d}S
+ \iint_{\partial V_i} \left[ \psi (\underline{n_i} \cdot \nabla) \phi_i - \phi_i (\underline{n_i} \cdot \nabla) \psi \right] \mathrm{d}S = 0
\xRightarrow{\underline{e}_n = - \underline{n} = \underline{n_i}}
\\
&\phi(\underline{r}_p)
=
\iint_{\partial V} \left[ \psi (\underline{e}_n \cdot \nabla) \phi - \phi (\underline{e}_n \cdot \nabla) \psi \right] \mathrm{d}S 
- \iint_{\partial V_i} \left[ \psi (\underline{e}_n \cdot \nabla) \phi_i - \phi_i (\underline{e}_n \cdot \nabla) \psi \right] \mathrm{d}S
\xRightarrow[\partial V_i =  S \cup S_w]{\partial V =  S_\infty \cup S \cup S_w}
\\
&\phi(\underline{r}_p)
=
\iint_{S \cup S_w} \left[ \psi (\underline{e}_n \cdot \nabla) (\phi - \phi_i) - (\phi - \phi_i) (\underline{e}_n \cdot \nabla) \psi \right] \mathrm{d}S
+ \iint_{S_\infty} \left[ \psi (\underline{e}_n \cdot \nabla) \phi - \phi (\underline{e}_n \cdot \nabla) \psi \right] \mathrm{d}S
\implies
\end{align*}
```

```math
\begin{align*}
&\phi(\underline{r}_p)
=
\iint_{S \cup S_w} \left[ \psi (\underline{e}_n \cdot \nabla) (\phi - \phi_i) - (\phi - \phi_i) (\underline{e}_n \cdot \nabla) \psi \right] \mathrm{d}S
+ \phi_\infty(\underline{r}_p)
\\
&\text{where} \quad \psi(\lVert \underline{r} - \underline{r}_p \rVert) = - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert}
\quad \text{and} \quad \phi_\infty(\underline{r}_p) = \iint_{S_\infty} \left[ \psi (\underline{e}_n \cdot \nabla) \phi - \phi (\underline{e}_n \cdot \nabla) \psi \right] \mathrm{d}S
\end{align*} 
```



#### Wake Assumption
By considering that the volume of space enclosed by the wake surface $` S_w `$ is very small, we can assume that $` (\underline{e}_n \cdot \nabla) (\phi - \phi_i) = 0 \quad \forall \underline{r} \in  S_w `$. 


```math
\begin{align*}
\phi(\underline{r}_p)
&=
\iint_S \left[ \psi (\underline{e}_n \cdot \nabla) (\phi - \phi_i) - (\phi - \phi_i) (\underline{e}_n \cdot \nabla) \psi \right] \mathrm{d}S
+ \iint_{S_w} - (\phi - \phi_i) (\underline{e}_n \cdot \nabla) \psi \mathrm{d}S
+ \phi_\infty(\underline{r}_p)
\implies
\\
\phi(\underline{r}_p)
&=
\iint_S \left[ \frac{1}{4 \pi} (\phi - \phi_i) (\underline{e}_n \cdot \nabla) \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} (\underline{e}_n \cdot \nabla) (\phi - \phi_i) \right] \mathrm{d}S
+ \iint_{S_w} \frac{1}{4 \pi} (\phi - \phi_i) (\underline{e}_n \cdot \nabla) \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \mathrm{d}S
+ \phi_\infty(\underline{r}_p)
\end{align*}
```


#### Velocity Potential $` \phi(\underline{r}_p) `$
```math
\begin{align*}
\phi(\underline{r}_p)
&=
\iint_S \left[ \frac{\mu}{4 \pi} (\underline{e}_n \cdot \nabla) \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} - \frac{\sigma}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \right] \mathrm{d}S
+ \iint_{S_w}  \frac{\mu}{4 \pi} (\underline{e}_n \cdot \nabla) \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \mathrm{d}S
+ \phi_\infty(\underline{r}_p)
\\
\phi(\underline{r}_p)
&=
\iint_S - \frac{\sigma}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \mathrm{d}S
+ \iint_{S \cup S_w}  \frac{\mu}{4 \pi} (\underline{e}_n \cdot \nabla) \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \mathrm{d}S
+ \phi_\infty(\underline{r}_p) 
\end{align*}
```

where:
   * $` \mu = \phi - \phi_i `$
   * $` \sigma = (\underline{e}_n \cdot \nabla)(\phi - \phi_i) `$
   * $` \phi_\infty(\underline{r}_p) = \iint_{S_\infty} \left[ \phi (\underline{e}_n \cdot \nabla) - \frac{1}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} (\underline{e}_n \cdot \nabla) \phi \right] \mathrm{d}S `$


if $` \mu = \phi - \phi_i = \phi - \phi_\infty `$ and $` \sigma = (\underline{e}_n \cdot \nabla)(\phi - \phi_i) = (\underline{e}_n \cdot \nabla)(\phi - \phi_\infty ) `$, then for $` P \in S^- `$ :

```math
\iint_S - \frac{\sigma}{4 \pi} \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \mathrm{d}S
+ \iint_{S \cup S_w}  \frac{\mu}{4 \pi} (\underline{e}_n \cdot \nabla) \frac{1}{\lVert \underline{r} - \underline{r}_p \rVert} \mathrm{d}S = 0, \qquad \forall (x_p, y_p, z_p) \in S: (\underline{r} - \underline{r}_p) \cdot \underline{e}_n \to 0
```

note: $` \sigma = (\underline{e}_n \cdot \nabla)(\phi - \phi_i) = (\underline{e}_n \cdot \nabla)(\phi - \phi_\infty ) = (\underline{e}_n \cdot \nabla)\phi - (\underline{e}_n \cdot \nabla)\phi_\infty = \underline{V} - \underline{e}_n \cdot \underline{V}_\infty = - \underline{e}_n \cdot \underline{V}_\infty `$


### Numerical Model (Panel Methods)

```math
\sum_{j=0}^{N_s - 1} B_{ij} \sigma_j + \sum_{j=0}^{N_s + N_w - 1}C_{ij} \mu_j = 0 , \qquad 0 \le i < N_s 
``` 

where:
   * $` B_{ij} = - \frac{1}{4 \pi} \iint_{S_j}  \frac{1}{\lVert \underline{r} - \underline{r}_{cp_i} \rVert} \mathrm{d}{S_j}  = - \frac{1}{4\pi} \iint_{S_j} \frac{1}{\sqrt{(l_{j} - l_{j_{cp_i}})^2 + (m_{j} - m_{j_{cp_i}})^2 + (n_{j} - n_{j_{cp_i}})^2}}  \mathrm{d}S_j = - \frac{1}{4\pi} \iint_{S_j} \frac{1}{\sqrt{(l_{j} - l_{j_{cp_i}})^2 + (m_{j} - m_{j_{cp_i}})^2 + n_{j_{cp_i}}^2}}  \mathrm{d}S_j `$

   * $` C_{ij} =  \frac{1}{4\pi} \iint_{S_j}  (\underline{e}_n \cdot \nabla) \frac{1}{\lVert \underline{r} - \underline{r}_{cp_i} \rVert} \mathrm{d}{S_j} = - \frac{1}{4\pi} \iint_{S_j}
\frac{n_{j} - n_{j_{cp_i}}}{\left( \sqrt{(l_{j} - l_{j_{cp_i}})^2 + (m_{j} - m_{j_{cp_i}})^2 + (n_{j} - n_{j_{cp_i}})^2} \right)^3} \mathrm{d}S_j =
\frac{1}{4\pi} \iint_{S_j}
\frac{n_{j_{cp_i}}}{\left( \sqrt{(l_{j} - l_{j_{cp_i}})^2 + (m_{j} - m_{j_{cp_i}})^2 + n_{j_{cp_i}}^2} \right)^3} \mathrm{d}S_j `$
   
   * $` \sigma_j = - \underline{e}_{n_j} \cdot \underline{V}_\infty `$

#### Steady Panel Method
from Kutta Condition: $` \mu_w = const = \mu_U - \mu_L `$
```math
A_{ij} \mu_j = - B_{ij} \sigma_j , \qquad A_{ij} = 
\begin{cases}
   C_{ij} + \sum\limits_{k=0}^{NWP_j} C_{if_j(k)} & \text{if the $j$-th panel is located on the upper side of the surface and adjoins the trailing edge}\\
   C_{ij} & \text{if the $j$-th panel do not adjoins the trailing edge}\\
   C_{ij} - \sum\limits_{k=0}^{NWP_j} C_{if_j(k)} & \text{if the $j$-th panel is located on the lower side of the surface and adjoins the trailing edge}
\end{cases} 
```

```math
0 \le i < N_s  \qquad 0 \le j < N_s  \qquad 0 \le k < NWP_j  
```
where $`NWP_j`$ is the number of panels that the wake row shedding from the $`j`$-th panel consists of, and $` f_j(k) `$ returns the id of the $`k`$-th panel of the wake row shedding from $`j`$-th panel

### Features
 1. Calculation of Non-Lifting Potential Flow about 3D arbitrarily-shaped rigid bodies
    1. Steady simulations
    2. Unsteady simulations

 2. Calculation of Lifting Pseudo-Potential Flow around 3D arbitrarily-shaped rigid bodies
      1. Steady state simulations with flat rigid wake model
      2. Steady state iterative simulations with flexible wake model 
      3. Unsteady simulations with a shedding wake model

## Simulation Results
