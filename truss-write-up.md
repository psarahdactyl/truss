title: Truss Optimization as a Linear Program

<script type="text/x-mathjax-config">
    MathJax.Hub.Config({ TeX: { equationNumbers: {autoNumber: "all"} } });
</script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>

<div style="display:none">
$\newcommand{\a}{\mathbf{a}}$
$\newcommand{\B}{\mathbf{B}}$
$\newcommand{\b}{\hat{\mathbf{b}}}$
$\newcommand{\E}{\mathbf{E}}$
$\newcommand{\e}{\mathbf{e}}$
$\newcommand{\f}{\mathbf{f}}$
$\newcommand{\I}{\mathbf{I}}$
$\newcommand{\Zero}{\mathbf{0}}$
$\newcommand{\l}{\mathbf{l}}$
$\newcommand{\t}{\mathbf{t}}$
$\newcommand{\R}{\mathbb{R}}$
$\newcommand{\V}{\mathbf{V}}$
$\renewcommand{\v}{\mathbf{v}}$
</div>

# Truss Optimization as a Linear Program

We consider the problem of finding minimal volume structure that can support a
given loading condition.
<!---->
We limit our consideration to truss structures: rigid networks of beams.
<!---->
For us, the beams of our truss maybe be different thicknesses, but each
individual beam must be straight and have a uniform thickness (i.e., constant
cross-section along the beam).
<!---->
We assume that beams are strong under compression and tension, but weak under
bending.
<!---->
Intuitively, we will imagine that beams perfectly resist compression (or
tension) until a breaking point, but would instantly fail in the presence of
forces causing them to bend.

A truss is a graph embedded $\R^d$ represented as a list of $n$ vertex locations
$\V = \{\v_1,\v_2,\dots, \v_n\} ∈ \R^{n×d}$, a set of $m$ edges $\E =
\{\e_1,\e_2,\dots, \e_m\} ∈ [1,n]^{m×2}$, and a list of non-negative
cross-sectional areas associated with each edge $\a = \{a_1,a_2,\dots,a_m\} ∈
\R_{≥0}^m$.

Treating the cross-sectional areas $\a$ as
unknowns, 
our truss optimization problem is to find the minimal volume truss that can
support a given set of external loads (forces) specified at vertices $\f
= \mathop{\text{vec}}(\{\f_1,\f_2,\dots,\f_n\}) ∈ \R^{dn}$.
<!---->
To pose this problem, we introduce a scalar _unknown_ per-edge to represent the tensile
force exerted on each of its endpoint vertices in the outward direction of the
edge: $\t = \{t_1,t_2,\dots,t_m\} ∈ \R^m$.
<!---->
If $i$ th edge $\e_i = \{j,k\}$ is under tension, then its tensile force is positive ($t_i > 0$), and
the edge _pulls_ inwardly on the vertex at $\v_j$ (in the
direction of $\b_i = (\v_k-\v_j)/\|\v_k-\v_j\| ∈ \R^d$) and on $\v_k$ (in the
opposite direction, $-\b_i$).
<!---->
If the tensile force is negative ($t_i<0$), then the edge is under compression,
and pulling becomes pushing.
<!---->
To account for the material failure of the beams, we assume we know the tensile
and compressive strength of the material given as a maximum ratio of tensile
force to cross-sectional area: $σ_t ≥ t_i/a_i$ and $σ_c ≥ -t_i/a_i$,
respectively. Without loss of generality, we'll assume these strengths to be the
same for all beams (i.e., the truss is made out of one material).
<!---->
Finally, let us introduce the list of edge lengths $\l = \{l_1,l_2,\dots,l_m\}
∈\R^m$, where $l_i = \|\v_k-\v_j\|$ for the edge $\e_i = \{j,k\}$.

With these definitions and variables defined, we can write the truss
optimization problem as minimization of total volume subject to $d\,n$ linear force
balance equalities and $2m$ linear inequalities imposing strength limits:

$$
\mathop{\text{minimize}}_{\a,\t} \sum\limits_{i=1}\^m a_i l_i \\\\
\text{ subject to } \sum\limits_{\e_i=\\{j,k\\}} t_i \b_i - \sum\limits_{\e_i=\\{k,j\\}} t_i \b_i = \f_j, \ ∀ j = 1\dots n  \\\\
%\text{ subject to } \f_j = \sum\limits_{i=1}\^m \begin{cases}
%  t_i \b_i & \text{if } \e_i = \\{j,k\\},\\\\
%  -t_i \b_i & \text{else if } \e_i = \\{k,j\\},\\\\
%  0         & \text{else.}
%\end{cases}, \ ∀ j = 1\dots n \\\\
\text{ and } σ_t a_i ≥ t_i, \ ∀ i = 1\dots m \\\\
\text{ and } σ_c a_i ≥ -t_i, \ ∀ i = 1\dots m
$$

Or gathering terms into appropriately sized matrices:
$$
\mathop{\text{minimize}}_{\begin{bmatrix}\a \\\\ \t \end{bmatrix}} \  \begin{bmatrix}\l^\top & \Zero \end{bmatrix} \begin{bmatrix}\a \\\\ \t \end{bmatrix} \\\\
\text{ subject to } \begin{bmatrix}\Zero & \B \end{bmatrix} \begin{bmatrix}\a \\\\ \t \end{bmatrix} = \f \\\\
\text{ and } \\begin{bmatrix} σ_t \I & -\I \\\\  σ_c \I & \I\end{bmatrix} \begin{bmatrix}\a \\\\ \t \end{bmatrix}≥ \Zero,
$$
where $\B ∈ \R^{nd × m}$ is a sparse matrix collecting the appropriately signed
edge directions and $\I ∈ \R^{m×m}$ is an identity matrix.

**Note on Constraints:**
Boundary conditions can be incoporated by removing or adjusting the rows of the
linear equality conditions ($\B\ \t = \f$). In particular, to _fix_ (_pin_) the
position of $j$ th vertex, the corresponding $d$ rows in $\B$ and $\f$ should be
_removed_.

This is a linear program and can be immediately solved using the `linprog`
function of libraries such as MATLAB or Mosek.

 - inequalities + linear term are a type of lasso, sparsity inducing
   combination: at the optimum most elements in $\a$ will be exactly zero.
 - not just low volume, but low number of beams.
 - initialize the topology of the truss with a dense, over-connected graph
 - careful not to include co-linearly overlapping edges


[][#zegard2015]

Could rewrite with total volume as a _constraint_ and treat minimal _complaince_
as an objective [][#freund2004].

Ground Structure

  - linprog of [a;n] write in Aeq, Aieq form...
  - units!
  - B: edge unit vectors, orientation matters

Pitfalls

[#freund2004]: Robert Freund. "Truss Design and Convex Optimization", 2004.
[#zegard2015]: Tomás Zegard, Glaucio H. Paulino. "GRAND3 — Ground structure
based topology optimization for arbitrary 3D domains using MATLAB", 2015.
