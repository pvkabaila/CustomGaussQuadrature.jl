# CustomGaussQuadrature

[![Build Status](https://github.com/pvkabaila/CustomGaussQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/pvkabaila/CustomGaussQuadrature.jl/actions/workflows/CI.yml?query=branch%3Amaster)

Suppose that we wish to evaluate

$$\int_{-\infty}^{\infty} g(x) f(x) dx, $$

where $f$ is a specified nonnegative integrable
weight function. The $n$-point Gauss quadrature approximation to this
integral is

$$\sum_{i=1}^n \lambda_i \ g(\tau_i),$$

where $\tau_1, \dots, \tau_n$ are called the nodes
and $\lambda_1, \dots, \lambda_n$ are called the
corresponding weights. The dependence of these 
nodes and weights on $n$ is implicit.
This approximation is exact
whenever $g$ is a polynomial of degree less
than or equal to $2n - 1$. 
We call these nodes and weights, taken together, as the Gauss rule with $n$ nodes.

For the case that $f$ leads to Gauss 
rules with nodes that are the roots of classical orthogonal
polynomials of a continuous variable then these rules,
such as Gauss Legendre, Gauss Hermite and Gauss Laguerre, are
readily accessible. We call this the 
classical case. In the non-classical case, however, the Gauss rule needs to be custom-made.

Gauss rules are usually
computed in the following two steps.

**Step 1**: Compute the recursion coefficients in the three-step recurrence 
relation, given by (1) and (2) below.

**Step 2**: Use these recursion coefficients to compute the Gauss quadrature nodes 
and weights.

In the classical case there are simple formulae for the recursion coefficients, making **Step 1** trivial. In the non-classical case Step 1 is difficult. The treatise of Gautschi (2004) describes several methods for this computation that differ widely in both complexity and sensitivity to roundoff errors. However, **Step 2** remains the same for both classical and non-classical cases. We carry out **Step 2** by computing the eigenvalues and eigenvectors of a symmetric tridiagonal matrix using the package `GenericLinearAlgebra.jl`.

## Installation

```julia
using Pkg
Pkg.add("CustomGaussQuadrature")
```

Then load the package:

```julia
using CustomGaussQuadrature
```

# Step 1 using either moment determinants or the Stieltjes procedure

## **Three-term recurrence relation** 


The norm of the function $u$ is denoted by $\lVert u \rVert$
and is defined to be the square root of 

$$\int_{-\infty}^{\infty} u^2(x) f(x) dx,$$

provided that this integral exists. For two functions $u$ and $v$ that satisfy $\lVert u \rVert <\infty$ and $\lVert v \rVert <\infty$, the inner product of these two functions is denoted $(u,v)$ and is defined to be

$$\int_{-\infty}^{\infty} u(x) v(x) f(x) dx.$$

If $(u,v) = 0$ then the functions $u$ and $v$ are said to be orthogonal. 

We consider only polynomials whose coefficients are real numbers. 
A polynomial in $x$ of degree $n$ is said to be monic if the coefficient of $x^n$ is 1. Let $\pi_k$ denote a monic polynomial of degree $k$. The monic polynomials $\pi_0, \pi_1, \pi_2, \ldots$ are called monic orthogonal polynomials with respect to the weight function $f$ if 

$$(\pi_k, \pi_{\ell}) = 0 \ \ \text{for all} \ k \ne \ell, \ \text{where} \ k, \ell \in \{0, 1, 2, \ldots\}$$

and 

$$\lVert \pi_k \rVert > 0 \ \ \text{for} \ k=0, 1, 2, \ldots.$$

The Gauss quadrature nodes $\tau_1, \ldots, \tau_n$ are the $n$ distinct roots of the polynomial $\pi_n$. 

The monic orthogonal polynomials with respect to the weight function $f$ satisfy the following three-term recurrence relation, see e.g.
Theorem 1.27 on p.10 of Gautschi (2004). Let 
$\pi_{-1} \equiv 0$ and $\pi_0 \equiv 1$. Then 

$$
\pi_{k+1}(x) = 	(x - \alpha_k)  \pi_k(x) - \beta_k  \pi_{k-1}(x) \ \ \text{for} \ k = 0, 1, 2, \dots, 
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (1)
$$ 

where 

$$
\alpha_k 
= \frac{\big(x \, \pi_k, \pi_k\big)}{\big(\pi_k, \pi_k\big)}
(k = 0, 1, \dots) \ \ \text{and} \ \ \
\beta_k 
= \frac{\big(\pi_k, \pi_k\big)}{\big(\pi_{k-1}, \pi_{k-1}\big)}
(k = 1, 2, \dots) 
\ \ \ \ \ \ (2)
$$


## **Computation of the recursion coefficients in the three-term recurrence relation using moment determinants** 

We carry out **Step 1** using moment determinants (Theorem 2.2 of Gautschi, 2004).
The map from the vector 
$\big(\mu_0, \mu_1, \dots, \mu_{2n-1} \big)$ of moments to the vector $\big(\alpha_0, \dots, \alpha_{n-1},  \beta_1, \dots, \beta_{n-1}\big)$ of recursion coefficients is severely ill-conditioned (Gautschi, 1983, 2004). However, it has long been recognized that this limitation can be overcome by the use of high-precision arithmetic, see e.g. Gautschi (1983), when applying 
 the method of moment determinants.

We use `BigFloat` arithmetic, with number of bits of precision `b`. For assurance that the chosen number of bits of precision `b` is sufficiently large, we use the following simple method. Suppose that $c_1$ and $c_2$ are numerical approximations to the same quantity computed using the same Julia function
that implements an ill-conditioned numerical method, 
using `BigFloat` arithmetic with number of bits of precision `b1` and `b2`, respectively, where `b2` is substantially larger than `b1`.
The approximation $c_2$ is therefore believed to be much more accurate than $c_1$. Then we take $c_2$ as our final approximation and assess 
the absolute error to be less than $|c_1 - c_2|$.


We require that there is a formula, which can be computed in `BigFloat`
arithmetic, for the $s$'th moment

$$
\mu_s = \int_{-\infty}^{\infty} x^s \ f(x) \ dx
$$

for all nonnegative integers $s \le 2 n - 1$. This formula 
is provided by the user, as illustrated in the Example 3 below.
The function `custom_gauss_quad_all_fn` is then used to compute the Gauss quadrature nodes and weights. 

The aim of the method implemented in `custom_gauss_quad_all_fn` is for (a) the *absolute errors* of each of the nodes to be less than
$10^{-17}$ and (b) the *relative errors* of each of the weights to be less than $10^{-17}$. In other words, the nodes and weights are computed with the purpose of being used in further extensive `Float64` computations.


We identify the weight function $f$ by the array 
`which_f` with components:\
(i) the name given to $f$ (a string),\
(ii) support interval of $f$ specified by a 2-vector of the endpoints and \
(iii) parameter vector (if any).

Two examples of this identification are the following:

### **Scaled chi pdf weight function**

Consider the "scaled chi pdf" weight function given by the probability density function (pdf)

$$
f(x) =
\begin{cases}
\dfrac{m^{m/2}}{\Gamma(m/2) \ 2^{(m/2) - 1}} \ x^{m-1} \, \exp\big(- m \ x^2 /2 \big) &\text{for } x > 0
\\
0 &\text{otherwise},
\end{cases}
$$

where $m$ is a positive integer parameter (the "degrees of freedom"). 
This weight function is identified by <br>
    `which_f = ["scaled chi pdf", [0,Inf], m]`, <br>
for some assigned value of the parameter `m`. 

### **Chemistry example weight function**

The "chemistry example" weight function is

$$
f(x) =
\begin{cases}
\exp(-x^3 / 3) &\ \text{for} \ x > 0
\\
0  &\ \text{otherwise.}
\end{cases}
$$

This weight function is considered by Gautschi (1983) and is identified by <br>
`which_f = ["chemistry example", [0,Inf]]`.

### **The inputs to `custom_gauss_quad_all_fn`**

Suppose that we want to compute the custom Gauss quadrature rule with n nodes for the weight function $f$ specified by `which_f`. The function `custom_gauss_quad_all_fn` does this by carrying out both 
**Step 1** and **Step 2**. It has the following inputs:

- moment_fn, <br /> 
which has inputs `T` (floating-point type), `which_f` and `s` (moment order). <br /> It computes the `s`'th moment for the weight function specified by `which_f`.

- `which_f` identifies the weight function $f$. 

- `n` is the number of Gauss quadrature nodes. 

- `upto_n` is a Boolean variable with default value `false`

- `extra_check` is a Boolean variable with default value `false`

For any of the following which_f's 

which_f = ["scaled chi pdf", [0,Inf], m]\
which_f = ["chemistry example", [0, Inf]]\
which_f = ["Hermite", [-Inf, Inf]]\
which_f = ["Generalized Laguerre", [0, Inf], α_GGL]\
which_f = ["Legendre", [-1, 1]],

we proceed as illustrated by the following two examples.

<!--say which_f = ["scaled chi pdf", [0,Inf], m],
we assign the value of positive integer parameter m and then use
```julia
	which_f = ["scaled chi pdf", [0,Inf], m]
	moment_fn = moment_stored_fn
-->


### **A note on the input of real-valued parameters**

Suppose that the exact value of the parameter k is 3.1. For convenience, we assign it using

```julia
k = 3.1
```

which sets k to the `Float64` approximation to 3.1. This package uses a range of `AbstractFloat` types `T`, including `BigFloat`, for its internal computations. It is therefore important to obtain the best approximation of type `T` to the exact value 3.1. Since k is of type `Float64`, the command `convert(BigFloat, k)` gives only `Float64`-precision:

```
3.100000000000000088817841970012523233890533447265625
```

A more precise `BigFloat` approximation is obtained via `parse(BigFloat, string(k))`, which gives:

```
3.099999999999999999999999999999999999999999999999999999999999999999999999999986
```

This is because `string(k)` uses Julia's `print`, which for `Float64` values produces the shortest correctly rounded decimal that round-trips identically. Consequently `string(3.1)` returns `"3.1"`.

The function `parse(T, string(...))` is used internally within this package wherever a user-supplied parameter needs to be converted to type `T`. The user only needs to assign parameters using standard Julia literals (e.g. `k = 3.0`).

### **Example 1** 

We identify this weight function by first assigning a value of the positive integer parameter `m` and then using Julia command 
`which_f = ["scaled chi pdf", [0,Inf], m]`.
Therefore the following command specifies this weight function with `m` set to 160.
 ```julia
 	m = 160
	which_f = ["scaled chi pdf", [0,Inf], m]
```



For this weight function, the $s$'th moment is 

$$
\left(\frac{2}{m} \right)^{s/2}
\frac{\Gamma\big((s+m)/2\big)}{\Gamma(m/2)}
    = \exp \left(\frac{s}{2} \log\left(\frac{2}{m} \right)  +
    \log \Gamma \left(\frac{s+m}{2}\right) 
    - \log \Gamma \left(\frac{m}{2}\right) \right)
$$

for $s = 0, 1, 2, \dots$.
This formula has been implemented in the function `moment_stored_fn`.

The following commands compute the nodes and weights, as 
`Double64` vectors for a 5-point Gauss quadrature rule. 
```julia
	moment_fn = moment_stored_fn
    n = 5
    nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n)
```


This is followed by conversion to 
`Float64` vectors and printing using the `@printf` macro from the package `Printf`.
```julia
    using Printf
    nodes = convert(Vector{Float64}, nodes)
    weights = convert(Vector{Float64}, weights)
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.16e     " nodes[i]
        @printf "%.16e  \n" weights[i]
    end
```

### **Example 2**
Consider the weight function identified by
```julia
    which_f = ["chemistry example", [0, Inf]]
```

For this weight function, the $s$'th moment is 

$$
3^{(s - 2) / 3}  \ \Gamma\big((s + 1) / 3\big)
$$

for $s = 0, 1, 2, \dots$.
This formula has been implemented in the function `moment_stored_fn`.


The following commands compute the nodes and weights, as 
`Double64` vectors, for a 15-point Gauss quadrature rule. 

```julia
	moment_fn = moment_stored_fn
    n = 15
    nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n)
```

This is followed by printing these nodes and weights in the same format as Table 2.2 of 
Gautschi (1983), using the `@printf` macro from the package `Printf`.

```julia
    using Printf
    nodes = convert(Vector{BigFloat}, nodes)
    weights = convert(Vector{BigFloat}, weights)
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.15e     " nodes[i]
        @printf "%.15e  \n" weights[i]
    end
```

These results agree with Table 2.2 of 
Gautschi (1983).


### **Example 3** 

Now suppose that, on the other hand, we want to specify a **new** weight function, say the Weibull pdf with 
shape parameter k > 0 and scale parameter set to 1.
In this case, the weight function is
$$
f(x) =
\begin{cases}
\ k \, x^{k - 1} \exp(-x^{k}) &\ \text{for} \ \ x > 0
\\
0  &\ \text{otherwise.}
\end{cases}
$$
For this weight function, the $s$'th moment is
$$
\Gamma\left(1 + \frac{s}{k}\right)
$$
for $s = 0, 1, 2, \dots$.


We identify this new weight function by first assigning
a value of the positive parameter $k$ and then using Julia command
`which_f = ["weibull pdf", [0, Inf], k]`.
Therefore the following command specifies this weight function with $k$ set to 2.0.
```julia
which_f = ["weibull pdf", [0, Inf], 2.0]
```
We provide the function for computing the s'th moment using 
```julia
  using SpecialFunctions
  function moment_weibull_pdf_fn(::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    k = which_f[3]
    @assert k > 0
    T_k = parse(T, string(k))
    @assert s ≥ 0
    if s == 0
      return(convert(T, 1))
    end
    T_1 = convert(T, 1)
    T_s = convert(T, s)
    gamma(T_1 + (T_s/T_k))
  end
```

The following commands compute the nodes and weights, as 
`Double64` vectors for a 9-point Gauss quadrature rule. 
```julia
	moment_fn = moment_weibull_pdf_fn
    n = 9
    nodes, weights = custom_gauss_quad_all_fn(moment_fn, which_f, n)
```


This is followed by conversion to 
`Float64` vectors using the commands
```julia
    nodes = convert(Vector{Float64}, nodes)
    weights = convert(Vector{Float64}, weights)
```

These nodes and weights are printed using the `@printf` macro from the package `Printf` as follows.
```julia
    using Printf
    nodes = convert(Vector{Float64}, nodes)
    weights = convert(Vector{Float64}, weights)
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.16e     " nodes[i]
        @printf "%.16e  \n" weights[i]
    end
```

## **Computation of the recursion coefficients in the three-term recurrence relation using the Stieltjes procedure**

The Stieltjes procedure, given in subsections 2.2.2 and 2.2.3 of Gautschi (2004), applies the three-term recurrence relation (1) and (2) iteratively through the following sequence:

$$\{\pi_{-1}, \pi_0\} \rightarrow
\{\alpha_0, \pi_1\} \rightarrow
(\pi_1, \pi_1) \rightarrow
\beta_1  \rightarrow
\{\alpha_1, \pi_2\} \rightarrow
(\pi_2, \pi_2) \rightarrow
\beta_2  \rightarrow \dots$$

The inner products used in the computation of the recursion coefficients are found using a high-quality quadrature rule with $r$ nodes. This is a discretization method that is expected to lead to
approximations to the recursion coefficients $\alpha_0, \alpha_1, \dots, \alpha_{n-1}$ and $\beta_1, \beta_2, \dots, \beta_{n-1}$
that converge to their exact values as $r \rightarrow \infty$.

### Implementation of the Stieltjes procedure
We describe the method used to compute the $w_i$'s and $x_i$'s in the $r$-node discrete approximation 

$$
\sum_{i =1}^r w_i  \ g(x_i)
$$

to $(u,v)$.

Suppose that $\{x: f(x) > 0\}$,  the support of the weight function $f$, is an interval with lower and upper endpoints $a$ and $b$, respectively. Here $-\infty \le a < b \le \infty$. The inner product of the functions $u$ and $v$ is therefore 

$$
\int_a^b g(x) \ f(x) \ dx,
$$

where $g = u  v$. To compute an approximation to this integral, we 
use the initial transformation described
on p.94 of Gautschi (2004). In other words, we transform the support interval with lower and upper endpoints $a$ and $b$, respectively, to the interval with lower and upper endpoints $-1$ and $1$, respectively, using the 
transformation

$$
\int_a^b g(x) \ f(x) \ dx	 
= \int_{-1}^1 g\big(\varphi(y)\big) \ f\big(\varphi(y)\big) \ \varphi'(y)  dy 
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ (3)
$$

where

$$
\varphi(y) = 
\begin{cases}
(1/2) (b - a) y + (1/2) (b + a) &\text{if } -\infty < a < b < \infty
\\
b - (1 - y) / (1 + y) &\text{if } -\infty = a < b < \infty
\\
a + (1 + y) / (1 - y) &\text{if } -\infty < a < b = \infty
\\
y / (1 - y^2)  &\text{if } -\infty = a < b = \infty.
\end{cases}
$$

It follows from this definition of the function $\varphi$ that
$$
\varphi'(y) = 
\begin{cases}
(1/2) (b - a) &\text{if } -\infty < a < b < \infty
\\
2 / (1 + y)^2 &\text{if } -\infty = a < b < \infty
\\
2 / (1 - y)^2 &\text{if } -\infty < a < b = \infty
\\
(1 + y^2) / \big(1 - y^2 \big)^2  &\text{if } -\infty = a < b = \infty.
\end{cases}
$$


A discrete approximation to the right-hand side of 	(3) can be found as follows.
Let 
$h(y) = g\big(\varphi(y)\big) \, f\big(\varphi(y)\big) \, \varphi'(y)$.
A high-quality quadrature rule is then used to provide a discrete approximation to the integral $\int_{-1}^1 h(y) \ dy$.
Gautschi (2004) uses a Fejer quadrature rule.
Instead, we use Gauss Legendre quadrature with $r$ nodes to approximate

$$
\int_{-1}^1 h(y) \ dy \quad \text{by} \quad \sum_{i =1}^r \ \xi_i \ h(y_i),
$$

where $y_1, \dots, y_r$ are the nodes and $\xi_1, \dots, \xi_r$ are the corresponding weights. Now 

$$
 \sum_{i =1}^r \ \xi_i \ h(y_i) = \sum_{i =1}^r \ w_i \ g(x_i),
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  (4)
$$

where $w_i =  \xi_i \ \varphi'(y_i) \ f\big(\varphi(y_i)\big)$
and $x_i = \varphi(y_i)$. To summarize, we 
approximate the inner product $(u, v)$ by the right-hand side of
(4).

Consider the case that either $-\infty = a$ or $b = \infty$, or both.
It was found that the computation of the $w_i$'s using 
`Double64` arithmetic, from the package `DoubleFloats`,
could result in some of these being `NaN`'s. To solve this problem, we compute the $\log(w_i)$'s instead using

$$
\log(w_i) = \log(\xi_i) + \log \big(\varphi'(y_i) \big) + \log\big(f\big(\varphi(y_i)\big) \big),
$$

where

$$
\log\big(\varphi'(y)\big) = 
\begin{cases}
\log(2) - 2 \log(1 + y) &\text{if } -\infty = a < b < \infty
\\
\log(2) - 2 \log(1 - y) &\text{if } -\infty < a < b = \infty
\\
\log(1 + y^2) - 2 \log \big(1 - y^2 \big)  &\text{if } -\infty = a < b = \infty.
\end{cases}
$$

The user needs to provide a Julia
function to evaluate $\log(f)$. Since $g = uv$, 

$$
w_i \ g(x_i)
= w_i \ u(x_i) \  v(x_i)
= \text{sign}\big(u(x_i)\big) \ \text{sign}\big(v(x_i)\big) \
w_i  \ \big|u(x_i)\big|  \ \big|v(x_i)\big|
$$

$$
= \text{sign}\big(u(x_i)\big) \ \text{sign}\big(v(x_i)\big) \
\exp\Big(\log(w_i) + \log\big(\big|u(x_i)\big|\big) + \log\big(\big|v(x_i)\big|\big) \Big),
$$

which is used to compute the right-hand side of (4).

### **Remark**
In the comments for his subroutine `qgp` (available at cs.purdue.edu/archives/2002/wxg/codes), Gautschi states of this method that
"It takes no account of the special nature of the weight function involved and hence may result in slow convergence of the discretization procedure."

For the "scaled chi pdf" weight function considered in **Example 1**, the graph of the weight function has a single peak which becomes narrower as the positive integer parameter $m$ increases. For a given number $n$ of Gauss quadrature nodes, this implies that the value of $r$, the number of nodes in the discrete approximation (4), required for 
sufficiently accurate results from the Stieltjes procedure
increases rapidly with increasing $m$.  

## **Two Examples**

For any of the following which_f's \
which_f = ["scaled chi pdf", [0,Inf], m]\
which_f = ["chemistry example", [0, Inf]]\
which_f = ["Hermite", [-Inf, Inf]]\
which_f = ["Generalized Laguerre", [0, Inf], α_GGL]\
which_f = ["Legendre", [-1, 1]],\
we use stieltjes_lnf_stored_scr.jl, which computes the $\log(f(x))$, as illustrated by the 
following example.

### **Example 4**

Consider the same weight function as in **Example 1**. For this weight function, with positive integer parameter `m` set to 160 and number of Gauss quadrature nodes `n` set to 33, we use the following
to carry out both **Step 1** and **Step 2**

<!--
The evaluation of the log-weight function $\log(f(x))$ **has been implemented** in the function `ln_scaled_chi_pdf_fn` which has inputs `T` (Floating Point type), `x` and `m` (positive integer parameter). 
-->
```julia
	m = 160
	which_f = ["scaled chi pdf", [0,Inf], m]
	n = 33
	pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
	include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))
```

This results in stieltjes_nodes and stieltjes_weights, which are
the custom Gauss quadrature nodes and weights, respectively, obtained using the Stieltjes procedure and the number of sampled function values `r` used in the discrete approximation to the 
inner product of functions used in this procedure. If `r` is of interest to us we use
```julia
	println("r = ", r)
```

<!--
This weight function is a probability density function. Therefore $\mu_0 = 1$. There is also a formula for $\mu_1$. We first compute $\mu_0$ and $\mu_1$ using

	T = BigFloat;
	which_f = ["scaled chi pdf", [0,Inf], m];
	μ₀, μ₁ = μ_offsetvec_fn(T, which_f, 1);

Carry out **Step 1**, 
the computation of the recursion coefficients in the three-step recurrence 
relation,
and print out the the chosen value of $r$, the number of nodes in the discrete approximation (4), using

	n = 5;
	a = convert(T,0);
	b = Inf;
	stieltjes_a_vec, stieltjes_b_vec, stieltjes_nbits, r = 
	stieltjes_a_vec_b_vec_final_fn(n, μ₀, lnf_fn, a, b);
	println("r = ", r)


Carry out **Step 2**, 
the computation of the Gauss quadrature nodes and weights from the recursion coefficients, using

	stieltjes_nodes, stieltjes_weights = 
	stieltjes_custom_gauss_quad_all_fn(n, μ₀, μ₁, stieltjes_a_vec, stieltjes_b_vec, a, b);
-->
Convert stieltjes_nodes and stieltjes_weights to `Float64` vectors using
```julia
	stieltjes_nodes = convert(Vector{Float64}, stieltjes_nodes)
	stieltjes_weights = convert(Vector{Float64}, stieltjes_weights)
```

Print these out in the same format as in **Example 1** using
```julia
	using Printf
	println("           stieltjes_nodes             stieltjes_weights")
	for i in 1:lastindex(stieltjes_nodes)
    	@printf "%2d     " i
    	@printf "%.16e     " stieltjes_nodes[i]
    @printf "%.16e  \n" stieltjes_weights[i]
	end
```
### **Example 5**

Consider the same weight function as in Example 3. We identify this new weight function by first assigning
a value of the positive parameter $k$ and then using Julia command
`which_f = ["weibull pdf", [0, Inf], k]`.
Therefore the following command specifies this weight function with $k$ set to 2.0, say.
```julia
which_f = ["weibull pdf", [0, Inf], 2.0]
```

In this case, 
$$
\log(f(x)) =
\log(k) + (k - 1) \log(x) -x^{k} \ \ \ \ \text{for} \ \ x > 0.
$$

We provide the following function for computing $\log(f(x))$ 
```julia
using SpecialFunctions
function lnf_weibull_pdf_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    @assert x > convert(T,0)
    k = which_f[3]
    T_k = parse(T, string(k))
    @assert T_k > convert(T, 0)
    log(T_k) + (T_k - convert(T,1)) * log(x) - x^T_k
end
```


For this weight function, with parameter $k$ set to 2.0 and number of Gauss quadrature nodes n set to 10, we use the following commands to specify the function that will be used for computing $\log(f)$, together with the value of $\mu_0$ which is 1 since the weight function is a pdf. In the following code,  $\mu_0$  is denoted by `mu0`.

```julia
	n = 10
	lnf_fn = x -> lnf_weibull_pdf_fn(T, which_f, x)
	mu0 = convert(Double64, 1)
```

Then we use the following command to carry out both **Step 1** and **Step 2**

<!-- According to Claude Opus 4.6 on 5 March 2026:
pathof(CustomGaussQuadrature) always returns the absolute path to CustomGaussQuadrature.jl, regardless of whether the package was:

activated locally (] activate .) — returns e.g. c:\Users\pkaba\...\CustomGaussQuadrature\src\CustomGaussQuadrature.jl

installed from the registry (Pkg.add(...)) — returns e.g. C:\Users\pkaba\.julia\packages\CustomGaussQuadrature\Ab1Cd\src\CustomGaussQuadrature.jl

In both cases, dirname(dirname(pathof(...))) gets you to the package root, and the include resolves to the correct file.

So this approach is universal — it works for both the local and registry-installed versions. 
-->
```julia
	pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
	include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))
```
This results in stieltjes_nodes and stieltjes_weights, which are
the custom Gauss quadrature nodes and weights, respectively, obtained using the Stieltjes procedure and the number of sampled function values `r` used in the discrete approximation to the 
inner product of functions used in this procedure. If `r` is of interest to us we use
```julia
	println("r = ", r)
```

The following commands convert the 
`Float64` vectors and print them using the `@printf` macro from the package `Printf`.
```julia
    using Printf
    nodes = convert(Vector{Float64}, stieltjes_nodes)
    weights = convert(Vector{Float64}, stieltjes_weights)
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.16e     " nodes[i]
        @printf "%.16e  \n" weights[i]
    end
```
# Step 2 using the eigenvalues and eigenvectors of the Jacobi matrix

Define the $n \times n$ Jacobi matrix 

$$
J_n = \left[
\begin{array}{cccccc}
\alpha_0 & \sqrt{\beta_1} & 0 & 0 & \dots & 0 \\
\sqrt{\beta_1} & \alpha_1 & \sqrt{\beta_2} & 0 & \ddots & \vdots \\
0 & \sqrt{\beta_2} & \alpha_2 & \sqrt{\beta_3} & 0  & 0\\
0 & \ddots & \ddots & \ddots & \ddots & 0\\
\vdots & \ddots & & \sqrt{\beta_{n-1}} & \alpha_{n-2} & \sqrt{\beta_{n-1}} \\
0 & \dots & 0 & 0 & \sqrt{\beta_{n-1}} & \alpha_{n-1} 
\end{array}
\right].
$$

The nodes $\tau_1, \dots, \tau_n$ are the eigenvalues of $J_n$, in increasing order. Let $\mathbf{x}_i$ denote an eigenvector corresponding to the eigenvalue $\tau_i$. The weight

$$
\lambda_i 
= \frac{\mu_0 \ (\text{first component of }\mathbf{x}_i)^2}
{(\mathbf{x}_i, \mathbf{x}_i)}.
$$

The  matrix $J_n$ is tridiagonal i.e. its nonzero elements are only on the subdiagonal, diagonal and superdiagonal. It is also symmetric.

We compute the eigenvalues and eigenvectors of a symmetric tridiagonal matrix using the package `GenericLinearAlgebra`.

---


### References ###

Gautschi, W. (1983). How and how not to check Gaussian quadrature formulae. BIT, 23, 209-216

Gautschi, W. (2004). Orthogonal Polynomials, 
Computation and Approximation. Oxford University Press, Oxford.

