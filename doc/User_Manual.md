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

`CustomGaussQuadrature` aims to deliver about double precision accuracy in the computed Gauss rule for subsequent extensive `Float64` computations. 

## Installation

```julia
using Pkg
Pkg.add("CustomGaussQuadrature")
```

Then load the package:

```julia
using CustomGaussQuadrature
```

In this manual and in the package source, names ending in `_fn` denote functions and filenames ending in `_scr` denote implementation scripts included in the module.

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

 We require that there is a formula, which can be computed in `BigFloat`
arithmetic, for the $s$'th moment

$$
\mu_s = \int_{-\infty}^{\infty} x^s \ f(x) \ dx
$$

for all nonnegative integers $s \le 2 n - 1$. This formula has been already implemented in code for several families of weight functions (described below). For "new" weight functions this formula must be provided by the user and implemented in code, as illustrated in  Example 3 below.


We use `BigFloat` arithmetic, with number of bits of precision `b`. For assurance that the chosen number of bits of precision `b` is sufficiently large, we use the following simple method. Suppose that $c_1$ and $c_2$ are numerical approximations to the same quantity computed using the same Julia function
that implements an ill-conditioned numerical method, 
using `BigFloat` arithmetic with number of bits of precision `b1` and `b2`, respectively, where `b2` is substantially larger than `b1`.
The approximation $c_2$ is therefore believed to be much more accurate than $c_1$. Then we take $c_2$ as our final approximation and assess 
the absolute error to be less than $|c_1 - c_2|$.

We compute approximations to the vectors $\big(\alpha_0, \dots, \alpha_{n-1} \big)$ and $\big(\sqrt{\beta_1}, \dots, \sqrt{\beta_{n-1}}\big)$ that are assessed by this method to have maximum absolute errors and maximum relative errors, respectively, bounded
above by $10^{-18}$. 

The function `custom_gauss_quad_all_fn` carries out both **Step 1** and **Step 2** to compute the Gauss quadrature nodes and weights. 
**Step 2** is usually very well conditioned and is carried out using quadruple precision arithmetic. Consequently, the computed Gauss quadrature nodes and weights are expected to have 
maximum absolute errors and maximum relative errors, respectively, bounded above by about $10^{-18}$.


We identify the weight function $f$ by the array 
`which_f` with components:\
(i) the name given to $f$ (a string),\
(ii) support interval of $f$ specified by a 2-vector of the endpoints and \
(iii) parameter vector (if any).

For user input of support interval endpoints and parameter values, use ordinary Julia values for exact quantities such as integers
and $\pm\infty$, but use strings for finite non-integer constants whose decimal
representation must be preserved. Thus `which_f = ["scaled chi pdf", [0,Inf], m]`
is natural when `m` is an integer, whereas a Weibull shape parameter should be
entered as a string such as `"2.1"`.

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

For any of the following which_f's 

which_f = ["scaled chi pdf", [0,Inf], m]\
which_f = ["chemistry example", [0, Inf]]\
which_f = ["hermite", [-Inf, Inf]]\
which_f = ["generalized laguerre", [0, Inf], α_gl]\
which_f = ["legendre", [-1, 1]],

we proceed as illustrated by the following two examples.

<!--say which_f = ["scaled chi pdf", [0,Inf], m],
we assign the value of positive integer parameter m and then use
```julia
	which_f = ["scaled chi pdf", [0,Inf], m]
	moment_fn = moment_stored_fn
-->


### **A note on the input of real-valued parameters**

Suppose that the exact value of the parameter k is 3.1. If we assign it using

```julia
k = 3.1
```

which sets k to the `Float64` approximation to 3.1. This package uses a range of `AbstractFloat` types `T`, including `BigFloat`, for its internal computations. It is therefore important to obtain the best approximation of type `T` to the exact value 3.1. Since k is of type `Float64`, the command `convert(BigFloat, k)` gives only `Float64`-precision:

```
3.100000000000000088817841970012523233890533447265625
```

Once the user enters `3.1`, Julia has already stored the nearby `Float64`
value rather than the exact decimal constant 3.1, so that lost precision
cannot be recovered later. For this reason, the recommended user-facing rule
is simple: use ordinary Julia values for exact quantities such as integers and
`±Inf`, but use strings for finite non-integer constants whose decimal
representation matters. Thus a Weibull parameter should be entered as `"3.1"`,
not `3.1`.

### **Example 1 - Scaled chi pdf weight function** 

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
`Double64` vectors for a 33-point Gauss quadrature rule. 
```julia
	moment_fn = moment_stored_fn
    n = 33
    nodes_momentdets_ex1, weights_momentdets_ex1 = custom_gauss_quad_all_fn(moment_fn, which_f, n)
```


This is followed by conversion to 
`Float64` vectors and printing using the `@printf` macro from the package `Printf`.
```julia
    using Printf
    nodes_momentdets_ex1 = convert(Vector{Float64}, nodes_momentdets_ex1)
    weights_momentdets_ex1 = convert(Vector{Float64}, weights_momentdets_ex1)
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.16e     " nodes_momentdets_ex1[i]
        @printf "%.16e  \n" weights_momentdets_ex1[i]
    end
```

### **Example 2 - Chemistry example weight function**
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
    nodes_momentdets_ex2, weights_momentdets_ex2 = custom_gauss_quad_all_fn(moment_fn, which_f, n)
```

This is followed by printing these nodes and weights in the same format as Table 2.2 of 
Gautschi (1983), using the `@printf` macro from the package `Printf`.

```julia
    using Printf
    nodes_momentdets_ex2 = convert(Vector{BigFloat}, nodes_momentdets_ex2)
    weights_momentdets_ex2 = convert(Vector{BigFloat}, weights_momentdets_ex2)
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.15e     " nodes_momentdets_ex2[i]
        @printf "%.15e  \n" weights_momentdets_ex2[i]
    end
```

These results agree with Table 2.2 of 
Gautschi (1983).


### **Example 3 - Weibull pdf (scale parameter = 1) weight function** 

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
Because finite non-integer parameters should be entered as strings, the following
command specifies this weight function with $k$ set to 3.1.
```julia
which_f = ["weibull pdf", [0, Inf], "3.1"]
```
We provide the function for computing the s'th moment using this `which_f`.

The helper function `materialize_scalar_spec_fn(T, value)` converts a user-supplied
scalar specification to a concrete value of type `T`. For example, it converts a
string such as `"3.1"` to a value of type `T`. This is useful when a parameter is
stored in `which_f` as a string so that its decimal value is preserved until the
working floating-point type `T` has been chosen.

```julia
using SpecialFunctions

function moment_weibull_pdf_fn(::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
        @assert which_f[1] == "weibull pdf"
    T_k = materialize_scalar_spec_fn(T, which_f[3])
        @assert T_k > convert(T, 0)
        @assert s ≥ 0
        if s == 0
                return convert(T, 1)
        end
        gamma(convert(T, 1) + convert(T, s) / T_k)
end
```

The following commands compute the nodes and weights, as 
`Double64` vectors for a 6-point Gauss quadrature rule. 
```julia
	moment_fn = moment_weibull_pdf_fn
    n = 6
    nodes_momentdets_ex3, weights_momentdets_ex3 = custom_gauss_quad_all_fn(moment_fn, which_f, n)
```


This is followed by conversion to 
`Float64` vectors using the commands
```julia
    nodes_momentdets_ex3 = convert(Vector{Float64}, nodes_momentdets_ex3)
    weights_momentdets_ex3 = convert(Vector{Float64}, weights_momentdets_ex3)
```

These nodes and weights are printed using the `@printf` macro from the package `Printf` as follows.
```julia
    using Printf
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.16e     " nodes_momentdets_ex3[i]
        @printf "%.16e  \n" weights_momentdets_ex3[i]
    end
```

### **Example 4 - Inverse gamma pdf (vector-valued parameter) weight function**

The following example shows how `which_f[3]` can itself be a vector of scalar
specifications. Consider the inverse gamma pdf

$$
f(x) =
\begin{cases}
\dfrac{\beta^{\alpha}}{\Gamma(\alpha)} x^{-(\alpha + 1)} \exp(-\beta / x) &\text{for } x > 0
\\
0 &\text{otherwise,}
\end{cases}
$$

where $\alpha > 0$ and $\beta > 0$. We identify this weight function by

```julia
which_f = ["inverse gamma pdf", [0, Inf], ["18.5", "3.2"]]
```

Here the parameter vector `which_f[3]` has two components, corresponding to
`α = "18.5"` and `β = "3.2"`.

For this weight function, the $s$'th moment is

$$
\beta^s \frac{\Gamma(\alpha - s)}{\Gamma(\alpha)},
$$

provided that $\alpha > s$. Since `custom_gauss_quad_all_fn` requires moments
up to order $2n - 1$, we must choose $n$ so that $\alpha > 2n - 1$. For
`n = 8` (as in the example below) this condition is satisfied.

```julia
using SpecialFunctions

function moment_inverse_gamma_pdf_fn(::Type{T}, which_f, s::Integer) where {T<:AbstractFloat}
    @assert which_f[1] == "inverse gamma pdf"
    @assert length(which_f[3]) == 2
    alpha_spec, beta_spec = which_f[3]
    T_alpha = materialize_scalar_spec_fn(T, alpha_spec)
    T_beta = materialize_scalar_spec_fn(T, beta_spec)
    @assert T_alpha > convert(T, 0)
    @assert T_beta > convert(T, 0)
    @assert s ≥ 0
    if s == 0
        return convert(T, 1)
    end
    T_s = convert(T, s)
    @assert T_alpha > T_s
    T_beta^T_s * gamma(T_alpha - T_s) / gamma(T_alpha)
end

n = 8
nodes_momentdets_ex4, weights_momentdets_ex4 = custom_gauss_quad_all_fn(moment_inverse_gamma_pdf_fn, which_f, n)
nodes_momentdets_ex4 = convert(Vector{Float64}, nodes_momentdets_ex4)
weights_momentdets_ex4 = convert(Vector{Float64}, weights_momentdets_ex4)
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
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ (3)
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
$h(y) = g\big(\varphi(y)\big) \thinspace f\big(\varphi(y)\big) \thinspace \varphi'(y)$.
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


We compute approximations to the vectors 
$\big(\alpha_0, \alpha_1, \dots, \alpha_{n-1}\big)$ and 
$\big(\sqrt{\beta_1}, \sqrt{\beta_2}, \dots, \sqrt{\beta_{n-1}} \big)$ 
that are assessed by the simple method described in the section 
**Computation of the recursion coefficients in the three-term recurrence relation using moment determinants** to have
maximum absolute errors and maximum relative errors, respectively, bounded above by $10^{-15}$. **Step 2** is usually very well conditioned and is carried out using quadruple precision arithmetic. Consequently, the computed Gauss quadrature nodes and weights are expected to have 
maximum absolute errors and maximum relative errors, respectively, bounded above by about $10^{-15}$.

### **Remark**
In the comments for his subroutine `qgp` (available at cs.purdue.edu/archives/2002/wxg/codes), Gautschi states of this method that
"It takes no account of the special nature of the weight function involved and hence may result in slow convergence of the discretization procedure."

For the "scaled chi pdf" weight function considered in **Example 1**, the graph of the weight function has a single peak which becomes narrower as the positive integer parameter $m$ increases. For a given number $n$ of Gauss quadrature nodes, this implies that the value of $r$, the number of nodes in the discrete approximation (4), required for 
sufficiently accurate results from the Stieltjes procedure
increases rapidly with increasing $m$.  

## **Examples**

For any of the following which_f's \
which_f = ["scaled chi pdf", [0,Inf], m]\
which_f = ["chemistry example", [0, Inf]]\
which_f = ["hermite", [-Inf, Inf]]\
which_f = ["generalized laguerre", [0, Inf], α_gl]\
which_f = ["legendre", [-1, 1]],\
we use stieltjes_lnf_stored_scr.jl, which computes the $\log(f(x))$, as illustrated by the 
following example.

### **Example 5 - Scaled chi pdf Stieltjes counterpart**

Consider the same weight function as in **Example 1**.
For this stored weight function, with positive integer parameter `m` set to 160
and number of Gauss quadrature nodes `n` set to 33,
we use the high-level stored Stieltjes driver
`stieltjes_lnf_stored_scr.jl`
to carry out both **Step 1** and **Step 2**.

```julia
m = 160
which_f = ["scaled chi pdf", [0, Inf], m]
n = 33

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))
```

This creates the following outputs in the caller:

- `nodes_stieltjes`
- `weights_stieltjes`
- `a_vec_stieltjes`
- `b_vec_stieltjes`
- `nbits_stieltjes`
- `r`

Here `nodes_stieltjes` and `weights_stieltjes`
are the custom Gauss quadrature nodes and weights obtained using the Stieltjes procedure,
and `r` is the number of sampled function values used in the discrete approximation to the inner product.
The older `stieltjes_*` names remain available as compatibility aliases.
If `r` is of interest, use

```julia
println("r = ", r)
```

Convert `nodes_stieltjes` and `weights_stieltjes` to `Float64` vectors using

```julia
nodes_stieltjes = convert(Vector{Float64}, nodes_stieltjes)
weights_stieltjes = convert(Vector{Float64}, weights_stieltjes)
```

Print these out in the same format as in **Example 1** using

```julia
using Printf
println("              nodes                     weights")
for i in 1:lastindex(nodes_stieltjes)
    @printf "%2d     " i
    @printf "%.16e     " nodes_stieltjes[i]
    @printf "%.16e  \n" weights_stieltjes[i]
end
```

If we also wish to compare the Stieltjes result with the result obtained by the
moment-determinants method, we refer to `nodes_momentdets_ex1` and
`weights_momentdets_ex1` computed earlier in **Example 1**.

We may then compare the two methods using

```julia
diff_nodes = nodes_stieltjes - nodes_momentdets_ex1
rel_diff_weights = (weights_stieltjes - weights_momentdets_ex1) ./ weights_momentdets_ex1

println("maximum(abs.(nodes_stieltjes - nodes_momentdets_ex1)) = ",
    convert(Float64, maximum(abs.(diff_nodes))))
println("maximum(abs.((weights_stieltjes - weights_momentdets_ex1) ./ weights_momentdets_ex1)) = ",
    convert(Float64, maximum(abs.(rel_diff_weights))))
```

### **Example 6 - Chemistry example Stieltjes counterpart**

Consider the same weight function as in **Example 2**.
For this stored weight function, with number of Gauss quadrature nodes
`n = 15`, we again use the high-level stored Stieltjes driver
`stieltjes_lnf_stored_scr.jl`.

```julia
which_f = ["chemistry example", [0, Inf]]
n = 15

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
include(joinpath(pkg_dir, "src", "stieltjes_lnf_stored_scr.jl"))
```

This creates the following outputs in the caller:

- `nodes_stieltjes`
- `weights_stieltjes`
- `a_vec_stieltjes`
- `b_vec_stieltjes`
- `nbits_stieltjes`
- `r`

Convert `nodes_stieltjes` and `weights_stieltjes` to `Float64` vectors using

```julia
nodes_stieltjes = convert(Vector{Float64}, nodes_stieltjes)
weights_stieltjes = convert(Vector{Float64}, weights_stieltjes)
```

Print these out in the same format as in **Example 2** using

```julia
using Printf
println("              nodes                     weights")
for i in 1:lastindex(nodes_stieltjes)
    @printf "%2d     " i
    @printf "%.15e     " nodes_stieltjes[i]
    @printf "%.15e  \n" weights_stieltjes[i]
end
```

If we also wish to compare the Stieltjes result with the result obtained by the
moment-determinants method, we refer to `nodes_momentdets_ex2` and
`weights_momentdets_ex2` computed earlier in **Example 2**.

We may then compare the two methods using

```julia
diff_nodes = nodes_stieltjes - nodes_momentdets_ex2
rel_diff_weights = (weights_stieltjes - weights_momentdets_ex2) ./ weights_momentdets_ex2

println("maximum(abs.(nodes_stieltjes - nodes_momentdets_ex2)) = ",
    convert(Float64, maximum(abs.(diff_nodes))))
println("maximum(abs.((weights_stieltjes - weights_momentdets_ex2) ./ weights_momentdets_ex2)) = ",
    convert(Float64, maximum(abs.(rel_diff_weights))))
```

### **Example 7 - Weibull pdf Stieltjes counterpart**

Consider the same weight function as in **Example 3**.
We identify this new weight function by first assigning
a value of the positive parameter `k`
and then using the Julia command
`which_f = ["weibull pdf", [0, Inf], k]`.
Therefore the following command specifies this weight function with `k` set to 3.1, say.

```julia
which_f = ["weibull pdf", [0, Inf], "3.1"]
```

In this case,

$$
\log(f(x)) =
\log(k) + (k - 1) \log(x) - x^k
\qquad \text{for } x > 0.
$$

We provide the following function for computing $\log(f(x))$:

Here `materialize_scalar_spec_fn(T, which_f[3])` converts the parameter stored in
`which_f[3]` to the floating-point type `T`.

```julia
using SpecialFunctions

function lnf_weibull_pdf_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
    @assert which_f[1] == "weibull pdf"
    @assert x > convert(T, 0)
    T_k = materialize_scalar_spec_fn(T, which_f[3])
    @assert T_k > convert(T, 0)
    log(T_k) + (T_k - convert(T, 1)) * log(x) - x^T_k
end
```

For a user-defined weight function,
the high-level Stieltjes driver
`stieltjes_lnf_new_scr.jl`
expects the following inputs in the caller:

- `which_f`
- `n`
- `lnf_typed_fn`
- `mu0`

For the Weibull pdf with shape parameter `k = "3.1"`
and number of Gauss quadrature nodes `n = 6`,
we therefore use

```julia
which_f = ["weibull pdf", [0, Inf], "3.1"]
n = 6

lnf_typed_fn = lnf_weibull_pdf_fn
mu0 = 1

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))
```

Later in this example, after defining `moment_weibull_pdf_fn`, we compare
`nodes_stieltjes` and `weights_stieltjes` with the corresponding Gauss rule
obtained using the moment-determinants method.

This creates the following outputs in the caller:

- `nodes_stieltjes`
- `weights_stieltjes`
- `a_vec_stieltjes`
- `b_vec_stieltjes`
- `nbits_stieltjes`
- `r`

Here `nodes_stieltjes` and `weights_stieltjes`
are the custom Gauss quadrature nodes and weights obtained using the Stieltjes procedure,
and `r` is the number of sampled function values used in the discrete approximation to the inner product.
The older `stieltjes_*` names remain available as compatibility aliases.
If `r` is of interest, use

```julia
println("r = ", r)
```

The following commands convert the nodes and weights to `Float64` vectors and print them using the `@printf` macro from `Printf`.

```julia
using Printf
nodes = convert(Vector{Float64}, nodes_stieltjes)
weights = convert(Vector{Float64}, weights_stieltjes)
println("              nodes                     weights")
for i in 1:n
    @printf "%2d     " i
    @printf "%.16e     " nodes[i]
    @printf "%.16e  \n" weights[i]
end
```

If we also wish to compare the Stieltjes result with the result obtained by the
moment-determinants method, we refer to `nodes_momentdets_ex3` and
`weights_momentdets_ex3` computed earlier in **Example 3**.

We may then compare the two methods using

```julia
diff_nodes = nodes_stieltjes - nodes_momentdets_ex3
rel_diff_weights = (weights_stieltjes - weights_momentdets_ex3) ./ weights_momentdets_ex3

println("maximum(abs.(nodes_stieltjes - nodes_momentdets_ex3)) = ",
    convert(Float64, maximum(abs.(diff_nodes))))
println("maximum(abs.((weights_stieltjes - weights_momentdets_ex3) ./ weights_momentdets_ex3)) = ",
    convert(Float64, maximum(abs.(rel_diff_weights))))
```

### **Example 8 - Inverse gamma pdf Stieltjes counterpart**

Consider the same weight function as in **Example 4**.
The same inverse gamma pdf can also be handled by the Stieltjes procedure.
Its log-weight function is

$$
\log(f(x)) =
\alpha \log(\beta) - \log \Gamma(\alpha) - (\alpha + 1)\log(x) - \beta / x
\qquad \text{for } x > 0.
$$

We again use

```julia
which_f = ["inverse gamma pdf", [0, Inf], ["18.5", "3.2"]]
```

and since this is a probability density function, `mu0 = 1`.

```julia
using SpecialFunctions

function lnf_inverse_gamma_pdf_fn(::Type{T}, which_f, x::AbstractFloat) where {T<:AbstractFloat}
    @assert which_f[1] == "inverse gamma pdf"
    @assert x > zero(T)
    @assert length(which_f[3]) == 2
    alpha_spec, beta_spec = which_f[3]
    T_alpha = materialize_scalar_spec_fn(T, alpha_spec)
    T_beta = materialize_scalar_spec_fn(T, beta_spec)
    @assert T_alpha > zero(T_alpha)
    @assert T_beta > zero(T_beta)
    T_one = one(T)
    T_alpha * log(T_beta) -
        logabsgamma(T_alpha)[1] -
        (T_alpha + T_one) * log(x) -
        T_beta / x
end

n = 8
lnf_typed_fn = lnf_inverse_gamma_pdf_fn
mu0 = 1

pkg_dir = dirname(dirname(pathof(CustomGaussQuadrature)))
include(joinpath(pkg_dir, "src", "stieltjes_lnf_new_scr.jl"))

nodes_stieltjes = convert(Vector{Float64}, nodes_stieltjes)
weights_stieltjes = convert(Vector{Float64}, weights_stieltjes)
println("r = ", r)
```

If we also wish to compare the Stieltjes result with the result obtained by the
moment-determinants method, we refer to `nodes_momentdets_ex4` and
`weights_momentdets_ex4` computed earlier in the inverse gamma moment-determinants
example.

```julia
diff_nodes = nodes_stieltjes - nodes_momentdets_ex4
rel_diff_weights = (weights_stieltjes - weights_momentdets_ex4) ./ weights_momentdets_ex4

println("maximum(abs.(nodes_stieltjes - nodes_momentdets_ex4)) = ",
    convert(Float64, maximum(abs.(diff_nodes))))
println("maximum(abs.((weights_stieltjes - weights_momentdets_ex4) ./ weights_momentdets_ex4)) = ",
    convert(Float64, maximum(abs.(rel_diff_weights))))
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

