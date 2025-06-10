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
readily accessible via several Julia packages. We call this the 
classical case. In the non-classical case, however, the Gauss rule needs to be custom-made.

Gauss rules are usually
computed in the following two steps.

**Step 1**: Compute the recursion coefficients in the three-step recurrence 
relation (Gautschi, 2004).

**Step 2**: Use these recursion coefficients to compute the Gauss quadrature nodes 
and weights.

In the classical case there are  simple formulae for the recursion coefficients, making **Step 1** trivial. In the non-classical case Step 1 is difficult. The treatise of Gautschi (2004) describes several methods for this computation that differ widely in both complexity and sensitivity to roundoff errors. However, **Step 2** remains the same for both classical and non-classical cases. We carry out **Step 2** by computing the eigenvalues and eigenvectors of a symmetric tridiagonal matrix using the package `GenericLinearAlgebra.jl`.

# Step 1 using either moment determinants or the Stjieltjes procedure

## **Three-term recurrence relation** 


The norm of the function $u$ is denoted by $\lVert u \rVert$
and is defined to be the square root of 

$$\int_{-\infty}^{\infty} u^2(x) f(x) dx,$$

provided that this integral exists. For two functions $u$ and $v$ that satisfy $\lVert u \rVert <\infty$ and $\lVert v \rVert <\infty$, the inner product of these two functions is denoted $(u,v)$ and is defined to be

$$\int_{-\infty}^{\infty} u(x) v(x) f(x) dx.$$

If $(u,v) = 0$ then the functions $u$ and $v$ are said to be orthogonal. 

We consider only polynomials whose coefficients are real numbers. 
A polynomial in $x$ of degree $n$ is said to be monic if the coefficient of $x^n$ is 1. Let $\pi_k$ denote a monic polynomial of degree $k$. The monic polynomials $\pi_0, \pi_1, \pi_2, \dots$ are called monic orthogonal polynomials with respect to the weight function $f$ if 

$$(\pi_k, \pi_{\ell}) = 0 \ \ \text{for all} \ k \ne \ell, \ \text{where} \ k, \ell \in \{0, 1, 2, \dots\}$$

and 

$$\lVert \pi_k \rVert > 0 \ \ \text{for} \ k=0, 1, 2, \dots.$$

The Gauss quadrature nodes $\tau_1, \dots, \tau_n$ are the $n$ distinct roots of the polynomial $\pi_n$. 

The monic orthogonal polynomials with respect to the weight function $f$ satisfy the following three-term recurrence relation, see e.g.
Theorem 1.27 on p.10 of Gautschi (2004). Let 
$\pi_{-1} \equiv 0$ and $\pi_0 \equiv 1$. Then 

$$\pi_{k+1}(x) = 	(x - \alpha_k)  \pi_k(x) - \beta_k  \pi_{k-1}(x) \ \ \mathrm{for} \ k = 0, 1, 2, \ldots, $$ 

where 

$$\alpha_k 
	= \frac{\big(x \, \pi_k, \pi_k\big)}{\big(\pi_k, \pi_k\big)}
	(k = 0, 1, ...) \ \text{and} \
	\beta_k 
	= \frac{\big(\pi_k, \pi_k\big)}{\big(\pi_{k-1}, \pi_{k-1}\big)}
	(k = 1, 2, ...)
    $$


## **Computation of the recursion coefficients in the three-term recurrence relation using moment determinants** 

We carry out **Step 1** using moment determinants (Theorem 2.2 of Gautschi, 2004).
The map from the vector 
$\big(\mu_0, \mu_1, \dots, \mu_{2n-1} \big)$ of moments to the vector $\big(\alpha_0, \dots, \alpha_{n-1},  \beta_1, \dots, \beta_{n-1}\big)$ of recursion coefficients is severely ill-conditioned (Gautschi, 1983, 2004). However, it has long been recognized that this limitation can be overcome by the use of high-precision arithmetic, see e.g. Gautschi (1983), when applying 
 the method of moment determinants.
We use `BigFloat` arithmetic, with arbitrary number of bits of precision `b`. For assurance that the chosen number of bits of precision `b` is sufficiently large, we use the following simple method. Suppose that $c$ and $\widetilde{c}$ are numerical approximations to the same quantity, where $\widetilde{c}$ is believed to be much more accurate. Then we take $\widetilde{c}$ as our final approximation and assess 
the absolute error to be less than $|c - \widetilde{c}|$.


We require that there is a formula, which can be computed in `BigFloat`
arithmetic, for the $r$'th moment

$$
\mu_r = \int_{-\infty}^{\infty} x^r \, f(x) \, dx
$$

for all nonnegative integers $r \le 2 n - 1$. This formula 
must be inserted by the user into the code for the function `moment_fn`, as illustrated in the examples below.
The function `custom_gauss_quad_all_fn` is then used to compute the Gauss quadrature nodes and weights. 

The aim of the method implemented in `custom_gauss_quad_all_fn` is for (a) the *absolute errors* of each of the nodes to be less than
$10^{-17}$ and (b) the *relative errors* of each of the weights to be less than $10^{-17}$. In other words, the nodes and weights are computed with the purpose of being used in further extensive `Float64` computations.



The inputs to   `custom_gauss_quad_all_fn` are
`which_f`, `n`, `upto_n` and `extra_check`. We describe this function for the default values of `upto_n` and `extra_check`, so that the inputs to this function are as follows.

- `which_f` specifies the weight function $f$. It has the following 
components:  (i)  name,
(ii) support specified by a two-vector of the endpoints 
of the interval, with finite endpoints specified by strings of numbers
that are later converted to the appropriate 
floating-point type using `parse()` and 
(iii) parameter vector (if any).

- `n` is the number of nodes. 

For `custom_gauss_quad_all_fn` the outputs are `nodes` and `weights`  in the form of `Double64` vectors (from the package `DoubleFloats`).


### **Example 1** 
Consider the "scaled chi pdf" weight function is

$$
	f(x) =
	\begin{cases}
		\dfrac{m^{m/2}}{\Gamma(m/2) \, 2^{(m/2) - 1}} \; x^{m-1} \, \exp\big(- m \, x^2 /2\big) &\text{for } x > 0
		\\
		0 &\text{otherwise},
	\end{cases}	
$$

where $m$ is a positive integer parameter. 
This weight function is specified by

    which_f = ["scaled chi pdf", [0,Inf], m]

for some assigned value of $m$. The first component is the name given to $f$, the second component is the support interval of $f$ specified by a 2-vector of the endpoints and the last component is the parameter $m$ 

For this weight function, the $r$'th moment is 

$$
	\left(\frac{2}{m} \right)^{r/2}
	\frac{\Gamma\big((r+m)/2\big)}{\Gamma(m/2)}
    = \exp \left(\frac{r}{2} \log\left(\frac{2}{m} \right)  +
    \log \Gamma \left(\frac{r+m}{2}\right) 
    - \log \Gamma \left(\frac{m}{2}\right) \right)
$$

for $r = 0, 1, 2, \dots$.
This formula **has been implemented** in the function `moment_fn` which has inputs `T` (Floating Point type), `which_f` and `r` (moment order), as follows.

    if r == 0
        return(convert(T, 1))
    end
    T_2 = convert(T, 2)
    m = which_f[3]
    T_m = convert(T, m)
    T_r = convert(T, r)
    term1 = (T_r/T_2) * log(T_2 / T_m)
    term2 = (logabsgamma((T_r + T_m)/T_2))[1]
    term3 = (logabsgamma(T_m/T_2))[1]
    moment = exp(term1 + term2 - term3)

The following commands compute the nodes and weights, as 
`Double64` vectors for a 5-point Gauss quadrature rule. This is followed by conversion to 
`Float64` vectors and printing using  using the `printf` function from the package `Printf`.

    using CustomGaussQuad
    which_f = ["scaled chi pdf", [0,Inf], 160]::Vector{Any};
    n = 5;
    nodes, weights = custom_gauss_quad_all_fn(which_f, n);
    nodes = convert(Vector{Float64}, nodes);
    weights = convert(Vector{Float64}, nodes);
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.16e     " nodes[i]
        @printf "%.16e  \n" weights[i]
    end


### **Example 2**

The "chemistry example" weight function is

$$
		f(x) =
	\begin{cases}
		\exp(-x^3 / 3) &\ \text{for} \ x > 0
		\\
		0  &\ \text{otherwise.}
	\end{cases}	
$$

This weight function is considered by Gautschi (1983) and is specified by

    which_f = ["chemistry example", [0,Inf]]::Vector{Any};

The first component is the name given to $f$ and the second component is the support interval of $f$ specified by a 2-vector of the endpoints.

For this weight function, the $r$'th moment is 

$$
	3^{(r - 2) / 3} \, \Gamma\big((r + 1) / 3\big)
$$

for $r = 0, 1, 2, \dots$.
This formula **has been implemented** in the function `moment_fn` which has inputs `T` (Floating Point type), `which_f` and `r` (moment order), as follows.




    T_1 = convert(T, 1)
    T_2 = convert(T, 2)
    T_3 = convert(T, 3)
    T_r = convert(T, r)
    term1 = T_3^((T_r - T_2) / T_3)
    term2 = gamma((T_r + T_1) / T_3)
    moment = term1 * term2


The following commands compute the nodes and weights, as 
`Double64` vectors for a 15-point Gauss quadrature rule. This is followed by printing these nodes and weights in the same format as Table 2.2 of 
Gautschi (1983), using the `printf` function from the package `Printf`.

    using CustomGaussQuad
    which_f = ["chemistry example", [0, Inf]]::Vector{Any}
    n= 15;
    nodes, weights = custom_gauss_quad_all_fn(which_f, n);  
    nodes = convert(Vector{BigFloat}, nodes);
    weights = convert(Vector{BigFloat}, weights);
    @printf "              nodes                     weights"
    for i in 1:n
        @printf "%2d     " i
        @printf "%.15e     " nodes[i]
        @printf "%.15e  \n" weights[i]
    end

These results agree with Table 2.2 of 
Gautschi (1983).


## **Computation of the recursion coefficients in the three-term recurrence relation using the Stjieltjes procedure**

The Stjeltjes procedure, given in subsections 2.2.2 and 2.2.3 of Gautschi (2004), applies the three-term recurrence relation (1) and (2) iteratively through the following sequence:

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

### Implementation of the Stjieltjes procedure
We describe the method used to compute the $w_i$'s and $x_i$'s in the $r$-node discrete approximation 

$$
\sum_{i =1}^r w_i \, g(x_i)
$$

to $(u,v)$.

Suppose that $\{x: f(x) > 0\}$,  the support of the weight function $f$, is an interval with lower and upper endpoints $a$ and $b$, respectively. Here $-\infty \le a < b \le \infty$. The inner product of the functions $u$ and $v$ is therefore 
$$
	\int_a^b g(x) \, f(x) \, dx,
$$
where $g = u \, v$. To compute an approximation to this integral, we 
we use the initial transformation described
on p.94 of Gautschi (2004). In other words, we transform the support interval with lower and upper endpoints $a$ and $b$, respectively, to the interval with lower and upper endpoints $-1$ and $1$, respectively, using the 
transformation

$$
	\int_a^b g(x) \, f(x) \, dx	 
	= \int_{-1}^1 g\big(\varphi(y)\big) \, f\big(\varphi(y)\big) \, \varphi'(y) \, dy
    \tag{3}
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


A discrete approximation to the right-hand side of 	(3) can be found as follows.
Let 
$h(y) = g\big(\varphi(y)\big) \, f\big(\varphi(y)\big) \, \varphi'(y)$, where

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

A high-quality quadrature rule then provides a discrete approximation to the integral $\int_{-1}^1 h(y) \, dy$.
Gautschi (2004) uses a Fejer quadrature rule.
Instead, we use Gauss Legendre quadrature with $r$ nodes to approximate

$$
	\int_{-1}^1 h(y) \, dy \quad \text{by} \quad \sum_{i =1}^r \xi_i \, h(y_i),	
$$

where $y_1, \dots, y_r$ are the nodes and $\xi_1, \dots, \xi_r$ are the corresponding weights. Now 

$$
	\sum_{i =1}^r \xi_i \, h(y_i) = \sum_{i =1}^r w_i \, g(x_i),	
    \tag{4}
$$

where $w_i =  \xi_i \, \varphi'(y_i) \, f\big(\varphi(y_i)\big)$
and $x_i = \varphi(y_i)$. To summarize, we 
approximate the inner product $(u, v)$ by
(4).
All that the user needs to provide is a Julia
function to evaluate $f$.


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

All that the user needs to provide is a Julia
function to evaluate $\log(f)$. Since $g = uv$, 

$$
w_i \, g(x_i)
= w_i \, u(x_i) \, v(x_i)
= \text{sign}\big(u(x_i)\big) \, \text{sign}\big(v(x_i)\big) \,
w_i \, \big|u(x_i)\big| \, \big|v(x_i)\big|
$$

$$
= \text{sign}\big(u(x_i)\big) \, \text{sign}\big(v(x_i)\big) \,
\exp\Big(\log(w_i) + \log\big(\big|u(x_i)\big|\big) + \log\big(\big|v(x_i)\big|\big) \Big),
$$

which is used to compute (4).


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

