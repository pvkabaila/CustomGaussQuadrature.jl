# CustomGaussQuadrature

[![Build Status](https://github.com/pvkabaila/CustomGaussQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/pvkabaila/CustomGaussQuadrature.jl/actions/workflows/CI.yml?query=branch%3Amaster)

Suppose that we wish to evaluate
$$
\int_{-\infty}^{\infty} g(x) f(x) dx,
$$
where $f$ is a specified nonnegative integrable
weight function. The $n$-point Gauss quadrature approximation to this
integral is
$$
\sum_{i=1}^n \lambda_i \, g(\tau_i),
$$
where $\tau_1, \dots, \tau_n$ are called the nodes
and $\lambda_1, \dots, \lambda_n$ are called the
corresponding weights. The dependence on these 
nodes and weights on $n$ is implicit.
This approximation is exact
whenever $g$ is a polynomial of degree less
than or equal to $2n - 1$. 

If $f$ takes a form that leads to Gauss quadrature
rules with nodes that are the roots of classical orthogonal
polynomials of a continuous variable then these rules
(such as Gauss Legendre, Gauss Hermite and Gauss Laguerre) are
readily accessible via several Julia packages.
If, however, $f$ does not take one of these
particular forms then the Gauss quadrature nodes and weights
need to be custom-made.


Usually Gauss quadrature nodes and weights are 
computed in the following two steps.

$\quad$ **Step 1**: Compute the recursion coefficients in the three-step recurrence 
relation (Gautschi, 2004).

$\quad$ **Step 2**: Use these recursion coefficients to compute the Gauss quadrature nodes 
and weights.

**Step 2** remains the same for any Gauss quadrature rule computation, 
including custom-made. **Step 1** is particularly simple when $f$ takes a form that leads to 
Gauss quadrature
rules with nodes that are the roots of classical orthogonal
polynomials of a continuous variable. However, for custom-made Gauss 
quadrature this step is complicated, see e.g. Gautschi (2004).  

We carry out **Step 1** using moment determinants (Theorem 2.2 of Gautschi, 2004). Although this method
is severely ill-conditioned, this 
limitation is overcome by the use of `BigFloat` arithmetic. 
We require that there is a formula, which can be computed in `BigFloat`
arithmetic, for the $r$'th moment
$$
\int_{-\infty}^{\infty} x^r \, f(x) \, dx
$$
for all nonnegative integers $r \le 2 n - 1$. This formula 
must be inserted by the user into the code for the function `moment_fn`, as illustrated in the examples below.
The function `custom_gauss_quad_all_fn` is then used to compute the Gauss quadrature nodes and weights. 

The aim of the method implemented in `custom_gauss_quad_all_fn` is for (a) the *absolute errors* of each of the nodes to be less than
$10^{-17}$ and (b) the *relative errors* of each of the weights to be less than $10^{-17}$. In other words, the nodes and weights are computed with the purpose of being used in further `Float64` computations.

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


## **Example 1** 
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

    which_f = ["scaled chi pdf", ["0.0",Inf], m]

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
    term2 = lgamma((T_r + T_m)/T_2)
    term3 = lgamma(T_m/T_2)
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


## **Example 2**

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

    which_f = ["chemistry example", ["0.0",Inf]]::Vector{Any};

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
    println("              nodes                     weights")
    for i in 1:n
        @printf "%2d     " i
        @printf "%.15e     " nodes[i]
        @printf "%.15e  \n" weights[i]
    end

These results agree with Table 2.2 of 
Gautschi (1983).

---


### References ###

Gautschi, W. (1983). How and how not to check Gaussian quadrature formulae. BIT, 23, 209-216

Gautschi, W. (2004). Orthogonal Polynomials, 
Computation and Approximation. Oxford University Press, Oxford.
Copied from earlier package
