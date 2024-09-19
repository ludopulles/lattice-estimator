# -*- coding: utf-8 -*-

from copy import copy
from dataclasses import dataclass

from sage.all import binomial, ceil, exp, log, oo, parent, pi, RealField, RR, sqrt


def stddevf(sigma):
    """
    Gaussian width parameter σ → standard deviation.

    :param sigma: Gaussian width parameter σ

    EXAMPLE::

        >>> from estimator import *
        >>> ND.stddevf(64.0)
        25.532...

        >>> ND.stddevf(64)
        25.532...

        >>> ND.stddevf(RealField(256)(64)).prec()
        256
    """
    try:
        prec = parent(sigma).prec()
    except AttributeError:
        prec = 0
    if prec > 0:
        FF = parent(sigma)
    else:
        FF = RR
    return FF(sigma) / FF(sqrt(2 * pi))


def sigmaf(stddev):
    """
    Standard deviation → Gaussian width parameter σ.

    :param stddev: standard deviation

    EXAMPLE::

        >>> from estimator import *
        >>> n = 64.0
        >>> ND.sigmaf(ND.stddevf(n))
        64.000...

        >>> ND.sigmaf(RealField(128)(1.0))
        2.5066282746310005024157652848110452530
        >>> ND.sigmaf(1.0)
        2.506628274631...
        >>> ND.sigmaf(1)
        2.506628274631...
    """
    RR = parent(stddev)
    #  check that we got ourselves a real number type
    try:
        if abs(RR(0.5) - 0.5) > 0.001:
            RR = RealField(53)  # hardcode something
    except TypeError:
        RR = RealField(53)  # hardcode something
    return RR(sqrt(2 * pi)) * stddev


@dataclass
class NoiseDistribution:
    """
    All noise distributions are instances of this class.
    """
    stddev: float = 0
    mean: float = 0
    n: int = None
    bounds: tuple = (-oo, oo)
    density: float = 1.0  # hamming_weight() / n.

    def __lt__(self, other):
        """
        We compare distributions by comparing their standard deviation.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.DiscreteGaussian(2.0) < ND.CenteredBinomial(18)
            True
            >>> ND.DiscreteGaussian(3.0) < ND.CenteredBinomial(18)
            False
            >>> ND.DiscreteGaussian(4.0) < ND.CenteredBinomial(18)
            False

        """
        try:
            return self.stddev < other.stddev
        except AttributeError:
            return self.stddev < other

    def __le__(self, other):
        """
        We compare distributions by comparing their standard deviation.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.DiscreteGaussian(2.0) <= ND.CenteredBinomial(18)
            True
            >>> ND.DiscreteGaussian(3.0) <= ND.CenteredBinomial(18)
            True
            >>> ND.DiscreteGaussian(4.0) <= ND.CenteredBinomial(18)
            False

        """
        try:
            return self.stddev <= other.stddev
        except AttributeError:
            return self.stddev <= other

    def __str__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> ND.DiscreteGaussianAlpha(0.01, 7681)
            D(σ=30.64)

        """
        if self.n:
            return f"D(σ={float(self.stddev):.2f}, μ={float(self.mean):.2f}, n={int(self.n)})"
        else:
            return f"D(σ={float(self.stddev):.2f}, μ={float(self.mean):.2f})"

    def __repr__(self):
        if self.mean == 0.0:
            return f"D(σ={float(self.stddev):.2f})"
        else:
            return f"D(σ={float(self.stddev):.2f}, μ={float(self.mean):.2f})"

    def __hash__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> hash(ND(3.0, 1.0)) == hash((3.0, 1.0, None))
            True

        """
        return hash((self.stddev, self.mean, self.n))

    def __len__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> D = ND.SparseTernary(1024, p=128, m=128)
            >>> len(D)
            1024
            >>> int(round(len(D) * float(D.density)))
            256

        """
        if self.n is None:
            raise ValueError("Distribution has no length.")
        return self.n

    def resize(self, new_n):
        """
        Return an altered distribution having a dimension `new_n`.
        """
        new_self = copy(self)
        new_self.n = new_n
        return new_self

    @property
    def hamming_weight(self):
        """
        The number of non-zero coefficients in this distribution
        """
        return round(len(self) * float(self.density))

    @property
    def is_bounded(self):
        """
        Whether the value of coefficients are bounded
        """
        return (self.bounds[1] - self.bounds[0]) < oo

    @property
    def is_Gaussian_like(self):
        return False

    @property
    def is_sparse(self):
        """
        Whether the density of the distribution is < 1/2.
        Note: 1/2 might be considered somewhat arbitrary.
        """
        # NOTE: somewhat arbitrary
        return self.density < 0.5

    def support_size(self, fraction=1.0):
        raise NotImplementedError("support_size")


"""
The follow noise distributions are implemented below:
- DiscreteGaussian
- DiscreteGaussianAlpha
- CenteredBinomial
- Uniform
- UniformMod
- SparseTernary
"""


class DiscreteGaussian(NoiseDistribution):
    """
    A discrete Gaussian distribution with standard deviation ``stddev`` per component.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.DiscreteGaussian(3.0, 1.0)
        D(σ=3.00, μ=1.00)
    """
    # cut-off for Gaussian distributions
    gaussian_tail_bound: int = 2
    # probability that a coefficient falls within the cut-off
    gaussian_tail_prob: float = 1 - 2 * exp(-4 * pi)

    def __init__(self, stddev, mean=0, n=None):
        super().__init__(stddev=stddev, mean=mean, n=n)

        b_val = oo if n is None else ceil(log(n, 2) * stddev)
        self.bounds = (-b_val, b_val)
        self.density = max(0.0, 1 - RR(1 / sigmaf(stddev)))

    @property
    def is_Gaussian_like(self):
        return True

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.DiscreteGaussian(1.0, n=128).support_size(0.99)
            2.686...e+174
        """
        # We will treat this noise distribution as bounded with failure probability `1 - fraction`.
        n = len(self)
        t = self.gaussian_tail_bound
        p = self.gaussian_tail_prob

        if p**n < fraction:
            raise NotImplementedError(
                f"TODO(DiscreteGaussian.support_size): raise t. {RR(p ** n)}, {n}, {fraction}"
            )

        b = 2 * t * sigmaf(self.stddev) + 1
        return (2 * b + 1)**n


def DiscreteGaussianAlpha(alpha, q, mean=0, n=None):
    """
    A discrete Gaussian distribution with standard deviation α⋅q/√(2π) per component.

    EXAMPLE::

        >>> from estimator import *
        >>> alpha, q = 0.001, 2048
        >>> ND.DiscreteGaussianAlpha(alpha, q)
        D(σ=0.82)
        >>> ND.DiscreteGaussianAlpha(alpha, q) == ND.DiscreteGaussian(ND.stddevf(alpha * q))
        True
    """
    return DiscreteGaussian(RR(stddevf(alpha * q)), RR(mean), n)


class CenteredBinomial(NoiseDistribution):
    """
    Sample a_1, …, a_η, b_1, …, b_η uniformly from {0, 1}, and return Σ(a_i - b_i).

    EXAMPLE::

        >>> from estimator import *
        >>> ND.CenteredBinomial(8)
        D(σ=2.00)
    """
    def __init__(self, eta, n=None):
        super().__init__(
            density=1 - binomial(2 * eta, eta) * 2 ** (-2 * eta),
            stddev=RR(sqrt(eta / 2.0)),
            bounds=(-eta, eta),
            n=n,
        )

    @property
    def is_Gaussian_like(self):
        return True

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.CenteredBinomial(3, 10).support_size()
            282475249
            >>> ND.CenteredBinomial(3, 10).support_size(0.99)
            279650497
        """
        # TODO: this might be suboptimal/inaccurate for binomial distribution
        a, b = self.bounds
        return ceil(RR(fraction) * (b - a + 1)**len(self))


class Uniform(NoiseDistribution):
    """
    Uniform distribution ∈ ``[a,b]``, endpoints inclusive.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.Uniform(-3, 3)
        D(σ=2.00)
        >>> ND.Uniform(-4, 3)
        D(σ=2.29, μ=-0.50)
    """
    def __init__(self, a, b, n=None):
        if b < a:
            raise ValueError(f"upper limit must be larger than lower limit but got: {b} < {a}")

        m = b - a + 1
        super().__init__(
            stddev=RR(sqrt((m**2 - 1) / 12)),
            mean=RR((a + b) / 2),
            bounds=(a, b),
            density=(1 - 1 / m if a <= 0 and b >= 0 else 1),
            n=n,
        )

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.Uniform(-3,3, 64).support_size(0.99)
            1207562882759477428726191443614714994252339953407098880
        """
        # TODO: this might be suboptimal/inaccurate for binomial distribution
        a, b = self.bounds
        return ceil(RR(fraction) * (b - a + 1)**len(self))


def UniformMod(q, n=None):
    """
    Uniform mod ``q``, with balanced representation.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.UniformMod(7)
        D(σ=2.00)
        >>> ND.UniformMod(8)
        D(σ=2.29, μ=0.50)
        >>> ND.UniformMod(2) == ND.Uniform(0, 1)
        True
    """
    half_q = q // 2
    return Uniform(half_q - q + 1, half_q, n=n)


class SparseTernary(NoiseDistribution):
    """
    Distribution of vectors of length ``n`` with ``p`` entries of 1 and ``m`` entries of -1, rest 0.

    EXAMPLE::
        >>> from estimator import *
        >>> ND.SparseTernary(100, p=10)
        D(σ=0.45)
        >>> ND.SparseTernary(100, p=10, m=10)
        D(σ=0.45)
        >>> ND.SparseTernary(100, p=10, m=8)
        D(σ=0.42, μ=0.02)
    """
    def __init__(self, n, p, m=None):
        if m is None:
            m = p
        self.p, self.m = p, m

        # Yes, n=0 might happen when estimating the cost of the dual attack!
        mean = 0 if n == 0 else RR((p - m) / n)
        density = 0 if n == 0 else RR((p + m) / n)
        stddev = sqrt(density - mean**2)

        super().__init__(
            stddev=stddev,
            mean=mean,
            density=density,
            bounds=(-1, 1),
            n=n
        )

    def resize(self, new_n):
        """
        Return an altered distribution having a dimension `new_n`.
        Assumes `p` and `m` stay the same.
        """
        return SparseTernary(new_n, self.p, self.m)

    @property
    def hamming_weight(self):
        return self.p + self.m

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.SparseTernary(64, 8).support_size()
            6287341680214194176
        """
        n, p, m = len(self), self.p, self.m
        return ceil(binomial(n, p) * binomial(n - p, m) * RR(fraction))
