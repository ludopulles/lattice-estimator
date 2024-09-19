# -*- coding: utf-8 -*-

from dataclasses import dataclass

from sage.all import binomial, ceil, exp, log, oo, parent, pi, RealField, RR, sqrt


def stddevf(sigma):
    """
    Gaussian width parameter σ → standard deviation.

    :param sigma: Gaussian width parameter σ

    EXAMPLE::

        >>> from estimator.nd import stddevf
        >>> stddevf(64.0)
        25.532...

        >>> stddevf(64)
        25.532...

        >>> stddevf(RealField(256)(64)).prec()
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

        >>> from estimator.nd import stddevf, sigmaf
        >>> n = 64.0
        >>> sigmaf(stddevf(n))
        64.000...

        >>> sigmaf(RealField(128)(1.0))
        2.5066282746310005024157652848110452530
        >>> sigmaf(1.0)
        2.506628274631...
        >>> sigmaf(1)
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
    stddev: float
    mean: float = 0
    n: int = None
    bounds: tuple = (-oo, oo)
    density: float = 1.0  # Hamming weight / dimension.

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
    def is_sparse(self):
        """
        Whether the density of the distribution is < 1/2.
        Note: 1/2 might be considered somewhat arbitrary.
        """
        # NOTE: somewhat arbitrary
        return self.density < 0.5


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

        >>> from estimator.nd import DiscreteGaussian
        >>> DiscreteGaussian(3.0, 1.0)
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

    def support_size(self, fraction=1.0, n=None):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator.nd import DiscreteGaussian
            >>> DiscreteGaussian(1.0, n=128).support_size(fraction=.99)
            ???
        """
        # We will treat this noise distribution as bounded with failure probability `1 - fraction`.
        t = self.gaussian_tail_bound
        p = self.gaussian_tail_prob

        if p**n < fraction:
            raise NotImplementedError(
                f"TODO(nd.support-size): raise t. {RR(p ** n)}, {n}, {fraction}"
            )

        b = 2 * t * sigmaf(self.stddev) + 1
        return (2 * b + 1)**n


class DiscreteGaussianAlpha(DiscreteGaussian):
    """
    A discrete Gaussian distribution with standard deviation α⋅q/√(2π) per component.

    EXAMPLE::

        >>> from estimator.nd import DiscreteGaussianAlpha
        >>> DiscreteGaussianAlpha(0.001, 2048)
        D(σ=0.82)
    """
    def __init__(self, alpha, q, mean=0, n=None):
        super().__init__(RR(stddevf(alpha * q)), RR(mean), n)


class CenteredBinomial(NoiseDistribution):
    """
    Sample a_1, …, a_η, b_1, …, b_η and return Σ(a_i - b_i).

    EXAMPLE::

        >>> import estimator.ND
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

    def support_size(self, fraction=1.0, n=None):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> CenteredBinomial(3, 10).support_size()
            282475249
            >>> ND.CenteredBinomial(3, 10).support_size(fraction=.99)
            279650497
        """
        if not n:
            n = len(self)

        # TODO: this might be suboptimal/inaccurate for binomial distribution
        a, b = self.bounds
        return ceil(RR(fraction) * (b - a + 1)**n)


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

    def support_size(self, fraction=1.0, n=None):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.Uniform(-3,3, 64).support_size(fraction=.99)
            1207562882759477428726191443614714994252339953407098880
        """
        if not n:
            n = len(self)

        # TODO: this might be suboptimal/inaccurate for binomial distribution
        a, b = self.bounds
        return ceil(RR(fraction) * (b - a + 1)**n)


class UniformMod(Uniform):
    """
    Uniform mod ``q``, with balanced representation.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.UniformMod(7)
        D(σ=2.00)
        >>> ND.UniformMod(8)
        D(σ=2.29, μ=-0.50)
    """
    def __init__(self, q, n=None):
        half_q = q // 2
        super().__init__(half_q - q + 1, half_q, n=n)


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

        # Yes, n=0 might happen in the dual attack!
        mean = 0 if n == 0 else RR((p - m) / n)
        density = 0 if n == 0 else RR((p + m) / n)
        stddev = sqrt(density - mean**2)
        # stddev = sqrt(p / n * (1 - mean)**2 +
        #               m / n * (-1 - mean)**2 +
        #               (n - (p + m)) / n * (mean)**2)

        super().__init__(
            stddev=stddev,
            mean=mean,
            density=density,
            bounds=(-1, 1),
            n=n
        )

    @property
    def hamming_weight(self):
        return self.p + self.m

    def support_size(self, fraction=1.0, n=None):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.SparseTernary(64, 8).support_size()
            32016101348447354880
        """
        if not n:
            n = len(self)

        h = self.hamming_weight
        # TODO: this is assuming that the non-zero entries are uniform over {-1,1}
        # need p and m for more accurate calculation
        return ceil(2**h * binomial(n, h) * RR(fraction))
