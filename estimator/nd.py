# -*- coding: utf-8 -*-

from dataclasses import dataclass

from sage.all import parent, RR, RealField, sqrt, pi, oo


def stddevf(sigma):
    """
    Gaussian width parameter σ → standard deviation.

    :param sigma: Gaussian width parameter σ

    EXAMPLE::

        sage: from estimator.nd import stddevf
        sage: stddevf(64.0)
        25.532...

        sage: stddevf(64)
        25.532...

        sage: stddevf(RealField(256)(64)).prec()
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

        sage: from estimator.nd import stddevf, sigmaf
        sage: n = 64.0
        sage: sigmaf(stddevf(n))
        64.000...

        sage: sigmaf(RealField(128)(1.0))
        2.5066282746310005024157652848110452530
        sage: sigmaf(1.0)
        2.50662827463100
        sage: sigmaf(1)
        2.50662827463100
        sage: sigmaf(1r)
        2.50662827463100

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
    bounds: tuple = None
    density: float = 1.0  # Hamming weight / dimension.
    tag: str = ""

    def __lt__(self, other):
        """
        We compare distributions by comparing their standard deviation.

        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.DiscreteGaussian(2.0) < ND.CentredBinomial(18)
            True
            sage: ND.DiscreteGaussian(3.0) < ND.CentredBinomial(18)
            False
            sage: ND.DiscreteGaussian(4.0) < ND.CentredBinomial(18)
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

            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.DiscreteGaussian(2.0) <= ND.CentredBinomial(18)
            True
            sage: ND.DiscreteGaussian(3.0) <= ND.CentredBinomial(18)
            True
            sage: ND.DiscreteGaussian(4.0) <= ND.CentredBinomial(18)
            False

        """
        try:
            return self.stddev <= other.stddev
        except AttributeError:
            return self.stddev <= other

    def __repr__(self):
        """
        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.DiscreteGaussianAlpha(0.01, 7681)
            D(σ=30.64, μ=0.00)

        """
        if self.n:
            return f"D(σ={float(self.stddev):.2f}, μ={float(self.mean):.2f}, n={int(self.n)})"
        else:
            return f"D(σ={float(self.stddev):.2f}, μ={float(self.mean):.2f})"

    def __hash__(self):
        """
        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: hash(ND(3.0, 1.0)) == hash((3.0, 1.0, None))
            True

        """
        return hash((self.stddev, self.mean, self.n))

    def __len__(self):
        """
        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: D = ND.SparseTernary(1024, p=128, m=128)
            sage: len(D)
            1024
            sage: round(len(D) * D.density)
            256

        """
        if self.n is not None:
            return self.n
        else:
            raise ValueError("Distribution has no length.")

    @staticmethod
    def DiscreteGaussian(stddev, mean=0, n=None):
        """
        A discrete Gaussian distribution with standard deviation ``stddev`` per component.

        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.DiscreteGaussian(3.0, 1.0)
            D(σ=3.00, μ=1.00)

        """
        return NoiseDistribution(
            stddev=RR(stddev), mean=RR(mean), n=n, bounds=(-oo, oo), tag="DiscreteGaussian"
        )

    @staticmethod
    def DiscreteGaussianAlpha(alpha, q, mean=0, n=None):
        """
        A discrete Gaussian distribution with standard deviation α⋅q/√(2π) per component.

        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.DiscreteGaussianAlpha(0.001, 2048)
            D(σ=0.82, μ=0.00)

        """
        stddev = stddevf(alpha * q)
        return NoiseDistribution.DiscreteGaussian(stddev=RR(stddev), mean=RR(mean), n=n)

    @staticmethod
    def CentredBinomial(eta, n=None):
        """
        Sample a_1, …, a_η, b_1, …, b_η and return Σ(a_i - b_i).

        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.CentredBinomial(8)
            D(σ=2.00, μ=0.00)

        """
        stddev = sqrt(eta / 2.0)
        # TODO: density
        return NoiseDistribution(
            stddev=RR(stddev), mean=RR(0), n=n, bounds=(-eta, eta), tag="CentredBinomial"
        )

    @staticmethod
    def Uniform(a, b, n=None):
        """
        Uniform distribution ∈ ``[a,b]``, endpoints inclusive.

        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.Uniform(-3, 3)
            D(σ=2.00, μ=0.00)
            sage: ND.Uniform(-4, 3)
            D(σ=2.29, μ=-0.50)

        """
        if b < a:
            raise ValueError(f"upper limit must be larger than lower limit but got: {b} < {a}")
        m = b - a + 1
        mean = (a + b) / RR(2)
        stddev = sqrt((m ** 2 - 1) / RR(12))

        if a <= 0 and 0 <= b:
            density = 1.0 / m
        else:
            density = 0.0

        return NoiseDistribution(
            n=n, stddev=stddev, mean=mean, bounds=(a, b), density=density, tag="Uniform"
        )

    @staticmethod
    def UniformMod(q, n=None):
        """
        Uniform mod ``q``, with balanced representation.

        EXAMPLE::

            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.UniformMod(7)
            D(σ=2.00, μ=0.00)
            sage: ND.UniformMod(8)
            D(σ=2.29, μ=-0.50)


        """
        a = -(q // 2)
        b = q // 2
        if q % 2 == 0:
            b -= 1
        return NoiseDistribution.Uniform(a, b, n=n)

    @staticmethod
    def SparseTernary(n, p, m=None):
        """
        Distribution of vectors of length ``n`` with ``p`` entries of 1 and ``m`` entries of -1, rest 0.

        EXAMPLE::
            sage: from estimator.nd import NoiseDistribution as ND
            sage: ND.SparseTernary(100, p=10)
            D(σ=0.45, μ=0.00, n=100)
            sage: ND.SparseTernary(100, p=10, m=10)
            D(σ=0.45, μ=0.00, n=100)
            sage: ND.SparseTernary(100, p=10, m=8)
            D(σ=0.42, μ=0.02, n=100)

        """
        if m is None:
            m = p
        mean = RR(p / n - m / n)
        stddev = RR(sqrt((p + m) / n))
        density = RR((p + m) / n)
        D = NoiseDistribution(
            stddev=stddev, mean=mean, density=density, bounds=(-1, 1), tag="SparseTernary", n=n
        )
        D.h = p + m
        return D