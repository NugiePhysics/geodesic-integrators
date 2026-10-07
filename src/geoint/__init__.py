"""Numerical integrators for black-hole geodesics.

The package compares explicit Runge-Kutta, Gauss-Legendre and Tao integrators on Schwarzschild
geodesics in two formulations: the second-order geodesic equation in ``(x^mu, u^mu)`` and
Hamilton's equations in ``(x^mu, p_mu)``. The theory is in ``docs/theory/formulations.md``.
"""

__version__ = "1.0.0"

__all__ = ["__version__"]
