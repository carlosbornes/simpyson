"""SimPYson: A Pythonic interface for SIMPSON solid-state NMR simulations."""

from __future__ import annotations

from importlib.metadata import version

from simpyson.calculator import SimpCalc, simulate_spectrum
from simpyson.io import read_simp, write_simp
from simpyson.simpy import Simpy

__version__ = version("simpyson")

__all__ = [
    "Simpy",
    "SimpCalc",
    "simulate_spectrum",
    "read_simp",
    "write_simp",
]
