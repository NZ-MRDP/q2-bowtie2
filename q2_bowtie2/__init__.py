"""Bowtie python library."""
from ._bowtie import align_paired, align_single, build
from ._type import BowtieReadStatistics

__version__ = "0.0.0"

__all__ = ["build", "align_single", "align_paired", "BowtieReadStatistics"]
