"""QIIME 2 plugin for running Bowtie2."""

from setuptools import find_packages, setup

setup(
    name="q2-bowtie2",
    version="0.0.0",
    packages=find_packages(),
    author="John Chase",
    author_email="jhch@novozymes.com",
    description="QIIME2 plugin for running Bowtie2",
    entry_points={"qiime2.plugins": ["q2-bowtie2=q2_bowtie2.plugin_setup:plugin"]},
)
