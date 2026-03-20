"""Nox sessions for q2-bowtie2."""

import nox

nox.options.sessions = ["lint", "typing"]


@nox.session(python=["3.8"])
def test(session):
    """Run the test suite."""
    session.run("uv", "sync", "--dev", external=True)
    session.run("pytest", "--cov-report=term-missing", "--cov=q2_bowtie2", "q2_bowtie2/tests")


@nox.session
def coverage(session):
    """Report test coverage."""
    session.run("uv", "sync", "--dev", external=True)
    session.run("coverage", "report")


@nox.session
def lint(session):
    """Run formatting and lint checks."""
    session.run("uv", "sync", "--dev", external=True)
    session.run("ruff", "format", "./")
    session.run("ruff", "check", "./")


@nox.session
def typing(session):
    """Run static type checks."""
    session.run("uv", "sync", "--dev", external=True)
    session.run("pyright", ".")
