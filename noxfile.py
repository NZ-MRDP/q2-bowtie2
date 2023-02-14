import nox

# nox.options.sessions = ["test", "lint", "coverage", "typing"]
nox.options.sessions = ["lint", "typing"]


@nox.session(python=["3.8"])
def test(session):
    session.install("poetry")
    session.run("poetry", "install")
    session.run("pytest", "--cov-report=term-missing", "--cov=app", "app/tests")


@nox.session
def coverage(session):
    session.install("poetry")
    session.run("poetry", "install")
    session.run("coverage", "report")


@nox.session
def lint(session):
    session.install("poetry")
    session.run("poetry", "lock", "--no-update")
    session.run("poetry", "install")
    session.run("black", "./")
    session.run("flake8")
    # session.run("pytest", "--isort")


@nox.session
def typing(session):
    session.install("poetry")
    session.run("poetry", "install")
    session.run("pyright", ".")
