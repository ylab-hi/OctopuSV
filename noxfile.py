import nox


@nox.session(python="3.10")
def tests(session):
    session.run("poetry", "install", external=True)
    session.run(
        "poetry",
        "run",
        "pytest",
        "--cov",
        "--cov-report=html",
        external=True,
    )


@nox.session(python="3.10")
def pre_commit(session):
    session.run("poetry", "install", external=True)
    session.run(
        "poetry",
        "run",
        "pre-commit",
        "run",
        "--all-files",
        external=True,
    )


@nox.session(python="3.10")
def mypy(session):
    session.run("poetry", "install", external=True)
    session.run(
        "poetry",
        "run",
        "mypy",
        ".",
        external=True,
    )


@nox.session(python="3.10")
def safety(session):
    session.run("poetry", "install", external=True)
    session.run(
        "poetry",
        "run",
        "safety",
        "check",
        external=True,
    )


@nox.session(name="docs-build", python="3.10")
def docs_build(session):
    session.run("poetry", "install", external=True)
    session.run(
        "poetry",
        "run",
        "sphinx-build",
        "-b",
        "html",
        "docs",
        "docs/_build",
        external=True,
    )