import nox

@nox.session(python="3.10")
def tests(session):
    session.run("poetry", "install", external=True)
    session.run("pytest", "--cov", "--cov-report=html")
