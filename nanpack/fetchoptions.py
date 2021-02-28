"""Main script to obtain available options from the command terminal."""

import nanpack.backend.fetchoptions as fopt


def fetch(what):
    """Show allowed user input options for various parameters.

    Parameter
    ---------
    what: str
        Provide a text string. Choose one from the below list.
        - "model" for allowed model inputs
        - "limiter-functions" for available TVD limiter function inputs
    """
    f = fopt._FetchOptions()
    if what.lower() == "model":
        print(*f.ModelOptions(), sep="\n")

    elif what.lower() == "limiter-functions":
        print(*f.TVDLimiterFunctionOptions(), sep="\n")


if __name__ == "__main__":
    import sys
    fetch(sys.argv[1])
