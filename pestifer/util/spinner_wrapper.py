# Author: ChatGPT 5
# Assistant Author: Cameron F. Abrams <cfa22@drexel.edu>

import sys, threading, time
from functools import wraps
from itertools import cycle

# guards against nested spinners writing to the same terminal line; only the outermost
# with_spinner animates, so an inner cache build inherits the outer spinner rather than
# garbling the line.  Toggled only from the main thread (the spinner runs in a daemon thread).
_spinner_active = False


def with_spinner(text: str = "Working..."):
    """Decorator: show a spinner until the wrapped function returns/raises.

    The animation is shown only when standard output is an interactive terminal, and never
    nested; otherwise (piped/redirected output, or already inside another spinner) the wrapped
    function simply runs with no spinner output.
    """
    def decorator(fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            global _spinner_active
            if _spinner_active or not sys.stdout.isatty():
                return fn(*args, **kwargs)
            _spinner_active = True
            stop = threading.Event()
            spinner = cycle("|/-\\")
            def run_spinner():
                start = time.time()
                while not stop.is_set():
                    s = next(spinner)
                    elapsed = time.time() - start
                    sys.stdout.write(f"\r{s} {text}  [{elapsed:5.1f}s]")
                    sys.stdout.flush()
                    time.sleep(0.09)
                # clear line
                sys.stdout.write("\r" + " " * 80 + "\r")
                sys.stdout.flush()

            t = threading.Thread(target=run_spinner, daemon=True)
            t.start()
            try:
                return fn(*args, **kwargs)
            finally:
                stop.set()
                t.join()
                _spinner_active = False
        return wrapper
    return decorator
