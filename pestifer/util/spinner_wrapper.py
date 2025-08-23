# Author: ChatGPT 5
# Assistant Author: Cameron F. Abrams <cfa22@drexel.edu>

import sys, threading, time
from functools import wraps
from itertools import cycle

def with_spinner(text: str = "Working..."):
    """Decorator: show a spinner until the wrapped function returns/raises."""
    def decorator(fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
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
        return wrapper
    return decorator
