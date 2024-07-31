import time
import logging


class timeit(object):
    """measure time spend in context

    :class:`timeit` is a context manager (to be used with the :keyword:`with`
    statement) that records the execution time for the enclosed context block
    in :attr:`elapsed`.

    Attributes
    ----------
    elapsed : float
        Time in seconds that elapsed between entering
        and exiting the context.

    Example
    -------
    Use as a context manager::

       with timeit() as total:
          # code to be timed

       print(total.elapsed, "seconds")

    See Also
    --------
    :func:`time.time`

    """

    def __enter__(self):
        self._start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        end_time = time.time()
        self.elapsed = end_time - self._start_time
        # always propagate exceptions forward
        return False


# NOTE: modify this so that it can accept "http://localhost:8888" strings
# think of other edge cases as well- should be robust
def parse_host_port(filename):
    # Check if the format is correct
    parts = filename.split(":")
    if len(parts) == 2:
        host = parts[0]  # Hostname part
        try:
            port = int(parts[1])  # Convert the port part to an integer
            return (host, port)
        except ValueError:
            # Handle the case where the port is not a valid integer
            raise ValueError("Port must be an integer")
    else:
        # Handle the case where the format does not match "host:port"
        raise ValueError("Filename must be in the format 'host:port'")
