__all__ = ["run_magical"]

def run_magical(*args, **kwargs):
    from .magical import run_magical as _run_magical
    return _run_magical(*args, **kwargs)
