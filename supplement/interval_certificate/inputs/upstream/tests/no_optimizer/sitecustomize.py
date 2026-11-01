"""Replay dependency guard: an optimizer import is an error, not a fallback."""
import builtins
_original_import=builtins.__import__
_FORBIDDEN={'scipy','highspy','cvxpy','pulp','ortools'}
def _no_optimizer(name,*args,**kwargs):
    if name.split('.')[0]in _FORBIDDEN:
        raise ImportError('Optimizer import deliberately disabled during exact replay: '+name)
    return _original_import(name,*args,**kwargs)
builtins.__import__=_no_optimizer
