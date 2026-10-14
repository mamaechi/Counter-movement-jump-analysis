"""Paired curve inference; samples are trials (subject) or subjects (group)."""
from itertools import product
import numpy as np
import spm1d


def paired_spm(right, left, alpha=.05, family_size=3):
    right,left=np.asarray(right,float),np.asarray(left,float)
    if right.shape != left.shape or right.ndim != 2 or right.shape[1] != 101:
        raise ValueError('Expected paired n × 101 arrays.')
    n=len(right)
    if n < 3:
        raise ValueError('Parametric SPM requires at least three pairs in this app; two trials support descriptive/exact comparisons only.')
    d=right-left
    if not np.isfinite(d).all() or np.any(d.std(0,ddof=1)<1e-10):
        raise ValueError('Zero difference variance at one or more nodes or non-finite values; SPM is undefined.')
    inference=spm1d.stats.ttest_paired(right,left).inference(alpha/family_size,two_tailed=True,interp=True)
    result=dict(method='paired parametric SPM{t}',n_pairs=n,alpha_family=alpha,
                family_size=family_size,alpha_plane=alpha/family_size,df=n-1,
                threshold=float(inference.zstar),reject=bool(inference.h0reject),
                clusters=[dict(start_percent=float(c.endpoints[0]),end_percent=float(c.endpoints[1]),
                               p_cluster=float(c.P)) for c in inference.clusters])
    return result,inference


def exact_max_t(right,left,alpha=.05,family_size=3):
    """Exhaustive paired sign flips with max-|t| correction across the 101 nodes.

    A separate test per plane with Bonferroni correction across planes. This is
    a permutation field test, not the parametric random-field SPM cluster test.
    Assumes exchangeable signs of independent trial differences under the null.
    """
    right,left=np.asarray(right,float),np.asarray(left,float)
    if right.shape != left.shape or right.ndim != 2:
        raise ValueError('Expected matching paired arrays.')
    d=right-left
    n=len(d)
    if d.shape != (n,101) or n<2 or n>12 or not np.isfinite(d).all():
        raise ValueError('Exact comparison supports 2–12 finite paired curves of 101 nodes.')
    def statistic(a):
        mean=a.mean(0); se=a.std(0,ddof=1)/np.sqrt(n)
        return np.divide(mean,se,out=np.where(np.abs(mean)>1e-12,np.copysign(np.inf,mean),0.),where=se>1e-12)
    observed=statistic(d)
    maxima=np.array([np.max(np.abs(statistic(d*np.array(signs)[:,None])))
                     for signs in product((-1,1),repeat=n)])
    p=np.mean(maxima[:,None]>=np.abs(observed)[None,:]-1e-12,axis=0)
    return dict(method='exact paired sign-flip max-|t|',n_pairs=n,permutations=2**n,
                minimum_two_sided_p=2/(2**n),alpha_plane=alpha/family_size,
                p_global=float(np.mean(maxima>=np.max(np.abs(observed))-1e-12)),
                reject=bool(np.any(p<=alpha/family_size))),p
