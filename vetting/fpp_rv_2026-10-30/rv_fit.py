"""Frozen ephemeris candidate RV fit; fixed published GP; linear Gaussian likelihood with bounded K priors."""
from pathlib import Path
import json,hashlib,sys,ast,time
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.stats import truncnorm
from scipy.special import ndtr
from scipy.integrate import cumulative_trapezoid
from astropy.constants import G,M_sun,M_earth
from net import space
P=18.908226;T=59904.5666 # RJD, source-note BJD minus 2400000
DATA=Path('evidence/rv/extracted/TOI-4342_RV/data/TOI4342_RV_Aug25.dat')

def sine(t,p,tc):return -np.sin(2*np.pi*(t-tc)/p)

def gp_qp(dt,A=8.6,O=.21,Prot=14.691,decay=79.):
 return A*A*np.exp(-np.sin(np.pi*dt/Prot)**2/(2*O*O)-dt*dt/(2*decay*decay))

def main():
 space();a=np.loadtxt(DATA);t,y,e=a[:,:3].T;order=np.argsort(t);t,y,e=t[order],y[order],e[order];dt=t[:,None]-t;C=gp_qp(dt)+np.diag(e*e+2.22**2)
 # ESP covariance replacement is activated only after source-defined exact kernel is locally available.
 kernel='Exact quasi-periodic kernel Eq. 1; ESP approximation comparison UNVERIFIED'
 if Path('evidence/rv/esp_dense.py').exists():
  env={'np':np};exec(Path('evidence/rv/esp_dense.py').read_text(),env);C=env['esp_dense'](dt,sig=8.6,P=14.691,rho=79.,eta=.21,nharm=4)+np.diag(e*e+2.22**2);kernel='S+LEAF ESP, n_harmonics=4, exact source-defined dense covariance'
 cf=cho_factor(C,lower=True);out=dict(data_source=str(DATA),data_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),n_rvs=len(t),baseline_d=float(t.max()-t.min()),candidate=dict(P_d=P,T0_RJD=T,e=0,K_prior_uniform_mps=[0,20]),GP=dict(amplitude_mps=8.6,coherence=.21,Prot_d=14.691,decay_d=79.,jitter_mps=2.22,kernel=kernel,fixed_published_hyperparameters=True),cases={},complete=False)
 for model in ['b_c_GP_candidate','b_c_47p5_GP_candidate']:
  cols=[np.ones(len(t)),sine(t,5.5382592,58654.53479),sine(t,10.688662,58659.3486),sine(t,P,T)];names=['gamma','K_b','K_c','K_18p9']
  if '47p5' in model:cols.extend([np.sin(2*np.pi*(t-59800)/47.5),np.cos(2*np.pi*(t-59800)/47.5)]);names.extend(['rv47_sin','rv47_cos'])
  X=np.column_stack(cols);ciX=cho_solve(cf,X);cov=np.linalg.inv(X.T@ciX);mean=cov@(X.T@cho_solve(cf,y));rng=np.random.default_rng(20261030);draws=rng.multivariate_normal(mean,cov,300000)
  # All noncandidate linear nuisance priors: uniform gamma, nonnegative b,c amplitudes in [0,10].
  good=(draws[:,1]>=0)&(draws[:,1]<=10)&(draws[:,2]>=0)&(draws[:,2]<=10)&(draws[:,3]>=0)&(draws[:,3]<=20);post=draws[good]
  if len(post)<10000:raise RuntimeError('Insufficient accepted posterior draws')
  k=post[:,3];percent=np.quantile(k,[.025,.16,.5,.84,.95,.975]);rr=dict(parameter_names=names,unconstrained_mean=mean.tolist(),unconstrained_covariance=cov.tolist(),n_proposed=len(draws),n_posterior=len(post),K_median_mps=float(percent[2]),K_68_mps=percent[[1,3]].tolist(),K_95_upper_mps=float(percent[4]),K_equal_tail_95_mps=percent[[0,5]].tolist(),Pr_K_below1=float(np.mean(k<1)),b_median=float(np.median(post[:,1])),c_median=float(np.median(post[:,2])),offset_median=float(np.median(post[:,0])),bounded_prior_note='Uniform b,c K [0,10], candidate K [0,20], gamma unbounded; all conditionally Gaussian. Sensitivity 47.5 d terms use broad signed linear quadratures, not the paper third-planet posterior.')
  np.savez_compressed('results/rv_'+model+'.npz',t_RJD=t,rv=y,err=e,covariance=C,design=X,unconstrained_mean=mean,unconstrained_covariance=cov,posterior_K=k[:50000]);out['cases'][model]=rr;print(model,rr['K_median_mps'],rr['K_68_mps'],rr['K_95_upper_mps'],rr['Pr_K_below1'],flush=True)
  # Check native package posterior planet-b model sign and time convention at fixed published parameters.
 pred=3.78*sine(t,5.5382592,58654.53479);p=Path('evidence/rv/extracted/TOI-4342_RV/model_files/RV_data_radial_velocities_b.dat');tab=np.loadtxt(p);out['sign_convention_check_max_abs_mps']=float(np.max(np.abs(pred-tab[:,6])))
 expected=[]
 for m in [4,8]:
  kval=((2*np.pi*G.value/(P*86400))**(1/3))*m*M_earth.value/((.587*M_sun.value+m*M_earth.value)**(2/3));expected.append(dict(mass_earth=m,K_mps=float(kval)))
 out['expected_K_circular_edge_on']=expected;out['bar']='FAIL real tension' if out['cases']['b_c_GP_candidate']['K_95_upper_mps']<1 else 'PASS no <1 m/s 95% upper-limit tension';out['complete']=True;out['verification_status']='VERIFIED conditional fixed-hyperparameter RV likelihood' if 'ESP,' in kernel else 'UNVERIFIED exact ESP approximation; QP cross-check only';Path('results/rv.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
