from pathlib import Path
import os,subprocess,resource,sys
benchmark_env=os.environ.copy()
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):benchmark_env[name]='1'
result=subprocess.run(['Rscript','--vanilla',str(Path(__file__).with_name('benchmark_tier.R'))]+sys.argv[1:],env=benchmark_env)
print('process_peak_rss_kib_on_linux',resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,flush=True)
sys.exit(result.returncode)
