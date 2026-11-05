import sys
from table2 import dev_full
T=float(sys.argv[1]); Bt=float(sys.argv[2]); N=int(sys.argv[3])
dm,dw,nf=dev_full(0.3,Bt,0.3,T,N=N)
print('T=%g K B=%.0e N=%d dM=%.4f dWF=%.3e nF=%d'%(T,Bt,N,dm,dw,nf),flush=True)
