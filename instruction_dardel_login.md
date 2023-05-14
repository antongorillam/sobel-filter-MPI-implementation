# Instruction for loggin in Dardel

Create a Kerberos ticket and then ssh to Dardel:

``` bash
kinit -f whafu@NADA.KTH.SE
ssh -o GSSAPIDelegateCredentials=yes -o GSSAPIKeyExchange=yes -o GSSAPIAuthentication=yes whafu@dardel.pdc.kth.se
```

To allocate time-slot, run `salloc -n <#core> -t <#time> -p shared -A edu23.sf2568`, an example might be:

``` bash
salloc -n 1 -t 00:15:00 -p shared -A edu23.sf2568
CC sobel_filter_mpi.cpp -o sobel_filter_mpi.x
srun sobel_filter_mpi.x images/mat_files/8049.csv -05
```
