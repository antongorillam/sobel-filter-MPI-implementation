# Instruction for loggin in Dardel

Create a Kerberos ticket and then ssh to Dardel:
```
kinit -f whafu@NADA.KTH.SE
ssh -o GSSAPIDelegateCredentials=yes -o GSSAPIKeyExchange=yes -o GSSAPIAuthentication=yes whafu@dardel.pdc.kth.se
```

To allocate time-slot, run `salloc -n <#core> -t <#time> -p shared -A edu23.sf2568`, an example might be:
```
salloc -n 8 -t 00:10:00 -p shared -A edu23.sf2568
cc odd_even_transposition_sort.c -o odd_even_transposition_sort.x
srun ./odd_even_transposition_sort.x 10000 -05
```

