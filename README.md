Generalized repeat complexity per chromosome.  Rk(c) = sum { f(x) : x \in c and f(x) > 1 } /(|c|-k+1), 
where f(x) is the frequency of x in the entire genome.  Note that x is a substring of chromosome c.

Example:

```
	cd usage
	go run compute_complexity.go chr1.fasta chr2.fasta chr3.fasta
```
