# Cryptanalysis of the Legendre PRF - Implementation

This repository contains an implementation of the attack from the paper

> [Cryptanalysis of the Legendre PRF and Generalizations](https://eprint.iacr.org/2019/1357)

by *Ward Beullens, Tim Beyne, Aleksei Udovenko, and Giuseppe Vitto*.

The code can be used to break **Challenge 2** of the [Legendre PRF Bounties](https://legendreprf.org/bountyinstances) in under 1500 CPU-hours. For more details, please refer to the paper.

The code can be run with the following command:

```
$ make threads=24 target=P74
```

- `threads` argument defines the amount of threads to be used in the second step of the attack;
- `target` argument can be one of **P40, P64, P74, P84**.

It requires a C++ compiler to be installed. Clang++ is recommended.
Furthermore, *libgmp* must be installed.
