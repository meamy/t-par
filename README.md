**_Note:_** _This tool has been mostly supplanted by [Feynman](https://github.com/meamy/feynman)_

# T-par
Written by Matthew Amy

## Description
T-par is a quantum circuit optimizer utilizing sum-over-paths representations
to analyze and re-synthesize Clifford+T circuits.

The default algorithm is described in the paper [Polynomial-time T-depth
Optimization of Clifford+T circuits via Matroid Partitioning](https://arXiv.org/abs/1303.2042).
It creates a semantic representation of a quantum circuit where the phase of a
basis state in the output is given by a polynomial over the input and "path"
variables. This representation causes individual phase gates to be combined or
cancelled in the polynomial, reducing the total number of phase gates (notably,
T gates). Matroid partitioning is then used to optimally parallelize the
remaining phase gates during resynthesis.

## Building

To build T-par, run make in the top level folder.

tpar requires the following libraries:

* Boost

Boost should be available through your package manager. Additionally, 
your compiler needs to support the c++0x/c++11 standard, or otherwise
the code will likely require some (minor) modifications.

## Usage
Run T-par with
```
  ./t-par [options]
```

tpar takes a circuit in the .qc format (a description can be found in
the [QCViewer](https://github.com/aparent/QCViewer) repository) 
from standard input and outputs the
resulting .qc circuit to standard output. The circuit can only contain the 
single qubit gates H, P, P*, T, T*, X, Y, Z, and the two qubit tof (CNOT) gate.
It also accepts doubly controlled Z gates, i.e. Z a b c.

Demos are also available in the subfolder "demos". For example, to compute an
optimized 6-bit cucarro adder, use the command

```
  ./t-par < demos/cucarro_adder.qc
```

### Options
```
  -ancillae [0.., unbounded] - Specify the number of ancillae to add during
                               resynthesis. A natural number argument adds that
                               number of ancillae, while "unbounded" allows the
                               synthesizer to use as many ancillae as needed to
                               maximally parallelize phase gates

  -no-hadamard - Perform the T-par algorithm only on {CNOT, T} subcircuits. It
                 may provide better T-parallelization in some circuits

  -synth=[ADHOC,GAUSS,PMH] - Specify the synthesis method for linear reversible
                             circuits (i.e. {CNOT} circuits). ADHOC is an
                             informal (but correct) Gaussian elimination type 
                             algorithm, GAUSS is formal Gaussian elimination,
                             and PMH is the asymptotically optimal algorithm
                             described in "Optimal synthesis of linear
                             reversible circuits"

  -no-post-process - Turns off post processing of the synthesized circuit to
                     remove swap gates and trivial identities. Turning this off
                     may speed up synthesis for very large circuits

  -log - Display a log of the algorithm's process
```

This README is far from complete, so please feel free to email me at 
matt.e.amy@gmail.com if you have any questions or if you find any bugs.
