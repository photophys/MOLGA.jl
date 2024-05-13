# Installation

MOLGA can be installed using [Julia's package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/). From the the Julia REPL, type `]` to start the package manager mode and run

```
pkg> update
pkg> add MOLGA
```

This package is an interplay between our genetic algorithm and third-party quantum-chemical calculations. As described [here](../features/interfaces.md), MOLGA currently supports Gaussian and XTB.

!!! info
    
    If the installation paths on your system differ from the defaults, you'll need to set up environment variables. This ensures MOLGA can locate the required executable files.
    
    You can specify these paths in an environment variable file (located in the same directory where you run MOLGA). Ensure you define the path(s) for the program(s) you intend to use, as specified in the [parameter file](../parameters/input-file.md).
    
    ```bash
    # .env file (showing defaults)
    GAUSSIAN=g16
    XTB=xtb
    ```
