## Install `gfortran`
It is possible to encounter errors when installing the R package `RtsneWithP` if `gfortran` is not already installed. The error is common among Mac with M-series chips. Below is the instructions to install `gfortran` for M-chip Mac users.
- Install `gcc` which includes `gfortran` with `brew install gcc`
- Create a file `~/.R/Makevars` (if it does not exist yet). You can run `mkdir -p ~/.R` to make the directory whether or not it exists. Then you can create the file using `touch ~/.R/Makevars`.
- Add the following lines to `~/.R/Makevars`.

`FC = /opt/homebrew/Cellar/gcc/11.3.0_2/bin/gfortran`
`F77 = /opt/homebrew/Cellar/gcc/11.3.0_2/bin/gfortran`
`FLIBS = -L/opt/homebrew/Cellar/gcc/11.3.0_2/lib/gcc/11`
