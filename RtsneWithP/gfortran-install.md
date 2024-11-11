## Install `gfortran`
It is possible to encounter errors when installing the R package `RtsneWithP` if `gfortran` is not already installed. The error is common among Mac with M-series chips. Below is the instruction to install `gfortran` for M-chip Mac users.
- Install `gcc` which includes `gfortran` with `brew install gcc`
- Create a file `~/.R/Makevars` (if it does not exist yet). You can run `mkdir -p ~/.R` to make the directory whether or not it exists. Then you can create the file using `touch ~/.R/Makevars`.
- Add the following lines to `~/.R/Makevars`.
>>`FC = /opt/homebrew/Cellar/gcc/11.3.0_2/bin/gfortran`

>>`F77 = /opt/homebrew/Cellar/gcc/11.3.0_2/bin/gfortran`

>>`FLIBS = -L/opt/homebrew/Cellar/gcc/11.3.0_2/lib/gcc/11`
- Run `brew list --versions gcc` to check your `gcc` version and change `gcc` version `11.3.0_2` above to whatever your `gcc` version is.
- Restart `R` and install `RtsneWithP` again.

### References
[1] https://stackoverflow.com/questions/69639782/installing-gfortran-on-macbook-with-apple-m1-chip-for-use-in-r

[2] https://github.com/Rdatatable/data.table/issues/4437
