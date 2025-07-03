# fastq_add_cell


This does open two or three fastq files and adds the read from the --cell fastq file as part of the r1 and r2 file's id, allowing any mapper to map the modified R1 and R2 fastq files and still recover the cell ID afterwards.

Should be pretty fast.

# Install

You need the Rust compiler and git.

## Download the repo and compile

```
git clone https://github.com/stela2502/fastq_add_cell.git
cd fastq_add_cell
cargo build -r
```

## Install from github

```
cargo install --git  https://github.com/stela2502/fastq_add_cell.git
```


# Usage


```
fastq_add_cell -h
Usage: fastq_add_cell [OPTIONS] --cell <CELL> --r1 <R1>

Options:
  -c, --cell <CELL>            FASTQ file for cell (required)
  -1, --r1 <R1>                FASTQ file for R1 (required)
  -2, --r2 <R2>                FASTQ file for R2 (optional)
      --from-char <FROM_CHAR>  Ommit some nulceotides from the start of the cell read?
      --to-char <TO_CHAR>      Stop reading from the cell read prematurely?
      --recomp                 Does the Cell read need to be reverse complemented?
  -h, --help                   Print help
  -V, --version                Print version
```

### Example

```
fastq_add_cell --cell <cell read fastq> --r1 <read 1> --r2 <read 2>
```

This will create new files from the read 1 and read 2 files and add a ``_cells_added`` to the fastq file names.
E.g. ``R1.fastq.gz`` will become ``R1_cells_added.fastq.gz``.
