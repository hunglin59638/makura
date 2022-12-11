# Makura: NCBI Genome downloader 

## Requirements

- rsync (linux command, required for downloading NCBI genoemes)
- python 3.8 (or greater)

## Installation

### Rsync
```
conda install -c conda-forge rsync
# or 
sudo apt install rsync
```

### Python packages

https://pypi.org/project/makura/

install from Pypi

```
pip install makura
```

install locally
```
python setup.py install
```

## Usage

Update the assembly summary and taxonomy information while first using.
```
makura update --assembly-source refseq
```
It's ok that you don't run the command, makura will automatically update if the assembly summary is not found.

Download bacteria and fungi genomes with complete assembly level in RefSeq database.  

```
makura download --group bacteria,fungi --assembly-level complete --assembly-source refseq --out_dir /path/to/dir
```


Print the records of genomes with JSON format
```
makura summary --accession GCF_016700215.2
```

Download genomes with selected taxids
```
makura download --taxid 2209
```

If you have many items to input, input a file contains lines is supported.
Example:
taxid_list.txt
```
61645
69218
550
```

```
makura download --taxid-list taxid_list.txt --out_dir /path/to/dir
```

Tips:

Running with multiple downloads in parallel is supported (Default: 4).  
We set the maximum is 8 to avoid NCBI blocks the downloads.  
```
makura download --group bacteria,fungi --parallel 4
```

While downloading the genomes, makura can check the MD5 checksum of them.
The MD5 values was stored named `md5checksums.txt` in output directory.


## Features in the future
- Creating minimap2 and bwa index using downloaded genomes.
- Downloading genomes by organism name, biosample, bioproject, etc.