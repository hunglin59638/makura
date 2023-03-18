# Makura: NCBI Genome downloader 

## Requirements

- rsync (linux command, required for downloading NCBI genomes)
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

install from docker
```
docker pull hunglin59638/makura
```

## Usage

Update the assembly summary and taxonomy information while first using.
```
makura update
```
It's ok that you don't run the command, makura will automatically update if the assembly summary is not found.

Download bacteria and fungi genomes with complete assembly level in RefSeq database.  

```
makura download --group bacteria,fungi --assembly-level complete --assembly-source refseq --out_dir /path/to/dir
```


Print the records of genomes with JSON lines format, default is TAB
```
makura summary --accession GCF_016700215.2 --as-json-lines
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
The MD5 values was stored to a file named `md5checksums.txt` in output directory.

## Developing function
Using the RESTful API to get assembly summary
1. run the API server
```
docker run --rm -p 5000:5000 hunglin59638/makura:1.1.0 makura api --port 5000
```
2. get the summary of assembly accessions
```
curl http://localhost:5000/summary?accessions=GCA_002287175.1,GCA_000762265.1
```
## Features in the future
- Creating minimap2 and bwa index using downloaded genomes.
- Downloading genomes by organism name, biosample, bioproject, etc.