# TP-Bio-Virus
This project studies an unknown virus genome in order to theorize what it does.

## Install
### Prerequisites
- You need Python 3
- To create a python virtual environment in your computer you need to have `virtualenv`.

### Creating virtual environment
You may need to add execution permissions to `initVenv.sh`. Then you run
```
./initVenv.sh
```
This will build the virtual environment with every dependency needed for this project

## Download swissprot
To download swissprot db needed for exercise 2, you need to run `download_swissprot.sh` in the `blast/data` file

## Run
To run the project you need to start the virtual environment. You can do this with the following command
```
source openVenv.sh
```

### Translate genome to proteins
```
python proteins_from_gen.py
```

#### Local Blast
```
./local_blast.sh
```

#### Remote Blast
```
python remote_blast.py
```

### Sequence Multialignment
```
python multialignment.py
```
