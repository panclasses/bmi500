# Environment Setup

1. either do 
```
conda create -n bmi500 python=3.12.4
```

or 

```
python3 -m venv bmi500
```

2. activate environment
```
conda activate bmi500
```

or 

```
source bmi500/bin/activate
```

3. Install required packages
```
pip install -r requirements.txt
```

It may be required to install g++ and cmake in order to install python louvain


4. run scanpy_pbmc.py
```
python scanpy_pbmc.py --data-dir {data_dir} --data-set {data_set} --out-dir {output_dir} --num-threads {num_thread}
```

`data_dir` is the top level directory in which the data_set subdirectories resides.  defaults to `data`

`output_dir` is where the output file is placed. defaults to `data`

`data_set` is one of 'pbmc3k', 'pbmc6k', and 'pbmc10k'.

