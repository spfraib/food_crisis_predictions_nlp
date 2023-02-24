#!/bin/bash
conda activate py27
mkdir /scratch/ab7325/frame-models/logs/$1 
python -m sesame.targetid --mode train --model_name targetid-manual --output_dir $1 
python -m sesame.frameid --mode train --model_name frameid-manual --output_dir $1
python -m sesame.argid --mode train --model_name argid-manual --output_dir $1
