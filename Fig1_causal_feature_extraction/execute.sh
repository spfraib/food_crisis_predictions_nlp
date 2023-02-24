#!/bin/bash
conda activate py2
mkdir /scratch/ab7325/frame-models/logs/v11_$1 
python -m sesame.targetid --mode predict --model_name targetid-manual --raw_input /scratch/ab7325/Factiva/extraction-v11/avro/all_sentences_with_lines_split_$1 --output_dir v11_$1 
python -m sesame.filter --mode filter --model_name frameid-manual --raw_input /scratch/ab7325/frame-models/logs/v11_$1/predicted-targets.conll --output_dir v11_$1
python -m sesame.frameid --mode predict --model_name frameid-manual --raw_input /scratch/ab7325/frame-models/logs/v11_$1/filtered-frames.conll --output_dir v11_$1
python -m sesame.argid --mode predict --model_name argid-manual --raw_input /scratch/ab7325/frame-models/logs/v11_$1/predicted-frames.conll --output_dir v11_$1
python -m sesame.filter --mode extract --model_name argid-manual --raw_input /scratch/ab7325/frame-models/logs/v11_$1/predicted-args.conll --output_dir v11_$1
python -m sesame.get_embedding  --output_dir v11_$1
