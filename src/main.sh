#!/bin/bash
    for alpha in 0.5 
    do
        for chunk in {1..100}
        do
                qsub -v alpha=$alpha,chunk=$chunk run_model.sh
        done
    done


