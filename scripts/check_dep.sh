#!/usr/bin/env bash

hash trim_galore 2>/dev/null || {
    echo "trim_galore error: no such command";
    echo "Please run: sudo apt install trim-galore";
    exit 1;
}

hash hisat2 2>/dev/null || {
    echo "hisat2 error: no such command";
    echo "Please run: sudo apt install hisat2";
    exit 1;
}

hash gffread 2>/dev/null || {
    echo "gffread error: no such command";
    echo "Please run: sudo apt install gffread";
    exit 1;
}

hash samtools 2>/dev/null || {
    echo "samtools error: no such command";
    echo "Please run: sudo apt isntall samtools";
    exit 1;
}

hash featureCounts 2>/dev/null || {
    echo "featureCounts error: no such command";
    echo "Please run: sudo apt install subread";
    exit 1;
}
