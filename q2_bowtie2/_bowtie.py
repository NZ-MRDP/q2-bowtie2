import os
import subprocess

import pandas as pd
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)


def build(reference_seqs: DNAFASTAFormat) -> Bowtie2IndexDirFmt:
    """build."""
    bowtie_database = Bowtie2IndexDirFmt()

    cmd = [
        "bowtie2-build",
        str(reference_seqs),
        os.path.join(str(bowtie_database), "bowtie2_index"),
    ]
    subprocess.run(cmd, check=True)
    return bowtie_database


def align(
    bowtie_database: Bowtie2IndexDirFmt,
    demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
) -> CasavaOneEightSingleLanePerSampleDirFmt:
    """align."""
    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    for _, fwd in df.itertuples():
        cmd = [
            "bowtie2",
            "-x",
            os.path.join(str(bowtie_database), "bowtie2_index"),
            "-U",
            str(fwd),
            "--un-gz",
            os.path.join(str(filtered_seqs), os.path.basename(fwd)),
        ]
        subprocess.run(cmd, check=True)
    return filtered_seqs
