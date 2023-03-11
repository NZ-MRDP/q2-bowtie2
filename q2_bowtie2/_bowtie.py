import os
import subprocess
import tempfile

import pandas as pd
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_types_genomics.per_sample_data import BAMDirFmt


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


# TODO: Refactor redundant code throughout
def align_paired(
    bowtie_database: Bowtie2IndexDirFmt,
    demultiplexed_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
) -> (CasavaOneEightSingleLanePerSampleDirFmt, CasavaOneEightSingleLanePerSampleDirFmt, BAMDirFmt,):  # type: ignore
    """align_paired.

    Parameters
    ----------
    bowtie_database : Bowtie2IndexDirFmt
        bowtie_database
    demultiplexed_sequences : SingleLanePerSampleSingleEndFastqDirFmt
        demultiplexed_sequences

    Returns
    -------
    (CasavaOneEightSingleLanePerSampleDirFmt, CasavaOneEightSingleLanePerSampleDirFmt)

    """
    """align."""
    aligned_filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    unaligned_filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    bowtie_alignment = BAMDirFmt()
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    for sample_id, fwd, rev in df.itertuples():
        # Bowtie renames fwd and rev sequences with .1 and 2.
        # so we give them a temporary name and change it back to what it should be later
        aligned_path = os.path.join(str(aligned_filtered_seqs), "temp_name")
        unaligned_path = os.path.join(str(unaligned_filtered_seqs), "temp_name")
        cmd = [
            "bowtie2",
            "-x",
            os.path.join(str(bowtie_database), "bowtie2_index"),
            "-1",
            str(fwd),
            "-2",
            str(rev),
            "--un-conc-gz",
            aligned_path,
            "--al-conc-gz",
            unaligned_path,
        ]

        with tempfile.NamedTemporaryFile() as temp:
            subprocess.run(cmd, check=True, stdout=temp)
            with open(os.path.join(str(bowtie_alignment), f"{sample_id}.bam"), "w") as bowtie_file:
                subprocess.run(["samtools", "view", "-bS", temp.name], check=True, stdout=bowtie_file)

        # rename the files what they should be
        os.rename(aligned_path + ".1", os.path.join(str(aligned_filtered_seqs), os.path.basename(fwd)))
        os.rename(aligned_path + ".2", os.path.join(str(aligned_filtered_seqs), os.path.basename(rev)))
        os.rename(unaligned_path + ".1", os.path.join(str(unaligned_filtered_seqs), os.path.basename(fwd)))
        os.rename(unaligned_path + ".2", os.path.join(str(unaligned_filtered_seqs), os.path.basename(rev)))
    return aligned_filtered_seqs, unaligned_filtered_seqs, bowtie_alignment


def align_single(
    bowtie_database: Bowtie2IndexDirFmt,
    demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
) -> (CasavaOneEightSingleLanePerSampleDirFmt, CasavaOneEightSingleLanePerSampleDirFmt, BAMDirFmt):  # type: ignore
    """align_single.

    Parameters
    ----------
    bowtie_database : Bowtie2IndexDirFmt
        bowtie_database
    demultiplexed_sequences : SingleLanePerSampleSingleEndFastqDirFmt
        demultiplexed_sequences

    Returns
    -------
    (CasavaOneEightSingleLanePerSampleDirFmt, CasavaOneEightSingleLanePerSampleDirFmt)

    """
    aligned_filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    unaligned_filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    bowtie_alignment = BAMDirFmt()
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    for sample_id, sample_path in df.itertuples():
        cmd = [
            "bowtie2",
            "-x",
            os.path.join(str(bowtie_database), "bowtie2_index"),
            "-U",
            str(sample_path),
            "--un-gz",
            os.path.join(str(unaligned_filtered_seqs), os.path.basename(sample_path)),
            "--al-gz",
            os.path.join(str(aligned_filtered_seqs), os.path.basename(sample_path)),
        ]
        with tempfile.NamedTemporaryFile() as temp:
            subprocess.run(cmd, check=True, stdout=temp)
            with open(os.path.join(str(bowtie_alignment), f"{sample_id}.bam"), "w") as bowtie_file:
                subprocess.run(["samtools", "view", "-bS", temp.name], check=True, stdout=bowtie_file)

    return aligned_filtered_seqs, unaligned_filtered_seqs, bowtie_alignment
