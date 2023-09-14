import gzip
import os
import re
import subprocess
import tempfile
from typing import Optional

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
    threads: int = 1,
    very_sensitive: bool = False,
) -> (CasavaOneEightSingleLanePerSampleDirFmt, CasavaOneEightSingleLanePerSampleDirFmt, BAMDirFmt,):  # type: ignore
    """align_paired.

    Parameters
    ----------
    bowtie_database : Bowtie2IndexDirFmt
        bowtie_database
    demultiplexed_sequences : SingleLanePerSampleSingleEndFastqDirFmt
        demultiplexed_sequences
    threads : int, optional
        number of alignment threads to launch. Defaults to 1.
    very_sensitive : bool, optional
        Whether to run bowtie2 with the very sensitive preset. Defaults to False.

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
            "--threads",
            str(threads),
        ]
        if very_sensitive:
            cmd.append("--very-sensitive")

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


def extract_first_number(text: str) -> Optional[int]:
    """
    Extract first number in string.

    Args:
        text (str): String wherefrom a number is extracted.

    Returns:
        Optional[int]: Number extracted.
    """
    match = re.search(r"\d+", text)
    if match:
        return int(match.group())
    else:
        return None


def align_single(
    bowtie_database: Bowtie2IndexDirFmt,
    demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    threads: int = 1,
    save_alignment: bool = False,
    very_sensitive: bool = False,
) -> (
    CasavaOneEightSingleLanePerSampleDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
    BAMDirFmt,
    pd.DataFrame,
):  # type: ignore
    """align_single.

    Parameters
    ----------
    bowtie_database : Bowtie2IndexDirFmt
        bowtie_database
    demultiplexed_sequences : SingleLanePerSampleSingleEndFastqDirFmt
        demultiplexed_sequences
    threads : int, optional
        number of alignment threads to launch. Defaults to 1.
    save_alignment : bool, optional
        Whether to save alignment files. Defaults to False.
    very_sensitive : bool, optional
        Whether to run bowtie2 with the very sensitive preset. Defaults to False.

    Returns
    -------
    (CasavaOneEightSingleLanePerSampleDirFmt, CasavaOneEightSingleLanePerSampleDirFmt, BAMDirFmt, pd.DataFrame)

    """
    aligned_filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    unaligned_filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    bowtie_alignment = BAMDirFmt()
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    read_depths = []
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
            "--threads",
            str(threads),
        ]
        
        if very_sensitive:
            cmd.append("--very-sensitive")
        
        with tempfile.NamedTemporaryFile() as temp:
            result = subprocess.run(cmd, stdout=temp, stderr=subprocess.PIPE)
            lines = result.stderr.decode().split("\n")
            total_reads = extract_first_number(lines[0])
            unaligned_reads = extract_first_number(lines[2])
            aligned_one_time_reads = extract_first_number(lines[3])
            aligned_multiple_times_reads = extract_first_number(lines[4])
            read_depths.append(
                (sample_id, total_reads, unaligned_reads, aligned_one_time_reads, aligned_multiple_times_reads)
            )
            if not save_alignment:
                file_path = os.path.basename(sample_path)
                # Making an empty fastq file per sample
                with gzip.open(os.path.join(str(aligned_filtered_seqs), file_path), "w") as f1:
                    pass
                # Making an empty bam file per sample
                tmp_sam_file = "empty.sam"
                header_content = """@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:1\tLN:8
@PG\tID:bowtie2\tPN:bowtie2\tVN:2.3.4.1\tCL:"/usr/bin/bowtie2-align-s --wrapper basic-0 -x dummyref.fa -f dummyquery.fa"
A\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTTTTT\tIIIIIIII\tYT:Z:UU\n"""
                with open(tmp_sam_file, "w") as header_file:
                    header_file.write(header_content)
                with open(os.path.join(str(bowtie_alignment), f"{sample_id}.bam"), "w") as bowtie_file:
                    subprocess.run(
                        ["samtools", "view", "-bS", tmp_sam_file],
                        check=True,
                        stdout=bowtie_file,
                    )
                os.remove(tmp_sam_file)

            else:
                with open(os.path.join(str(bowtie_alignment), f"{sample_id}.bam"), "w") as bowtie_file:
                    subprocess.run(["samtools", "view", "-bS", temp.name], check=True, stdout=bowtie_file)
    reads_df = pd.DataFrame(
        read_depths,
        columns=["sample_id", "reads_total", "reads_unaligned", "reads_aligned_once", "reads_aligned_multiple"],
    )
    reads_df = reads_df.set_index("sample_id")
    return aligned_filtered_seqs, unaligned_filtered_seqs, bowtie_alignment, reads_df
