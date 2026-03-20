"""QIIME 2 plugin for bowtie2."""

import importlib

import qiime2.plugin
from q2_types.bowtie2 import Bowtie2Index
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable
from q2_types.per_sample_sequences import (
    AlignmentMap,
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from q2_types_variant import GenBankSequence
from qiime2.plugin import Bool, Choices, Int, Str

import q2_bowtie2

from ._format import BowtieReadStatsDirFormat
from ._type import BowtieReadStatistics

plugin = qiime2.plugin.Plugin(
    name="bowtie2",
    version="0.1.4",
    description="QIIME 2 plugin for bowtie2",
    website="https://huttenhower.sph.harvard.edu/humann/",
    package="q2_bowtie2",
    user_support_text=("bowtie around my neck is why they call me the gangster mac"),
    citation_text=None,
)


plugin.methods.register_function(
    function=q2_bowtie2.build,
    inputs={"reference_seqs": FeatureData[Sequence | GenBankSequence]},  # type: ignore
    parameters={},
    outputs=[("bowtie_index", Bowtie2Index)],  # type: ignore
    input_descriptions={
        "reference_seqs": "Reference sequences to index for downstream Bowtie2 alignment.",
    },
    parameter_descriptions={},
    output_descriptions={
        "bowtie_index": "Bowtie2 FM-index built from the provided reference sequences.",
    },
    name="Build Bowtie2 index",
    description=(
        "Build a Bowtie2 index from reference sequences. Bowtie2 is an ultrafast and"
        " memory-efficient tool for aligning sequencing reads to long reference sequences."
        " It is particularly good at aligning reads of about 50 up to 100s or 1,000s of"
        " characters, and particularly good at aligning to relatively long (e.g. mammalian)"
        " genomes."
    ),
)


plugin.methods.register_function(
    function=q2_bowtie2.align_single,
    inputs={
        "bowtie_database": Bowtie2Index,
        "demultiplexed_sequences": SampleData[SequencesWithQuality],  # type: ignore
    },
    parameters={
        "threads": Int,
        "mode": Str % Choices(["end-to-end", "local"]),
        "sensitivity": Str % Choices(["very-fast", "fast", "sensitive", "very-sensitive"]),
        "save_alignment": Bool,
        "very_sensitive": Bool,
    },
    outputs=[
        ("aligned_reads", SampleData[SequencesWithQuality]),  # type: ignore
        ("unaligned_reads", SampleData[SequencesWithQuality]),  # type: ignore
        ("bowtie2_alignment", SampleData[AlignmentMap]),  # type: ignore
        ("read_features", FeatureTable[BowtieReadStatistics]),  # type: ignore
    ],  # type: ignore
    input_descriptions={
        "bowtie_database": "Bowtie2 index built from the target reference sequences.",
        "demultiplexed_sequences": "Single-end demultiplexed reads to align.",
    },
    parameter_descriptions={
        "threads": "Number of alignment threads to launch.",
        "mode": "Alignment mode. `end-to-end` aligns the full read; `local` allows soft clipping at the ends.",
        "sensitivity": "Bowtie2 preset controlling the speed/sensitivity tradeoff.",
        "save_alignment": "Whether to retain the BAM alignment output and aligned-read FASTQ records.",
        "very_sensitive": "Legacy convenience flag that overrides `sensitivity` with the `very-sensitive` preset.",
    },
    output_descriptions={
        "aligned_reads": "Reads reported by Bowtie2 as aligned to the reference.",
        "unaligned_reads": "Reads reported by Bowtie2 as unaligned.",
        "bowtie2_alignment": "Per-sample BAM alignment files emitted by Bowtie2.",
        "read_features": "Per-sample read alignment counts parsed from the Bowtie2 summary output.",
    },
    name="Align single-end reads with Bowtie2",
    description=(
        "Align single-end demultiplexed sequences to a Bowtie2 index. Returns aligned"
        " reads, unaligned reads, a BAM alignment file, and a table of per-sample read"
        " alignment statistics."
    ),
)

plugin.methods.register_function(
    function=q2_bowtie2.align_paired,
    inputs={
        "bowtie_database": Bowtie2Index,
        "demultiplexed_sequences": SampleData[PairedEndSequencesWithQuality],  # type: ignore
    },
    parameters={
        "threads": Int,
        "mode": Str % Choices(["end-to-end", "local"]),
        "sensitivity": Str % Choices(["very-fast", "fast", "sensitive", "very-sensitive"]),
        "very_sensitive": Bool,
    },
    outputs=[
        ("aligned_reads", SampleData[PairedEndSequencesWithQuality]),  # type: ignore
        ("unaligned_reads", SampleData[PairedEndSequencesWithQuality]),  # type: ignore
        ("bowtie2_alignment", SampleData[AlignmentMap]),  # type: ignore
    ],
    input_descriptions={
        "bowtie_database": "Bowtie2 index built from the target reference sequences.",
        "demultiplexed_sequences": "Paired-end demultiplexed reads to align.",
    },
    parameter_descriptions={
        "threads": "Number of alignment threads to launch.",
        "mode": "Alignment mode. `end-to-end` aligns the full read pair; `local` allows soft clipping.",
        "sensitivity": "Bowtie2 preset controlling the speed/sensitivity tradeoff.",
        "very_sensitive": "Legacy convenience flag that overrides `sensitivity` with the `very-sensitive` preset.",
    },
    output_descriptions={
        "aligned_reads": "Read pairs reported by Bowtie2 as aligned concordantly or discordantly to the reference.",
        "unaligned_reads": "Read pairs reported by Bowtie2 as unaligned.",
        "bowtie2_alignment": "Per-sample BAM alignment files emitted by Bowtie2.",
    },
    name="Align paired-end reads with Bowtie2",
    description=(
        "Align paired-end demultiplexed sequences to a Bowtie2 index. Returns aligned"
        " reads, unaligned reads, and a BAM alignment file."
    ),
)

plugin.register_formats(BowtieReadStatsDirFormat)
plugin.register_semantic_type_to_format(
    FeatureTable[BowtieReadStatistics], artifact_format=BowtieReadStatsDirFormat  # type: ignore
)
importlib.import_module("q2_bowtie2._transformer")
