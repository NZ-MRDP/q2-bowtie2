"""QIIME 2 plugin for bowtie2."""

import qiime2.plugin
from q2_types.bowtie2 import Bowtie2Index
from q2_types.feature_data import FeatureData, Sequence
from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from q2_types_genomics.per_sample_data._type import AlignmentMap

import q2_bowtie2

plugin = qiime2.plugin.Plugin(
    name="bowtie2",
    version="0.0.0",
    description="QIIME 2 plugin for bowtie2",
    website="https://huttenhower.sph.harvard.edu/humann/",
    package="q2_bowtie2",
    user_support_text=("bowtie around my neck is why they call me the gangster mac"),
    citation_text=None,
)

plugin.methods.register_function(
    function=q2_bowtie2.build,
    inputs={"reference_seqs": FeatureData[Sequence]},  # type: ignore
    parameters={},
    outputs=[("bowtie_index", Bowtie2Index)],  # type: ignore
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name="bowtie2 qiime plugin",
    description=(
        "bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing"
        " reads to long reference sequences. it is particularly good at aligning reads"
        " of about 50 up to 100s or 1,000s of characters, and particularly good at aligning"
        " to relatively long (e.g. mammalian) genomes."
    ),
)

plugin.methods.register_function(
    function=q2_bowtie2.align_single,
    inputs={
        "bowtie_database": Bowtie2Index,
        "demultiplexed_sequences": SampleData[SequencesWithQuality],  # type: ignore
    },
    parameters={},
    outputs=[
        ("aligned_reads", SampleData[SequencesWithQuality]),  # type: ignore
        ("unaligned_reads", SampleData[SequencesWithQuality]),  # type: ignore
        ("bowtie2_alignment", SampleData[AlignmentMap]),  # type: ignore
    ],  # type: ignore
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name="bowtie2 qiime plugin",
    description=("Description of bowtie2.align"),
)

plugin.methods.register_function(
    function=q2_bowtie2.align_paired,
    inputs={
        "bowtie_database": Bowtie2Index,
        "demultiplexed_sequences": SampleData[PairedEndSequencesWithQuality],  # type: ignore
    },
    parameters={},
    outputs=[
        ("aligned_reads", SampleData[PairedEndSequencesWithQuality]),  # type: ignore
        ("unaligned_reads", SampleData[PairedEndSequencesWithQuality]),  # type: ignore
        ("bowtie2_alignment", SampleData[AlignmentMap]),  # type: ignore
    ],
    input_descriptions={},
    parameter_descriptions={},
    output_descriptions={},
    name="bowtie2 qiime plugin",
    description=("Description of bowtie2.align"),
)
