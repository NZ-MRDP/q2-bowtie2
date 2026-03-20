import qiime2.plugin.model as model


class BowtieReadStatsFileFormat(model.TextFileFormat):
    """TSV file format for Bowtie2 per-sample read alignment statistics."""

    # TODO: Add validation
    def _validate_(self, *args):
        pass


BowtieReadStatsDirFormat = model.SingleFileDirectoryFormat(
    "BowtieReadStatsFileFormat", "bowtie2_read_stats.tsv", BowtieReadStatsFileFormat
)
