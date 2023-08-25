import qiime2.plugin.model as model


class BowtieReadStatsFileFormat(model.TextFileFormat):
    """BowtieReadStatsFileFormat."""

    # TODO: Add validation
    def _validate_(self, *args):
        pass


BowtieReadStatsDirFormat = model.SingleFileDirectoryFormat(
    "BowtieReadStatsFileFormat", "bowtie2_read_stats.tsv", BowtieReadStatsFileFormat
)
