from q2_types.feature_table import FeatureTable
from qiime2.plugin import SemanticType

BowtieReadStatistics = SemanticType("BowtieReadStatistics", variant_of=FeatureTable.field["content"])
