from cana.boolean_node import BooleanNode
from binarize import to_binary
from display_tables import plot_look_up_table, plot_schemata

rules = [2, 3, 4, 6]

# rule 7
rule = 7
digits = 4
arr2 = to_binary(rule, digits)
bn = BooleanNode.from_output_list(arr2)
bn.input_symmetry()

plot_look_up_table(rule, bn, snakemake.output.redlut)
plot_schemata(rule, bn, snakemake.output.redwc)

# rule 3
rule = 3
digits = 4
arr2 = to_binary(rule, digits)
bn = BooleanNode.from_output_list(arr2)
bn.input_symmetry()

plot_look_up_table(rule, bn, snakemake.output.unqlut)
plot_schemata(rule, bn, snakemake.output.unqwc)

# rule 1
rule = 1
digits = 4
arr2 = to_binary(rule, digits)
bn = BooleanNode.from_output_list(arr2)
bn.input_symmetry()

plot_look_up_table(rule, bn, snakemake.output.alllut)
plot_schemata(rule, bn, snakemake.output.allwc)

# rule 6
rule = 6
digits = 4
arr2 = to_binary(rule, digits)
bn = BooleanNode.from_output_list(arr2)
bn.input_symmetry()

plot_look_up_table(rule, bn, snakemake.output.synlut)
plot_schemata(rule, bn, snakemake.output.synwc)
