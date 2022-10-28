{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

#SBATCH -N 1
#SBATCH --mail-type=ALL

echo  "Running on host" hostname
echo  "Time is" date

conda activate mosdef-study38

{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}
