{% extends "slurm.sh" %}

{% block header %}
{{- super () -}}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
{% set cpus = operations|map(attribute='directives.np')|sum %}

{% if gpus %}
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --nodelist=potoff3x
{%- else %}
#SBATCH -N 1
#SBATCH --mail-type=ALL
{%- endif %}


echo  "Running on host" hostname
echo  "Time is" date

conda activate mosdef-study38

{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}
