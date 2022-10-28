{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

{% if gpus %}
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --gres gpu:{{ gpus }}
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
