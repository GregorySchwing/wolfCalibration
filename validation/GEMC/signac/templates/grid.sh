{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

{% if gpus %}
#SBATCH -q gpu
#SBATCH --gres gpu:{{ gpus }}
{%- else %}
#SBATCH -q primary
{%- endif %}

#SBATCH -N 1


echo  "Running on host" hostname
echo  "Time is" date
conda activate nobel_gas

module load python/3.8

{% if gpus %}
module load cuda/10.0
nvidia-smi
{%- endif %}

{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}
